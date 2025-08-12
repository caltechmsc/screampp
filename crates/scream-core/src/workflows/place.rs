use crate::core::forcefield::parameterization::Parameterizer;
use crate::core::forcefield::params::Forcefield;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::core::rotamers::library::RotamerLibrary;
use crate::core::topology::registry::TopologyRegistry;
use crate::engine::cache::ELCache;
use crate::engine::config::PlacementConfig;
use crate::engine::context::{OptimizationContext, resolve_selection_to_ids};
use crate::engine::error::EngineError;
use crate::engine::placement::place_rotamer_on_system;
use crate::engine::progress::{Progress, ProgressReporter};
use crate::engine::state::{
    InitialState, OptimizationState, PlacementResult, Solution, SolutionState,
};
use crate::engine::tasks;
use rand::prelude::IteratorRandom;
use rand::{Rng, seq::SliceRandom, thread_rng};
use std::collections::{HashMap, HashSet};
use tracing::{info, instrument};

const CLASH_THRESHOLD_KCAL_MOL: f64 = 25.0;

#[instrument(skip_all, name = "placement_workflow")]
pub fn run(
    initial_system: &MolecularSystem,
    config: &PlacementConfig,
    reporter: &ProgressReporter,
) -> Result<PlacementResult, EngineError> {
    // === Phase 0: Preparation and Parameterization ===
    reporter.report(Progress::PhaseStart {
        name: "Preparation",
    });
    info!("Starting workflow setup: loading resources and parameterizing system.");

    let forcefield = Forcefield::load(
        &config.forcefield.forcefield_path,
        &config.forcefield.delta_params_path,
    )?;
    let topology_registry = TopologyRegistry::load(&config.topology_registry_path)?;
    let mut rotamer_library = RotamerLibrary::load(
        &config.sampling.rotamer_library_path,
        &topology_registry,
        &forcefield,
        config.forcefield.s_factor,
    )?;

    let active_residues = prepare_context(
        initial_system,
        config,
        &mut rotamer_library,
        &topology_registry,
    )?;

    let parameterizer =
        Parameterizer::new(&forcefield, &topology_registry, config.forcefield.s_factor);
    let mut working_system = initial_system.clone();
    parameterizer.parameterize_system(&mut working_system)?;

    let context = OptimizationContext::new(
        &working_system,
        &forcefield,
        reporter,
        config,
        &rotamer_library,
        &topology_registry,
    );

    reporter.report(Progress::PhaseFinish);

    // === Phase 1: Calculate constant energy and initial state ===
    let (initial_state, energy_offset_constant) =
        calculate_initial_state(&context, &active_residues)?;

    // === Phase 2: Precompute empty lattice energy (EL Energy) ===
    let el_cache = tasks::el_energy::run(&context)?;

    // === Phase 3: Initialize optimization state (place ground state) ===
    let mut state = initialize_optimization_state(&context, &active_residues, &el_cache)?;

    // === Phase 4: Clash resolution ===
    resolve_clashes(&mut state, &active_residues, &context, &el_cache)?;

    // === Phase 5: Simulated annealing (optional) ===
    if config.optimization.simulated_annealing.is_some() {
        run_simulated_annealing(&mut state, &active_residues, &context, &el_cache)?;
    }

    // === Phase 6: Final refinement ===
    final_refinement(&mut state, &active_residues, &context, &el_cache)?;

    // === Phase 7: Organize and return results ===
    let result = finalize_results(
        state,
        initial_state,
        energy_offset_constant,
        config.optimization.num_solutions,
    );

    info!(
        "Workflow complete. Returning {} solution(s).",
        result.solutions.len()
    );
    Ok(result)
}

fn prepare_context(
    initial_system: &MolecularSystem,
    config: &PlacementConfig,
    rotamer_library: &mut RotamerLibrary,
    topology_registry: &TopologyRegistry,
) -> Result<HashSet<ResidueId>, EngineError> {
    let active_residues = resolve_selection_to_ids(
        initial_system,
        &config.residues_to_optimize,
        rotamer_library,
    )?;

    if config.optimization.include_input_conformation {
        info!("Including original side-chain conformations in the rotamer library.");
        rotamer_library.include_system_conformations(
            initial_system,
            &active_residues,
            topology_registry,
        );
    }
    Ok(active_residues)
}

fn calculate_initial_state(
    context: &OptimizationContext<PlacementConfig>,
    active_residues: &HashSet<ResidueId>,
) -> Result<(InitialState, EnergyTerm), EngineError> {
    context.reporter.report(Progress::PhaseStart {
        name: "Calculating Initial State",
    });
    info!("Calculating energy of the initial input conformation.");

    let energy_offset_constant = tasks::fixed_energy::run(context)?;
    let initial_interaction_energy =
        tasks::interaction_energy::run(context.system, context.forcefield, active_residues)?;
    let initial_el_energy = tasks::el_energy::calculate_current(context)?;

    let initial_optimization_score_term = initial_interaction_energy + initial_el_energy;
    let initial_optimization_score = initial_optimization_score_term.total();
    let initial_total_energy = initial_optimization_score + energy_offset_constant.total();

    let initial_state = InitialState {
        system: context.system.clone(),
        total_energy: initial_total_energy,
        optimization_score: initial_optimization_score,
    };

    info!(
        total_energy = initial_total_energy,
        optimization_score = initial_optimization_score,
        "Initial state calculated."
    );
    context.reporter.report(Progress::PhaseFinish);
    Ok((initial_state, energy_offset_constant))
}

fn initialize_optimization_state(
    context: &OptimizationContext<PlacementConfig>,
    active_residues: &HashSet<ResidueId>,
    el_cache: &ELCache,
) -> Result<OptimizationState, EngineError> {
    context.reporter.report(Progress::PhaseStart {
        name: "Initializing Ground State",
    });
    info!("Placing ground-state rotamers to initialize optimization.");

    let mut ground_state_system = context.system.clone();
    let mut ground_state_rotamers = HashMap::new();
    let mut ground_state_el_energy = EnergyTerm::default();

    for &residue_id in active_residues {
        let residue = context.system.residue(residue_id).unwrap();
        if let Some(residue_type) = residue.residue_type {
            if let Some((idx, energy)) = el_cache.find_ground_state_for(residue_id, residue_type) {
                ground_state_rotamers.insert(residue_id, idx);
                ground_state_el_energy += *energy;

                let rotamer = &context
                    .rotamer_library
                    .get_rotamers_for(residue_type)
                    .unwrap()[idx];
                let res_name = residue_type.to_three_letter();
                let topology = context.topology_registry.get(res_name).unwrap();
                place_rotamer_on_system(&mut ground_state_system, residue_id, rotamer, topology)?;
            }
        }
    }

    let ground_state_interaction =
        tasks::interaction_energy::run(&ground_state_system, context.forcefield, active_residues)?;

    let ground_state_optimization_score =
        (ground_state_el_energy + ground_state_interaction).total();

    info!(
        score = ground_state_optimization_score,
        "Ground state optimization score calculated."
    );
    context.reporter.report(Progress::PhaseFinish);

    Ok(OptimizationState::new(
        ground_state_system,
        ground_state_rotamers,
        ground_state_optimization_score,
        context.config.optimization.num_solutions,
    ))
}

fn resolve_clashes(
    state: &mut OptimizationState,
    active_residues: &HashSet<ResidueId>,
    context: &OptimizationContext<'a, PlacementConfig>,
    el_cache: &ELCache,
) -> Result<(), EngineError> {
    context.reporter.report(Progress::PhaseStart {
        name: "Clash Resolution",
    });
    info!("Starting iterative clash resolution loop.");

    let max_iterations = context.config.optimization.max_iterations;
    let convergence_energy_threshold = context.config.optimization.convergence.energy_threshold;
    let convergence_patience_iterations =
        context.config.optimization.convergence.patience_iterations;

    let mut last_total_energy = state.best_solution().map(|s| s.energy).unwrap_or(f64::MAX);
    let mut iterations_without_significant_improvement = 0;

    for iteration in 0..max_iterations {
        let clashes = tasks::clash_detection::run(
            &state.working_state.system,
            context.forcefield,
            active_residues,
            CLASH_THRESHOLD,
            context.reporter,
        )?;

        context.reporter.report(Progress::StatusUpdate {
            text: format!(
                "Pass {}/{}, Clashes Found: {}",
                iteration + 1,
                max_iterations,
                clashes.len()
            ),
        });

        if clashes.is_empty() {
            info!(
                "System converged after {} clash resolution iterations (no clashes found).",
                iteration
            );
            context.reporter.report(Progress::Message(format!(
                "Converged after {} iterations (no clashes).",
                iteration
            )));
            break;
        }

        let worst_clash = &clashes[0];

        let doublet_result = tasks::doublet_optimization::run(
            worst_clash.residue_a,
            worst_clash.residue_b,
            &state.working_state.system,
            el_cache,
            context,
        )?;

        let res_a_id = worst_clash.residue_a;
        let res_b_id = worst_clash.residue_b;
        let res_a_type = state
            .working_state
            .system
            .residue(res_a_id)
            .unwrap()
            .residue_type
            .unwrap();
        let rotamer_a = &context
            .rotamer_library
            .get_rotamers_for(res_a_type)
            .unwrap()[doublet_result.rotamer_idx_a];

        let res_name_a = res_a_type.to_three_letter();
        let topology_a = context.topology_registry.get(res_name_a).ok_or_else(|| {
            EngineError::TopologyNotFound {
                residue_name: res_name_a.to_string(),
            }
        })?;
        place_rotamer_on_system(
            &mut state.working_state.system,
            res_a_id,
            rotamer_a,
            topology_a,
        )?;

        state
            .working_state
            .rotamers
            .insert(res_a_id, doublet_result.rotamer_idx_a);

        let res_b_type = state
            .working_state
            .system
            .residue(res_b_id)
            .unwrap()
            .residue_type
            .unwrap();
        let rotamer_b = &context
            .rotamer_library
            .get_rotamers_for(res_b_type)
            .unwrap()[doublet_result.rotamer_idx_b];

        let res_name_b = res_b_type.to_three_letter();
        let topology_b = context.topology_registry.get(res_name_b).ok_or_else(|| {
            EngineError::TopologyNotFound {
                residue_name: res_name_b.to_string(),
            }
        })?;
        place_rotamer_on_system(
            &mut state.working_state.system,
            res_b_id,
            rotamer_b,
            topology_b,
        )?;

        state
            .working_state
            .rotamers
            .insert(res_b_id, doublet_result.rotamer_idx_b);

        let energy_after_update = tasks::total_energy::run(
            &state.working_state.system,
            context.forcefield,
            active_residues,
            &state.working_state.rotamers,
            el_cache,
        )?
        .total();
        state.current_energy = energy_after_update;
        state.submit_current_solution();

        let current_best_energy = state.best_solution().unwrap().energy;
        let energy_improvement = last_total_energy - current_best_energy;

        if energy_improvement < convergence_energy_threshold {
            iterations_without_significant_improvement += 1;
        } else {
            iterations_without_significant_improvement = 0;
        }

        last_total_energy = current_best_energy;

        if iterations_without_significant_improvement >= convergence_patience_iterations {
            info!(
                "System converged after {} clash resolution iterations (energy stabilized below {}).",
                iteration, convergence_energy_threshold
            );
            context.reporter.report(Progress::Message(format!(
                "Converged after {} iterations (energy stabilized).",
                iteration
            )));
            break;
        }
    }
    context.reporter.report(Progress::PhaseFinish);
    Ok(())
}

#[instrument(skip_all, name = "final_refinement_loop")]
fn final_refinement<'a>(
    state: &mut OptimizationState,
    active_residues: &HashSet<ResidueId>,
    context: &OptimizationContext<'a, PlacementConfig>,
    el_cache: &ELCache,
) -> Result<(), EngineError> {
    context.reporter.report(Progress::PhaseStart {
        name: "Final Refinement",
    });
    info!("Starting final refinement (singlet optimization).");

    let final_refinement_iterations = context.config.optimization.final_refinement_iterations;
    if final_refinement_iterations == 0 {
        info!("Final refinement skipped as configured (0 iterations).");
        context.reporter.report(Progress::PhaseFinish);
        return Ok(());
    }

    for i in 0..final_refinement_iterations {
        let mut changed_in_cycle = false;

        info!(
            "Beginning refinement pass {}/{}.",
            i + 1,
            final_refinement_iterations
        );

        context.reporter.report(Progress::StatusUpdate {
            text: format!("Pass {}/{}", i + 1, final_refinement_iterations),
        });

        context.reporter.report(Progress::TaskStart {
            total: active_residues.len() as u64,
        });

        let mut residues_to_process: Vec<ResidueId> = active_residues.iter().cloned().collect();
        residues_to_process.shuffle(&mut thread_rng());

        for (res_idx, &residue_id) in residues_to_process.iter().enumerate() {
            debug!(
                "Refining residue {}/{} (ID: {:?}, Type: {})",
                res_idx + 1,
                residues_to_process.len(),
                residue_id,
                state.working_state.system.residue(residue_id).unwrap().name
            );

            let residue_type = state
                .working_state
                .system
                .residue(residue_id)
                .unwrap()
                .residue_type
                .unwrap();
            let rotamers = context
                .rotamer_library
                .get_rotamers_for(residue_type)
                .unwrap();

            let residue_name = residue_type.to_three_letter();
            let topology = context.topology_registry.get(residue_name).ok_or_else(|| {
                EngineError::TopologyNotFound {
                    residue_name: residue_name.to_string(),
                }
            })?;

            let current_rot_idx = *state.working_state.rotamers.get(&residue_id).unwrap();
            let mut best_rotamer_idx = current_rot_idx;
            let mut best_energy = state.current_energy;

            let mut temp_system_for_eval = state.working_state.system.clone();
            let mut temp_rotamers_for_eval = state.working_state.rotamers.clone();
            let energy_of_current_rotamer = state.current_energy;

            for (idx, rotamer) in rotamers.iter().enumerate() {
                if idx == current_rot_idx {
                    continue;
                }

                place_rotamer_on_system(&mut temp_system_for_eval, residue_id, rotamer, topology)?;
                temp_rotamers_for_eval.insert(residue_id, idx);

                let new_energy = tasks::total_energy::run(
                    &temp_system_for_eval,
                    context.forcefield,
                    active_residues,
                    &temp_rotamers_for_eval,
                    el_cache,
                )?
                .total();

                if new_energy < best_energy {
                    best_energy = new_energy;
                    best_rotamer_idx = idx;
                }

                place_rotamer_on_system(
                    &mut temp_system_for_eval,
                    residue_id,
                    &rotamers[current_rot_idx],
                    topology,
                )?;
                temp_rotamers_for_eval.insert(residue_id, current_rot_idx);
            }

            if best_rotamer_idx != current_rot_idx {
                changed_in_cycle = true;

                debug!(
                    "  Found better rotamer for {:?}. Index: {} -> {}. Global energy: {:.2} -> {:.2}",
                    residue_id,
                    current_rot_idx,
                    best_rotamer_idx,
                    energy_of_current_rotamer,
                    best_energy
                );

                let best_rotamer = &rotamers[best_rotamer_idx];
                place_rotamer_on_system(
                    &mut state.working_state.system,
                    residue_id,
                    best_rotamer,
                    topology,
                )?;
                state
                    .working_state
                    .rotamers
                    .insert(residue_id, best_rotamer_idx);
                state.current_energy = best_energy;
                state.submit_current_solution();
            }

            context
                .reporter
                .report(Progress::TaskIncrement { amount: 1 });
        }

        context.reporter.report(Progress::TaskFinish);

        if changed_in_cycle {
            info!(
                "Refinement pass {} complete. At least one rotamer was changed. Current best energy: {:.4}",
                i + 1,
                state.best_solution().map(|s| s.energy).unwrap_or(f64::NAN)
            );
        } else {
            info!(
                "Refinement converged after {} pass(es). No further improvements in this pass.",
                i + 1
            );
            break;
        }
    }

    state.submit_current_solution();

    context.reporter.report(Progress::PhaseFinish);
    Ok(())
}

#[instrument(skip_all, name = "simulated_annealing_loop")]
fn run_simulated_annealing<'a>(
    state: &mut OptimizationState,
    active_residues: &HashSet<ResidueId>,
    context: &OptimizationContext<'a, PlacementConfig>,
    el_cache: &ELCache,
) -> Result<(), EngineError> {
    context.reporter.report(Progress::PhaseStart {
        name: "Simulated Annealing",
    });
    info!("Starting Simulated Annealing exploration.");

    let sa_config = context
        .config
        .optimization
        .simulated_annealing
        .as_ref()
        .unwrap();

    let initial_temperature = sa_config.initial_temperature;
    let final_temperature = sa_config.final_temperature;
    let cooling_rate = sa_config.cooling_rate;
    let steps_per_temperature = sa_config.steps_per_temperature;

    let mut rng = thread_rng();
    let mut current_temperature = initial_temperature;

    let active_residues_vec: Vec<ResidueId> = active_residues.iter().cloned().collect();

    let mut temp_step_current = 0;
    while current_temperature > final_temperature {
        temp_step_current += 1;
        context.reporter.report(Progress::StatusUpdate {
            text: format!(
                "SA Temp: {:.4} (Step {})",
                current_temperature, temp_step_current
            ),
        });

        context.reporter.report(Progress::TaskStart {
            total: steps_per_temperature as u64,
        });

        for _step in 0..steps_per_temperature {
            let residue_to_perturb_id = *active_residues_vec.choose(&mut rng).unwrap();

            let residue_type = state
                .working_state
                .system
                .residue(residue_to_perturb_id)
                .unwrap()
                .residue_type
                .unwrap();
            let rotamers = context
                .rotamer_library
                .get_rotamers_for(residue_type)
                .unwrap();

            let residue_name = residue_type.to_three_letter();
            let topology = context.topology_registry.get(residue_name).ok_or_else(|| {
                EngineError::TopologyNotFound {
                    residue_name: residue_name.to_string(),
                }
            })?;

            if rotamers.len() <= 1 {
                context
                    .reporter
                    .report(Progress::TaskIncrement { amount: 1 });
                continue;
            }

            let current_rot_idx = *state
                .working_state
                .rotamers
                .get(&residue_to_perturb_id)
                .unwrap();
            let new_rot_idx = loop {
                let r_idx = rng.gen_range(0..rotamers.len());
                if r_idx != current_rot_idx {
                    break r_idx;
                }
            };
            let new_rotamer = &rotamers[new_rot_idx];

            let current_total_energy = state.current_energy;

            let mut temp_system_for_eval = state.working_state.system.clone();
            let mut temp_rotamers_for_eval = state.working_state.rotamers.clone();

            place_rotamer_on_system(
                &mut temp_system_for_eval,
                residue_to_perturb_id,
                new_rotamer,
                topology,
            )?;
            temp_rotamers_for_eval.insert(residue_to_perturb_id, new_rot_idx);

            let energy_after_change = tasks::total_energy::run(
                &temp_system_for_eval,
                context.forcefield,
                active_residues,
                &temp_rotamers_for_eval,
                el_cache,
            )?
            .total();

            let delta_e = energy_after_change - current_total_energy;

            if delta_e < 0.0 {
                info!("  SA accepted: ΔE={:.4} (downhill)", delta_e);
                state.working_state.system = temp_system_for_eval;
                state.working_state.rotamers = temp_rotamers_for_eval;
                state.current_energy = energy_after_change;
                state.submit_current_solution();
            } else {
                let acceptance_probability = (-delta_e / current_temperature).exp();
                if rng.r#gen::<f64>() < acceptance_probability {
                    info!(
                        "  SA accepted: ΔE={:.4} (uphill, prob={:.4})",
                        delta_e, acceptance_probability
                    );
                    state.working_state.system = temp_system_for_eval;
                    state.working_state.rotamers = temp_rotamers_for_eval;
                    state.current_energy = energy_after_change;
                    state.submit_current_solution();
                } else {
                    info!(
                        "  SA rejected: ΔE={:.4} (uphill, prob={:.4})",
                        delta_e, acceptance_probability
                    );
                }
            }
            context
                .reporter
                .report(Progress::TaskIncrement { amount: 1 });
        }

        context.reporter.report(Progress::TaskFinish);

        current_temperature *= cooling_rate;
    }
    info!("Simulated Annealing finished. Final temperature reached.");
    context.reporter.report(Progress::PhaseFinish);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::parameterization::Parameterizer;
    use crate::core::forcefield::params::Forcefield;
    use crate::core::forcefield::scoring::Scorer;
    use crate::core::models::atom::Atom;
    use crate::core::models::chain::ChainType;
    use crate::core::models::ids::{AtomId, ChainId};
    use crate::core::models::residue::ResidueType;
    use crate::core::topology::registry::TopologyRegistry;
    use crate::engine::config::{
        ConvergenceConfig, PlacementConfigBuilder, ResidueSelection, ResidueSpecifier,
        SimulatedAnnealingConfig,
    };
    use nalgebra::{Point3, Rotation3, Vector3};
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::Write;
    use std::path::{Path, PathBuf};
    use tempfile::{TempDir, tempdir};

    mod setup {
        use super::*;
        use crate::core::models::topology::BondOrder;
        use phf::{Map, phf_map};

        struct TemplateAtom {
            name: &'static str,
            ff_type: &'static str,
            charge: f64,
            pos: Point3<f64>,
        }

        #[derive(Clone, Copy)]
        struct ResidueTemplate {
            atoms: &'static [TemplateAtom],
            bonds: &'static [(&'static str, &'static str)],
            n_anchor: &'static str,
            ca_anchor: &'static str,
            c_anchor: &'static str,
        }

        static ALA_TEMPLATE: ResidueTemplate = ResidueTemplate {
            atoms: &[
                TemplateAtom {
                    name: "N",
                    ff_type: "N_R",
                    charge: -0.47,
                    pos: Point3::new(-0.52, 1.36, 0.0),
                },
                TemplateAtom {
                    name: "H",
                    ff_type: "H_",
                    charge: 0.31,
                    pos: Point3::new(-1.31, 1.98, 0.0),
                },
                TemplateAtom {
                    name: "CA",
                    ff_type: "C_31",
                    charge: 0.07,
                    pos: Point3::new(0.0, 0.0, 0.0),
                },
                TemplateAtom {
                    name: "HA",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(0.63, -0.47, -0.76),
                },
                TemplateAtom {
                    name: "CB",
                    ff_type: "C_33",
                    charge: -0.27,
                    pos: Point3::new(-0.76, -0.8, -1.08),
                },
                TemplateAtom {
                    name: "HB1",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(-0.21, -1.74, -1.16),
                },
                TemplateAtom {
                    name: "HB2",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(-1.6, -0.4, -1.67),
                },
                TemplateAtom {
                    name: "HB3",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(-1.2, -1.2, -0.19),
                },
                TemplateAtom {
                    name: "C",
                    ff_type: "C_R",
                    charge: 0.51,
                    pos: Point3::new(1.2, -0.1, 0.9),
                },
                TemplateAtom {
                    name: "O",
                    ff_type: "O_2",
                    charge: -0.51,
                    pos: Point3::new(1.8, 0.7, 1.4),
                },
            ],
            bonds: &[
                ("N", "H"),
                ("N", "CA"),
                ("CA", "HA"),
                ("CA", "CB"),
                ("CA", "C"),
                ("CB", "HB1"),
                ("CB", "HB2"),
                ("CB", "HB3"),
                ("C", "O"),
            ],
            n_anchor: "N",
            ca_anchor: "CA",
            c_anchor: "C",
        };

        static GLY_TEMPLATE: ResidueTemplate = ResidueTemplate {
            atoms: &[
                TemplateAtom {
                    name: "N",
                    ff_type: "N_R",
                    charge: -0.47,
                    pos: Point3::new(-0.52, 1.36, 0.0),
                },
                TemplateAtom {
                    name: "H",
                    ff_type: "H_",
                    charge: 0.31,
                    pos: Point3::new(-1.31, 1.98, 0.0),
                },
                TemplateAtom {
                    name: "CA",
                    ff_type: "C_32",
                    charge: -0.02,
                    pos: Point3::new(0.0, 0.0, 0.0),
                },
                TemplateAtom {
                    name: "HA1",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(0.63, -0.47, -0.76),
                },
                TemplateAtom {
                    name: "HA2",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(-0.63, -0.47, 0.76),
                },
                TemplateAtom {
                    name: "C",
                    ff_type: "C_R",
                    charge: 0.51,
                    pos: Point3::new(1.2, -0.1, 0.9),
                },
                TemplateAtom {
                    name: "O",
                    ff_type: "O_2",
                    charge: -0.51,
                    pos: Point3::new(1.8, 0.7, 1.4),
                },
            ],
            bonds: &[
                ("N", "H"),
                ("N", "CA"),
                ("CA", "HA1"),
                ("CA", "HA2"),
                ("CA", "C"),
                ("C", "O"),
            ],
            n_anchor: "N",
            ca_anchor: "CA",
            c_anchor: "C",
        };

        static LEU_TEMPLATE: ResidueTemplate = ResidueTemplate {
            atoms: &[
                TemplateAtom {
                    name: "N",
                    ff_type: "N_R",
                    charge: -0.47,
                    pos: Point3::new(-0.52, 1.36, 0.0),
                },
                TemplateAtom {
                    name: "CA",
                    ff_type: "C_31",
                    charge: 0.07,
                    pos: Point3::new(0.0, 0.0, 0.0),
                },
                TemplateAtom {
                    name: "C",
                    ff_type: "C_R",
                    charge: 0.51,
                    pos: Point3::new(1.2, -0.1, 0.9),
                },
                TemplateAtom {
                    name: "O",
                    ff_type: "O_2",
                    charge: -0.51,
                    pos: Point3::new(1.8, 0.7, 1.4),
                },
                TemplateAtom {
                    name: "H",
                    ff_type: "H_",
                    charge: 0.31,
                    pos: Point3::new(-1.31, 1.98, 0.0),
                },
                TemplateAtom {
                    name: "HA",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(0.6, -0.5, -0.7),
                },
                TemplateAtom {
                    name: "CB",
                    ff_type: "C_32",
                    charge: -0.18,
                    pos: Point3::new(-0.8, -0.8, -1.1),
                },
                TemplateAtom {
                    name: "HB1",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(-0.2, -1.7, -1.2),
                },
                TemplateAtom {
                    name: "HB2",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(-1.6, -0.4, -1.7),
                },
                TemplateAtom {
                    name: "CG",
                    ff_type: "C_31",
                    charge: -0.09,
                    pos: Point3::new(-1.5, -1.5, 0.0),
                },
                TemplateAtom {
                    name: "HG",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(-1.0, -2.5, 0.0),
                },
                TemplateAtom {
                    name: "CD1",
                    ff_type: "C_33",
                    charge: -0.27,
                    pos: Point3::new(-2.5, -1.0, 0.8),
                },
                TemplateAtom {
                    name: "HD11",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(-2.2, -0.5, 1.7),
                },
                TemplateAtom {
                    name: "HD12",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(-3.0, -1.8, 1.2),
                },
                TemplateAtom {
                    name: "HD13",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(-3.3, -0.5, 0.2),
                },
                TemplateAtom {
                    name: "CD2",
                    ff_type: "C_33",
                    charge: -0.27,
                    pos: Point3::new(-2.0, -2.5, -1.0),
                },
                TemplateAtom {
                    name: "HD21",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(-1.5, -3.3, -1.5),
                },
                TemplateAtom {
                    name: "HD22",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(-2.8, -2.9, -0.5),
                },
                TemplateAtom {
                    name: "HD23",
                    ff_type: "H_",
                    charge: 0.09,
                    pos: Point3::new(-2.5, -2.0, -1.8),
                },
            ],
            bonds: &[
                ("N", "H"),
                ("N", "CA"),
                ("CA", "HA"),
                ("CA", "C"),
                ("C", "O"),
                ("CA", "CB"),
                ("CB", "HB1"),
                ("CB", "HB2"),
                ("CB", "CG"),
                ("CG", "HG"),
                ("CG", "CD1"),
                ("CG", "CD2"),
                ("CD1", "HD11"),
                ("CD1", "HD12"),
                ("CD1", "HD13"),
                ("CD2", "HD21"),
                ("CD2", "HD22"),
                ("CD2", "HD23"),
            ],
            n_anchor: "N",
            ca_anchor: "CA",
            c_anchor: "C",
        };

        static RESIDUE_TEMPLATES: Map<&'static str, ResidueTemplate> = phf_map! {
            "ALA" => ALA_TEMPLATE,
            "GLY" => GLY_TEMPLATE,
            "LEU" => LEU_TEMPLATE,
        };

        pub struct TestSystemBuilder {
            system: MolecularSystem,
            last_c_id: Option<AtomId>,
            chain_id: ChainId,
        }

        impl TestSystemBuilder {
            pub fn new() -> Self {
                let mut system = MolecularSystem::new();
                let chain_id = system.add_chain('A', ChainType::Protein);
                Self {
                    system,
                    last_c_id: None,
                    chain_id,
                }
            }

            pub fn add_residue(mut self, residue_type: ResidueType, residue_number: isize) -> Self {
                let template = RESIDUE_TEMPLATES
                    .get(residue_type.to_three_letter())
                    .expect("Residue template not found");

                let (rotation, translation) = if let Some(prev_c_id) = self.last_c_id {
                    let prev_c_atom = self.system.atom(prev_c_id).unwrap();
                    let prev_residue = self.system.residue(prev_c_atom.residue_id).unwrap();

                    let prev_ca_id = prev_residue.get_first_atom_id_by_name("CA").unwrap();
                    let prev_ca_atom = self.system.atom(prev_ca_id).unwrap();

                    let target_c = prev_c_atom.position;
                    let target_ca = prev_ca_atom.position;
                    let target_axis = (target_c - target_ca).normalize();
                    let target_origin = target_c + target_axis * 1.33;

                    let template_n_pos = template
                        .atoms
                        .iter()
                        .find(|a| a.name == template.n_anchor)
                        .unwrap()
                        .pos;
                    let template_ca_pos = template
                        .atoms
                        .iter()
                        .find(|a| a.name == template.ca_anchor)
                        .unwrap()
                        .pos;
                    let source_axis = (template_n_pos - template_ca_pos).normalize();

                    let rotation = Rotation3::rotation_between(&source_axis, &-target_axis)
                        .unwrap_or_else(Rotation3::identity);
                    let translation = target_origin - rotation * template_n_pos;
                    (rotation, translation)
                } else {
                    (Rotation3::identity(), Vector3::zeros())
                };

                let res_id = self
                    .system
                    .add_residue(
                        self.chain_id,
                        residue_number,
                        residue_type.to_three_letter(),
                        Some(residue_type),
                    )
                    .unwrap();

                let mut name_to_id = HashMap::new();
                for atom_template in template.atoms {
                    let new_pos = rotation * atom_template.pos + translation;
                    let mut atom = Atom::new(atom_template.name, res_id, new_pos);
                    atom.force_field_type = atom_template.ff_type.to_string();
                    atom.partial_charge = atom_template.charge;

                    let new_id = self.system.add_atom_to_residue(res_id, atom).unwrap();
                    name_to_id.insert(atom_template.name, new_id);
                }

                for &(name1, name2) in template.bonds {
                    self.system
                        .add_bond(name_to_id[name1], name_to_id[name2], BondOrder::Single)
                        .unwrap();
                }

                if let Some(prev_c_id) = self.last_c_id {
                    self.system
                        .add_bond(prev_c_id, name_to_id[template.n_anchor], BondOrder::Single)
                        .unwrap();
                }

                self.last_c_id = Some(name_to_id[template.c_anchor]);
                self
            }

            pub fn build(self) -> MolecularSystem {
                self.system
            }
        }

        pub struct TestEnvironment {
            pub temp_dir: TempDir,
            pub forcefield: Forcefield,
            pub topology_registry: TopologyRegistry,
            pub rotamer_lib_path: PathBuf,
            pub initial_system: MolecularSystem,
            forcefield_path: PathBuf,
            delta_path: PathBuf,
            topology_reg_path: PathBuf,
        }

        impl TestEnvironment {
            pub fn new() -> Self {
                let temp_dir = tempdir().unwrap();

                let forcefield_path = create_dummy_forcefield_file(temp_dir.path());
                let delta_path = create_dummy_delta_file(temp_dir.path());
                let topology_reg_path = create_dummy_topology_reg_file(temp_dir.path());
                let rotamer_lib_path = create_dummy_rotamer_lib_file(temp_dir.path());

                let forcefield = Forcefield::load(&forcefield_path, &delta_path)
                    .expect("Failed to load dummy forcefield");
                let topology_registry = TopologyRegistry::load(&topology_reg_path)
                    .expect("Failed to load dummy topology");

                let mut system = TestSystemBuilder::new()
                    .add_residue(ResidueType::Alanine, 1)
                    .add_residue(ResidueType::Glycine, 2)
                    .add_residue(ResidueType::Leucine, 3)
                    .build();

                let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 0.0);
                parameterizer
                    .parameterize_system(&mut system)
                    .expect("Parameterization failed in test setup");

                Self {
                    temp_dir,
                    forcefield,
                    topology_registry,
                    rotamer_lib_path,
                    initial_system: system,
                    forcefield_path,
                    delta_path,
                    topology_reg_path,
                }
            }

            pub fn create_default_config_builder(&self) -> PlacementConfigBuilder {
                PlacementConfigBuilder::new()
                    .forcefield_path(&self.forcefield_path)
                    .delta_params_path(&self.delta_path)
                    .rotamer_library_path(&self.rotamer_lib_path)
                    .topology_registry_path(&self.topology_reg_path)
                    .s_factor(0.8)
                    .max_iterations(20)
                    .num_solutions(1)
                    .include_input_conformation(true)
                    .final_refinement_iterations(2)
                    .convergence_config(ConvergenceConfig {
                        energy_threshold: 0.01,
                        patience_iterations: 5,
                    })
            }
        }

        fn create_dummy_forcefield_file(dir: &Path) -> PathBuf {
            let path = dir.join("test.ff");
            let mut file = File::create(&path).unwrap();
            write!(
                file,
                r#"
                [globals]
                dielectric_constant = 6.0
                potential_function = "lennard-jones-12-6"
                [vdw]
                N_R = {{ radius = 1.6, well_depth = 0.1 }}
                H_ = {{ radius = 1.0, well_depth = 0.02 }}
                C_31 = {{ radius = 1.8, well_depth = 0.1 }}
                C_32 = {{ radius = 1.8, well_depth = 0.1 }}
                C_33 = {{ radius = 1.8, well_depth = 0.1 }}
                C_R = {{ radius = 1.8, well_depth = 0.1 }}
                O_2 = {{ radius = 1.5, well_depth = 0.2 }}
                [hbond]
            "#
            )
            .unwrap();
            path
        }

        fn create_dummy_delta_file(dir: &Path) -> PathBuf {
            let path = dir.join("test.delta.csv");
            let mut file = File::create(&path).unwrap();
            write!(
                file,
                "residue_type,atom_name,mu,sigma\nALA,CB,0.1,0.05\nLEU,CB,0.1,0.05\n"
            )
            .unwrap();
            path
        }

        fn create_dummy_rotamer_lib_file(dir: &Path) -> PathBuf {
            let path = dir.join("test.rotlib.toml");
            let mut file = File::create(&path).unwrap();
            write!(file, r#"
                [[ALA]]
                atoms = [
                    {{ serial = 1, atom_name = "N", position = [-0.52, 1.36, 0.0], partial_charge = -0.47, force_field_type = "N_R" }},
                    {{ serial = 2, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.07, force_field_type = "C_31" }},
                    {{ serial = 3, atom_name = "C", position = [1.2, -0.1, 0.9], partial_charge = 0.51, force_field_type = "C_R" }},
                    {{ serial = 4, atom_name = "CB", position = [-0.76, -0.8, -1.08], partial_charge = -0.27, force_field_type = "C_33" }},
                    {{ serial = 5, atom_name = "HB1", position = [-0.21, -1.74, -1.16], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 6, atom_name = "HB2", position = [-1.6, -0.4, -1.67], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 7, atom_name = "HB3", position = [-1.2, -1.2, -0.19], partial_charge = 0.09, force_field_type = "H_" }}
                ]
                bonds = [ [1,2], [2,3], [2,4], [4,5], [4,6], [4,7] ]

                [[ALA]]
                atoms = [
                    {{ serial = 1, atom_name = "N", position = [-0.52, 1.36, 0.0], partial_charge = -0.47, force_field_type = "N_R" }},
                    {{ serial = 2, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.07, force_field_type = "C_31" }},
                    {{ serial = 3, atom_name = "C", position = [1.2, -0.1, 0.9], partial_charge = 0.51, force_field_type = "C_R" }},
                    {{ serial = 4, atom_name = "CB", position = [-0.8, -0.9, 1.2], partial_charge = -0.27, force_field_type = "C_33" }},
                    {{ serial = 5, atom_name = "HB1", position = [-0.3, -1.8, 1.3], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 6, atom_name = "HB2", position = [-1.7, -0.5, 1.8], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 7, atom_name = "HB3", position = [-1.3, -1.3, 0.3], partial_charge = 0.09, force_field_type = "H_" }}
                ]
                bonds = [ [1,2], [2,3], [2,4], [4,5], [4,6], [4,7] ]

                [[GLY]]
                atoms = [
                    {{ serial = 1, atom_name = "N", position = [-0.52, 1.36, 0.0], partial_charge = -0.47, force_field_type = "N_R" }},
                    {{ serial = 2, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.02, force_field_type = "C_32" }},
                    {{ serial = 3, atom_name = "C", position = [1.2, -0.1, 0.9], partial_charge = 0.51, force_field_type = "C_R" }},
                    {{ serial = 4, atom_name = "HA1", position = [0.63, -0.47, -0.76], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 5, atom_name = "HA2", position = [-0.63, -0.47, 0.76], partial_charge = 0.09, force_field_type = "H_" }}
                ]
                bonds = [ [1,2], [2,3], [2,4], [2,5] ]

                [[LEU]]
                atoms = [
                    {{ serial = 1, atom_name = "N", position = [-0.52, 1.36, 0.0], partial_charge = -0.47, force_field_type = "N_R" }},
                    {{ serial = 2, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.07, force_field_type = "C_31" }},
                    {{ serial = 3, atom_name = "C", position = [1.2, -0.1, 0.9], partial_charge = 0.51, force_field_type = "C_R" }},
                    {{ serial = 4, atom_name = "CB", position = [-0.8, -0.8, -1.1], partial_charge = -0.18, force_field_type = "C_32" }},
                    {{ serial = 5, atom_name = "HB1", position = [-0.2, -1.7, -1.2], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 6, atom_name = "HB2", position = [-1.6, -0.4, -1.7], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 7, atom_name = "CG", position = [-1.5, -1.5, 0.0], partial_charge = -0.09, force_field_type = "C_31" }},
                    {{ serial = 8, atom_name = "HG", position = [-1.0, -2.5, 0.0], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 9, atom_name = "CD1", position = [-2.5, -1.0, 0.8], partial_charge = -0.27, force_field_type = "C_33" }},
                    {{ serial = 10, atom_name = "HD11", position = [-2.2, -0.5, 1.7], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 11, atom_name = "HD12", position = [-3.0, -1.8, 1.2], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 12, atom_name = "HD13", position = [-3.3, -0.5, 0.2], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 13, atom_name = "CD2", position = [-2.0, -2.5, -1.0], partial_charge = -0.27, force_field_type = "C_33" }},
                    {{ serial = 14, atom_name = "HD21", position = [-1.5, -3.3, -1.5], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 15, atom_name = "HD22", position = [-2.8, -2.9, -0.5], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 16, atom_name = "HD23", position = [-2.5, -2.0, -1.8], partial_charge = 0.09, force_field_type = "H_" }}
                ]
                bonds = [ [2,1], [2,3], [2,4], [4,5], [4,6], [4,7], [7,8], [7,9], [7,13], [9,10], [9,11], [9,12], [13,14], [13,15], [13,16] ]

                [[LEU]]
                atoms = [
                    {{ serial = 1, atom_name = "N", position = [-0.52, 1.36, 0.0], partial_charge = -0.47, force_field_type = "N_R" }},
                    {{ serial = 2, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.07, force_field_type = "C_31" }},
                    {{ serial = 3, atom_name = "C", position = [1.2, -0.1, 0.9], partial_charge = 0.51, force_field_type = "C_R" }},
                    {{ serial = 4, atom_name = "CB", position = [0.5, 1.0, -0.8], partial_charge = -0.18, force_field_type = "C_32" }},
                    {{ serial = 5, atom_name = "HB1", position = [0.9, 1.8, -0.4], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 6, atom_name = "HB2", position = [0.2, 1.2, -1.8], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 7, atom_name = "CG", position = [1.5, 0.0, -0.5], partial_charge = -0.09, force_field_type = "C_31" }},
                    {{ serial = 8, atom_name = "HG", position = [2.0, 0.2, 0.4], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 9, atom_name = "CD1", position = [2.5, -0.5, -1.5], partial_charge = -0.27, force_field_type = "C_33" }},
                    {{ serial = 10, atom_name = "HD11", position = [2.2, -1.5, -1.8], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 11, atom_name = "HD12", position = [3.5, -0.5, -1.2], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 12, atom_name = "HD13", position = [2.8, 0.2, -2.2], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 13, atom_name = "CD2", position = [1.0, -1.0, 0.5], partial_charge = -0.27, force_field_type = "C_33" }},
                    {{ serial = 14, atom_name = "HD21", position = [1.5, -1.8, 0.8], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 15, atom_name = "HD22", position = [0.0, -1.2, 0.8], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 16, atom_name = "HD23", position = [0.8, -0.5, -0.4], partial_charge = 0.09, force_field_type = "H_" }}
                ]
                bonds = [ [2,1], [2,3], [2,4], [4,5], [4,6], [4,7], [7,8], [7,9], [7,13], [9,10], [9,11], [9,12], [13,14], [13,15], [13,16] ]
            "#).unwrap();
            path
        }

        fn create_dummy_topology_reg_file(dir: &Path) -> PathBuf {
            let path = dir.join("test.topology.toml");
            let mut file = File::create(&path).unwrap();
            write!(file, r#"
                [ALA]
                anchor_atoms = ["N", "CA", "C"]
                sidechain_atoms = ["CB", "HB1", "HB2", "HB3"]
                [GLY]
                anchor_atoms = ["N", "CA", "C", "HA2"]
                sidechain_atoms = ["HA1"]
                [LEU]
                anchor_atoms = ["N", "CA", "C"]
                sidechain_atoms = ["CB", "HB1", "HB2", "CG", "HG", "CD1", "HD11", "HD12", "HD13", "CD2", "HD21", "HD22", "HD23"]
            "#).unwrap();
            path
        }
    }

    fn get_sidechain_centroid(system: &MolecularSystem, res_id: ResidueId) -> Point3<f64> {
        let residue = system.residue(res_id).unwrap();
        let mut sidechain_atoms = Vec::new();
        for &atom_id in residue.atoms() {
            if let Some(atom) = system.atom(atom_id) {
                if atom.role == crate::core::models::atom::AtomRole::Sidechain {
                    sidechain_atoms.push(atom);
                }
            }
        }

        if sidechain_atoms.is_empty() {
            return Point3::origin();
        }

        let sum_vec: Vector3<f64> = sidechain_atoms
            .iter()
            .map(|atom| atom.position.coords)
            .sum();
        Point3::from(sum_vec / sidechain_atoms.len() as f64)
    }

    #[test]
    fn test_run_workflow_successfully_on_simple_case() {
        let env = setup::TestEnvironment::new();
        let config = env
            .create_default_config_builder()
            .residues_to_optimize(ResidueSelection::All)
            .build()
            .unwrap();
        let reporter = ProgressReporter::new();

        let system_to_run = env.initial_system.clone();
        let chain_a_id = system_to_run.find_chain_by_id('A').unwrap();
        let initial_ala_centroid = get_sidechain_centroid(
            &system_to_run,
            system_to_run.find_residue_by_id(chain_a_id, 1).unwrap(),
        );
        let initial_gly_centroid = get_sidechain_centroid(
            &system_to_run,
            system_to_run.find_residue_by_id(chain_a_id, 2).unwrap(),
        );

        let result = run(&system_to_run, &config, &reporter);

        assert!(result.is_ok(), "Workflow failed: {:?}", result.err());
        let solutions = result.unwrap();
        assert!(!solutions.is_empty());

        let final_system = &solutions[0].state.system;
        let final_chain_a_id = final_system.find_chain_by_id('A').unwrap();
        let final_ala_centroid = get_sidechain_centroid(
            final_system,
            final_system
                .find_residue_by_id(final_chain_a_id, 1)
                .unwrap(),
        );
        let final_gly_centroid = get_sidechain_centroid(
            final_system,
            final_system
                .find_residue_by_id(final_chain_a_id, 2)
                .unwrap(),
        );

        assert!(
            (final_ala_centroid - initial_ala_centroid).norm() > 1e-6,
            "Alanine sidechain should have moved"
        );
        assert!(
            (final_gly_centroid - initial_gly_centroid).norm() < 1e-6,
            "Glycine sidechain should not have moved"
        );
        assert!(solutions[0].energy < 0.0, "Final energy should be negative");
    }

    #[test]
    fn test_workflow_respects_num_solutions_config() {
        let env = setup::TestEnvironment::new();
        let config = env
            .create_default_config_builder()
            .residues_to_optimize(ResidueSelection::All)
            .num_solutions(3)
            .build()
            .unwrap();
        let reporter = ProgressReporter::new();

        let system_to_run = env.initial_system.clone();
        let solutions = run(&system_to_run, &config, &reporter).unwrap();

        assert!(
            solutions.len() <= 3 && !solutions.is_empty(),
            "Expected between 1 and 3 solutions"
        );
        assert!(solutions.windows(2).all(|w| w[0].energy <= w[1].energy));
    }

    #[test]
    fn test_workflow_respects_residue_selection_list() {
        let env = setup::TestEnvironment::new();
        let selection = ResidueSelection::List {
            include: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 3,
            }],
            exclude: vec![],
        };
        let config = env
            .create_default_config_builder()
            .residues_to_optimize(selection)
            .build()
            .unwrap();
        let reporter = ProgressReporter::new();

        let system_to_run = env.initial_system.clone();
        let chain_a_id = system_to_run.find_chain_by_id('A').unwrap();
        let initial_ala_centroid = get_sidechain_centroid(
            &system_to_run,
            system_to_run.find_residue_by_id(chain_a_id, 1).unwrap(),
        );
        let initial_leu_centroid = get_sidechain_centroid(
            &system_to_run,
            system_to_run.find_residue_by_id(chain_a_id, 3).unwrap(),
        );

        let solutions = run(&system_to_run, &config, &reporter).unwrap();

        let final_system = &solutions[0].state.system;
        let final_chain_a_id = final_system.find_chain_by_id('A').unwrap();
        let final_ala_centroid = get_sidechain_centroid(
            final_system,
            final_system
                .find_residue_by_id(final_chain_a_id, 1)
                .unwrap(),
        );
        let final_leu_centroid = get_sidechain_centroid(
            final_system,
            final_system
                .find_residue_by_id(final_chain_a_id, 3)
                .unwrap(),
        );

        assert!(
            (final_ala_centroid - initial_ala_centroid).norm() < 1e-6,
            "Alanine (not selected) should not move"
        );
        assert!(
            (final_leu_centroid - initial_leu_centroid).norm() > 1e-6,
            "Leucine (selected) should move"
        );
    }

    #[test]
    fn test_clash_resolution_improves_energy() {
        let env = setup::TestEnvironment::new();
        let mut clashing_system = env.initial_system.clone();

        let chain_id = clashing_system.find_chain_by_id('A').unwrap();
        let ala_res_id = clashing_system.find_residue_by_id(chain_id, 1).unwrap();
        let leu_res_id = clashing_system.find_residue_by_id(chain_id, 3).unwrap();

        let ala_cb_id = clashing_system
            .residue(ala_res_id)
            .unwrap()
            .get_first_atom_id_by_name("CB")
            .unwrap();
        let leu_cg_id = clashing_system
            .residue(leu_res_id)
            .unwrap()
            .get_first_atom_id_by_name("CG")
            .unwrap();

        let ala_pos = clashing_system.atom(ala_cb_id).unwrap().position;
        clashing_system.atom_mut(leu_cg_id).unwrap().position = ala_pos;

        let parameterizer = Parameterizer::new(&env.forcefield, &env.topology_registry, 0.0);
        parameterizer
            .parameterize_system(&mut clashing_system)
            .unwrap();

        let config = env
            .create_default_config_builder()
            .residues_to_optimize(ResidueSelection::All)
            .build()
            .unwrap();
        let reporter = ProgressReporter::new();

        let scorer = Scorer::new(&clashing_system, &env.forcefield);
        let initial_energy = scorer
            .score_interaction(
                clashing_system.residue(ala_res_id).unwrap().atoms(),
                clashing_system.residue(leu_res_id).unwrap().atoms(),
            )
            .unwrap();

        assert!(
            initial_energy.vdw > 10.0,
            "Test setup should create a severe clash. Got VdW energy: {}",
            initial_energy.vdw
        );

        let result = run(&clashing_system, &config, &reporter);
        assert!(
            result.is_ok(),
            "Workflow failed with error: {:?}",
            result.err()
        );
        let solutions = result.unwrap();

        assert!(
            !solutions.is_empty(),
            "Workflow should produce at least one solution"
        );
        assert!(
            solutions[0].energy < initial_energy.total(),
            "Final energy ({:.4}) should be much lower than the initial clashing energy ({:.4})",
            solutions[0].energy,
            initial_energy.total()
        );
    }

    #[test]
    fn test_simulated_annealing_runs_when_configured() {
        let env = setup::TestEnvironment::new();
        let sa_config = SimulatedAnnealingConfig {
            initial_temperature: 10.0,
            final_temperature: 0.1,
            cooling_rate: 0.9,
            steps_per_temperature: 2,
        };
        let config = env
            .create_default_config_builder()
            .residues_to_optimize(ResidueSelection::All)
            .simulated_annealing_config(Some(sa_config))
            .build()
            .unwrap();

        let sa_started = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
        let sa_started_clone = sa_started.clone();
        let reporter = ProgressReporter::with_callback(Box::new(move |progress| {
            if let Progress::PhaseStart { name } = progress {
                if name == "Simulated Annealing" {
                    sa_started_clone.store(true, std::sync::atomic::Ordering::Relaxed);
                }
            }
        }));

        let system_to_run = env.initial_system.clone();
        let result = run(&system_to_run, &config, &reporter);

        assert!(result.is_ok());
        assert!(
            sa_started.load(std::sync::atomic::Ordering::Relaxed),
            "Simulated Annealing phase should have started"
        );
    }

    #[test]
    fn test_run_with_no_active_residues_is_a_no_op() {
        let env = setup::TestEnvironment::new();
        let selection = ResidueSelection::List {
            include: vec![],
            exclude: vec![],
        };
        let config = env
            .create_default_config_builder()
            .residues_to_optimize(selection)
            .build()
            .unwrap();
        let reporter = ProgressReporter::new();

        let system_to_run = env.initial_system.clone();
        let solutions = run(&system_to_run, &config, &reporter).unwrap();

        assert_eq!(solutions.len(), 1);
        let final_system = &solutions[0].state.system;

        assert_eq!(
            final_system.atoms_iter().count(),
            env.initial_system.atoms_iter().count()
        );
        for (id, initial_atom) in env.initial_system.atoms_iter() {
            let final_atom = final_system.atom(id).unwrap();
            assert!((final_atom.position - initial_atom.position).norm() < 1e-9);
        }
    }
}
