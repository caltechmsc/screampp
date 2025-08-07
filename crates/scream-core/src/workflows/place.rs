use crate::core::forcefield::parameterization::Parameterizer;
use crate::core::forcefield::params::Forcefield;
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
use crate::engine::state::{OptimizationState, Solution};
use crate::engine::tasks;
use rand::Rng;
use rand::seq::SliceRandom;
use rand::thread_rng;
use std::collections::{HashMap, HashSet};
use tracing::{debug, info, instrument, warn};

const CLASH_THRESHOLD: f64 = 25.0; // Clash threshold (kcal/mol)

#[instrument(skip_all, name = "placement_workflow")]
pub fn run(
    initial_system: &MolecularSystem,
    config: &PlacementConfig,
    reporter: &ProgressReporter,
) -> Result<Vec<Solution>, EngineError> {
    // --- Phase 0: Setup, Resource Loading, and System Initialization ---
    reporter.report(Progress::PhaseStart { name: "Setup" });
    info!("Starting workflow setup: loading resources and parameterizing system.");

    let forcefield = Forcefield::load(
        &config.forcefield.forcefield_path,
        &config.forcefield.delta_params_path,
    )?;
    let topology_registry = TopologyRegistry::load(&config.topology_registry_path)?;

    let parameterizer =
        Parameterizer::new(&forcefield, &topology_registry, config.forcefield.s_factor);

    let mut rotamer_library = RotamerLibrary::load(
        &config.sampling.rotamer_library_path,
        &topology_registry,
        &forcefield,
        config.forcefield.s_factor,
    )?;
    let mut working_system = initial_system.clone();

    let active_residues = resolve_selection_to_ids(
        &working_system,
        &config.residues_to_optimize,
        &rotamer_library,
    )?;

    if config.optimization.include_input_conformation {
        info!("Including original side-chain conformations in the rotamer library.");
        rotamer_library.include_system_conformations(
            &working_system,
            &active_residues,
            &topology_registry,
        );
    }

    info!("Parameterizing the input molecular system...");
    parameterizer.parameterize_system(&mut working_system)?;

    reporter.report(Progress::PhaseFinish);

    // --- Phase 1: Pre-computation ---
    let context = OptimizationContext::new(
        &working_system,
        &forcefield,
        reporter,
        config,
        &rotamer_library,
        &topology_registry,
    );
    let el_cache = precompute_el_energies(&context)?;

    // --- Phase 2: Initialization ---
    let mut state = initialize_state(&working_system, &active_residues, &context, &el_cache)?;

    // --- Phase 3: Clash-Driven Optimization ---
    resolve_clashes_iteratively(&mut state, &active_residues, &context, &el_cache)?;

    // --- Phase 4: Optional Global Search (Simulated Annealing) ---
    if context.config.optimization.simulated_annealing.is_some() {
        run_simulated_annealing(&mut state, &active_residues, &context, &el_cache)?;
    }

    // --- Phase 5: Final Refinement (Singlet Optimization) ---
    final_refinement(&mut state, &active_residues, &context, &el_cache)?;

    // --- Phase 6: Finalization ---
    let sorted_solutions = state.into_sorted_solutions();
    info!(
        "Workflow complete. Returning {} sorted solutions.",
        sorted_solutions.len()
    );
    Ok(sorted_solutions)
}

#[instrument(skip_all, name = "el_energy_precomputation")]
fn precompute_el_energies<'a>(
    context: &OptimizationContext<'a, PlacementConfig>,
) -> Result<ELCache, EngineError> {
    context.reporter.report(Progress::PhaseStart {
        name: "EL Pre-computation",
    });
    let el_cache = tasks::el_energy::run(context)?;
    context.reporter.report(Progress::PhaseFinish);
    Ok(el_cache)
}

#[instrument(skip_all, name = "state_initialization")]
fn initialize_state<'a>(
    working_system: &MolecularSystem,
    active_residues: &HashSet<ResidueId>,
    context: &OptimizationContext<'a, PlacementConfig>,
    el_cache: &ELCache,
) -> Result<OptimizationState, EngineError> {
    context.reporter.report(Progress::PhaseStart {
        name: "Initialization",
    });
    info!("Initializing system with ground-state rotamers.");

    let mut initial_rotamers = HashMap::new();
    let mut working_system = working_system.clone();

    for &residue_id in active_residues {
        let residue_type = working_system
            .residue(residue_id)
            .and_then(|r| r.residue_type)
            .ok_or(EngineError::Internal(format!(
                "Active residue {:?} has no residue type.",
                residue_id
            )))?;

        if let Some((ground_state_idx, _)) =
            el_cache.find_ground_state_for(residue_id, residue_type)
        {
            let rotamer = &context
                .rotamer_library
                .get_rotamers_for(residue_type)
                .unwrap()[ground_state_idx];

            let residue_name = residue_type.to_three_letter();
            let topology = context.topology_registry.get(residue_name).ok_or_else(|| {
                EngineError::TopologyNotFound {
                    residue_name: residue_name.to_string(),
                }
            })?;

            place_rotamer_on_system(&mut working_system, residue_id, rotamer, topology)?;
        } else {
            warn!(
                "No ground state found in EL cache for residue {:?}. It may not be placed correctly initially.",
                residue_id
            );
        }
    }

    let initial_energy_term = tasks::total_energy::run(
        &working_system,
        context.forcefield,
        active_residues,
        &initial_rotamers,
        el_cache,
    )?;
    let initial_energy = initial_energy_term.total();
    info!(initial_energy, "Initial system energy calculated.");

    let state = OptimizationState::new(
        working_system,
        initial_rotamers,
        initial_energy,
        context.config.optimization.num_solutions,
    );

    context.reporter.report(Progress::PhaseFinish);
    Ok(state)
}

#[instrument(skip_all, name = "clash_resolution_loop")]
fn resolve_clashes_iteratively<'a>(
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
                ("CB", "CG"),
                ("CG", "HG"),
                ("C", "O"),
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
            parameterizer: Parameterizer,
            last_c_id: Option<AtomId>,
            chain_id: ChainId,
        }

        impl TestSystemBuilder {
            pub fn new(forcefield: Forcefield) -> Self {
                let mut system = MolecularSystem::new();
                let chain_id = system.add_chain('A', ChainType::Protein);
                let parameterizer = Parameterizer::new(forcefield, 0.0);
                Self {
                    system,
                    parameterizer,
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

                    self.parameterizer
                        .parameterize_atom(&mut atom, residue_type.to_three_letter())
                        .unwrap_or_else(|e| {
                            panic!(
                                "Failed to parameterize atom {} in test setup: {}",
                                atom_template.name, e
                            )
                        });

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
            pub forcefield_path: PathBuf,
            pub delta_path: PathBuf,
            pub rotamer_lib_path: PathBuf,
            pub placement_reg_path: PathBuf,
            pub initial_system: MolecularSystem,
        }

        impl TestEnvironment {
            pub fn new() -> Self {
                let temp_dir = tempdir().unwrap();
                let forcefield_path = create_dummy_forcefield_file(temp_dir.path());
                let delta_path = create_dummy_delta_file(temp_dir.path());
                let rotamer_lib_path = create_dummy_rotamer_lib_file(temp_dir.path());
                let placement_reg_path = create_dummy_placement_reg_file(temp_dir.path());

                let forcefield = Forcefield::load(&forcefield_path, &delta_path)
                    .expect("Failed to load dummy forcefield in test setup");

                let system = TestSystemBuilder::new(forcefield)
                    .add_residue(ResidueType::Alanine, 1)
                    .add_residue(ResidueType::Glycine, 2)
                    .add_residue(ResidueType::Leucine, 3)
                    .build();

                Self {
                    temp_dir,
                    forcefield_path,
                    delta_path,
                    rotamer_lib_path,
                    placement_reg_path,
                    initial_system: system,
                }
            }

            pub fn create_default_config_builder(&self) -> PlacementConfigBuilder {
                PlacementConfigBuilder::new()
                    .forcefield_path(&self.forcefield_path)
                    .delta_params_path(&self.delta_path)
                    .rotamer_library_path(&self.rotamer_lib_path)
                    .placement_registry_path(&self.placement_reg_path)
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
            write!(
                file,
                r#"
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

                [[LEU]]
                atoms = [
                    {{ serial = 1, atom_name = "N", position = [-0.52, 1.36, 0.0], partial_charge = -0.47, force_field_type = "N_R" }},
                    {{ serial = 2, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.07, force_field_type = "C_31" }},
                    {{ serial = 3, atom_name = "C", position = [1.2, -0.1, 0.9], partial_charge = 0.51, force_field_type = "C_R" }},
                    {{ serial = 4, atom_name = "CB", position = [-0.8, -0.8, -1.1], partial_charge = -0.18, force_field_type = "C_32" }},
                    {{ serial = 5, atom_name = "HB1", position = [-0.2, -1.7, -1.2], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 6, atom_name = "HB2", position = [-1.6, -0.4, -1.7], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 7, atom_name = "CG", position = [-1.5, -1.5, 0.0], partial_charge = -0.09, force_field_type = "C_31" }},
                    {{ serial = 8, atom_name = "HG", position = [-1.0, -2.5, 0.0], partial_charge = 0.09, force_field_type = "H_" }}
                ]
                bonds = [ [1,2], [2,3], [2,4], [4,5], [4,6], [4,7], [7,8] ]

                [[LEU]]
                atoms = [
                    {{ serial = 1, atom_name = "N", position = [-0.52, 1.36, 0.0], partial_charge = -0.47, force_field_type = "N_R" }},
                    {{ serial = 2, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.07, force_field_type = "C_31" }},
                    {{ serial = 3, atom_name = "C", position = [1.2, -0.1, 0.9], partial_charge = 0.51, force_field_type = "C_R" }},
                    {{ serial = 4, atom_name = "CB", position = [0.5, 1.0, -0.8], partial_charge = -0.18, force_field_type = "C_32" }},
                    {{ serial = 5, atom_name = "HB1", position = [0.9, 1.8, -0.4], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 6, atom_name = "HB2", position = [0.2, 1.2, -1.8], partial_charge = 0.09, force_field_type = "H_" }},
                    {{ serial = 7, atom_name = "CG", position = [1.5, 0.0, -0.5], partial_charge = -0.09, force_field_type = "C_31" }},
                    {{ serial = 8, atom_name = "HG", position = [2.0, 0.2, 0.4], partial_charge = 0.09, force_field_type = "H_" }}
                ]
                bonds = [ [1,2], [2,3], [2,4], [4,5], [4,6], [4,7], [7,8] ]
                "#
            ).unwrap();
            path
        }

        fn create_dummy_placement_reg_file(dir: &Path) -> PathBuf {
            let path = dir.join("test.placement.toml");
            let mut file = File::create(&path).unwrap();
            write!(
                file,
                r#"
                [ALA]
                anchor_atoms = ["N", "CA", "C"]
                sidechain_atoms = ["CB", "HB1", "HB2", "HB3"]
                exact_match_atoms = []
                connection_points = []
                [LEU]
                anchor_atoms = ["N", "CA", "C"]
                sidechain_atoms = ["CB", "HB1", "HB2", "CG", "HG"]
                exact_match_atoms = []
                connection_points = []
            "#
            )
            .unwrap();
            path
        }
    }

    fn get_sidechain_centroid(system: &MolecularSystem, res_id: ResidueId) -> Point3<f64> {
        let residue = system.residue(res_id).unwrap();
        let sc_atom_ids: Vec<_> = residue
            .atoms()
            .iter()
            .filter(|id| {
                let atom_name = &system.atom(**id).unwrap().name;
                !["N", "H", "CA", "HA", "HA1", "HA2", "C", "O"].contains(&atom_name.as_str())
            })
            .copied()
            .collect();

        if sc_atom_ids.is_empty() {
            return Point3::origin();
        }

        let sum_vec: Vector3<f64> = sc_atom_ids
            .iter()
            .map(|id| system.atom(*id).unwrap().position.coords)
            .sum();
        Point3::from(sum_vec / sc_atom_ids.len() as f64)
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
        let initial_ala_centroid = get_sidechain_centroid(
            &system_to_run,
            system_to_run
                .find_residue_by_id(system_to_run.find_chain_by_id('A').unwrap(), 1)
                .unwrap(),
        );
        let initial_gly_centroid = get_sidechain_centroid(
            &system_to_run,
            system_to_run
                .find_residue_by_id(system_to_run.find_chain_by_id('A').unwrap(), 2)
                .unwrap(),
        );

        let result = run(&system_to_run, &config, &reporter);

        assert!(result.is_ok(), "Workflow failed: {:?}", result.err());
        let solutions = result.unwrap();
        assert!(!solutions.is_empty());

        let final_system = &solutions[0].state.system;
        let final_ala_centroid = get_sidechain_centroid(
            final_system,
            final_system
                .find_residue_by_id(final_system.find_chain_by_id('A').unwrap(), 1)
                .unwrap(),
        );
        let final_gly_centroid = get_sidechain_centroid(
            final_system,
            final_system
                .find_residue_by_id(final_system.find_chain_by_id('A').unwrap(), 2)
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
        let initial_ala_centroid = get_sidechain_centroid(
            &system_to_run,
            system_to_run
                .find_residue_by_id(system_to_run.find_chain_by_id('A').unwrap(), 1)
                .unwrap(),
        );
        let initial_leu_centroid = get_sidechain_centroid(
            &system_to_run,
            system_to_run
                .find_residue_by_id(system_to_run.find_chain_by_id('A').unwrap(), 3)
                .unwrap(),
        );

        let solutions = run(&system_to_run, &config, &reporter).unwrap();

        let final_system = &solutions[0].state.system;
        let final_ala_centroid = get_sidechain_centroid(
            final_system,
            final_system
                .find_residue_by_id(final_system.find_chain_by_id('A').unwrap(), 1)
                .unwrap(),
        );
        let final_leu_centroid = get_sidechain_centroid(
            final_system,
            final_system
                .find_residue_by_id(final_system.find_chain_by_id('A').unwrap(), 3)
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
        let leu_cb_id = clashing_system
            .residue(leu_res_id)
            .unwrap()
            .get_first_atom_id_by_name("CB")
            .unwrap();

        let ala_pos = clashing_system.atom(ala_cb_id).unwrap().position;
        clashing_system.atom_mut(leu_cb_id).unwrap().position = ala_pos;

        let config = env
            .create_default_config_builder()
            .residues_to_optimize(ResidueSelection::All)
            .build()
            .unwrap();
        let reporter = ProgressReporter::new();
        let forcefield = Forcefield::load(&env.forcefield_path, &env.delta_path).unwrap();
        let scorer = Scorer::new(&clashing_system, &forcefield);

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
