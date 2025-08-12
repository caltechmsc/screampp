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
use crate::engine::state::{InitialState, OptimizationState, Solution, SolutionState};
use crate::engine::tasks;
use rand::prelude::IteratorRandom;
use rand::{Rng, seq::SliceRandom, thread_rng};
use std::collections::{HashMap, HashSet};
use tracing::{info, instrument};

const CLASH_THRESHOLD_KCAL_MOL: f64 = 25.0;

#[derive(Debug, Clone)]
pub struct PlacementResult {
    pub initial_state: InitialState,
    pub solutions: Vec<Solution>,
}

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
    context: &OptimizationContext<PlacementConfig>,
    el_cache: &ELCache,
) -> Result<(), EngineError> {
    context.reporter.report(Progress::PhaseStart {
        name: "Clash Resolution",
    });
    info!("Starting iterative clash resolution loop.");

    let max_iter = context.config.optimization.max_iterations;
    let energy_threshold = context.config.optimization.convergence.energy_threshold;
    let patience_iterations = context.config.optimization.convergence.patience_iterations;

    let mut last_best_score = state.best_energy();
    let mut iterations_without_improvement = 0;

    for iter in 0..max_iter {
        let clashes = tasks::clash_detection::run(
            &state.working_state.system,
            context.forcefield,
            active_residues,
            CLASH_THRESHOLD_KCAL_MOL,
            context.reporter,
        )?;

        context.reporter.report(Progress::StatusUpdate {
            text: format!("Pass {}/{}, Clashes: {}", iter + 1, max_iter, clashes.len()),
        });

        if clashes.is_empty() {
            info!(
                iteration = iter + 1,
                "Convergence reached: no clashes found."
            );
            context.reporter.report(Progress::Message(format!(
                "Converged after {} iterations (no clashes).",
                iter + 1
            )));
            break;
        }

        let worst_clash = &clashes[0];
        let (res_a_id, res_b_id) = (worst_clash.residue_a, worst_clash.residue_b);

        let doublet_result = tasks::doublet_optimization::run(
            res_a_id,
            res_b_id,
            &state.working_state.system,
            el_cache,
            context,
        )?;

        update_rotamers_in_state(state, res_a_id, doublet_result.rotamer_idx_a, context)?;
        update_rotamers_in_state(state, res_b_id, doublet_result.rotamer_idx_b, context)?;

        let new_score = calculate_optimization_score(state, active_residues, context, el_cache)?;
        state.current_optimization_score = new_score;
        state.submit_current_solution();

        let current_best_score = state.best_energy();
        let improvement = last_best_score - current_best_score;

        if improvement < energy_threshold {
            iterations_without_improvement += 1;
        } else {
            iterations_without_improvement = 0;
        }

        last_best_score = current_best_score;

        if iterations_without_improvement >= patience_iterations {
            info!(
                iteration = iter + 1,
                patience = patience_iterations,
                threshold = energy_threshold,
                "Convergence reached: energy stabilized."
            );
            context.reporter.report(Progress::Message(format!(
                "Converged after {} iterations (energy stabilized).",
                iter + 1
            )));
            break;
        }
    }

    context.reporter.report(Progress::PhaseFinish);
    Ok(())
}

fn run_simulated_annealing(
    state: &mut OptimizationState,
    active_residues: &HashSet<ResidueId>,
    context: &OptimizationContext<PlacementConfig>,
    el_cache: &ELCache,
) -> Result<(), EngineError> {
    let sa_config = match &context.config.optimization.simulated_annealing {
        Some(config) => config,
        None => return Ok(()),
    };
    context.reporter.report(Progress::PhaseStart {
        name: "Simulated Annealing",
    });
    info!("Starting Simulated Annealing.");

    let mut rng = thread_rng();
    let mut current_temp = sa_config.initial_temperature;
    let active_residue_vec: Vec<_> = active_residues.iter().cloned().collect();

    while current_temp > sa_config.final_temperature {
        context.reporter.report(Progress::StatusUpdate {
            text: format!("SA Temp: {:.2}", current_temp),
        });
        for _ in 0..sa_config.steps_per_temperature {
            let res_id = *active_residue_vec.choose(&mut rng).unwrap();
            let res_type = state
                .working_state
                .system
                .residue(res_id)
                .unwrap()
                .residue_type
                .unwrap();
            let rotamers = context.rotamer_library.get_rotamers_for(res_type).unwrap();
            if rotamers.len() <= 1 {
                continue;
            }

            let current_rot_idx = state.working_state.rotamers[&res_id];
            let new_rot_idx = (0..rotamers.len())
                .filter(|&i| i != current_rot_idx)
                .choose(&mut rng)
                .unwrap();

            let original_system = state.working_state.system.clone();
            let original_rotamers = state.working_state.rotamers.clone();
            let original_score = state.current_optimization_score;

            update_rotamers_in_state(state, res_id, new_rot_idx, context)?;
            let new_score =
                calculate_optimization_score(state, active_residues, context, el_cache)?;

            let delta_e = new_score - original_score;
            if delta_e < 0.0 || rng.r#gen::<f64>() < (-delta_e / current_temp).exp() {
                state.current_optimization_score = new_score;
                state.submit_current_solution();
            } else {
                state.working_state.system = original_system;
                state.working_state.rotamers = original_rotamers;
                state.current_optimization_score = original_score;
            }
        }
        current_temp *= sa_config.cooling_rate;
    }

    context.reporter.report(Progress::PhaseFinish);
    Ok(())
}

fn final_refinement(
    state: &mut OptimizationState,
    active_residues: &HashSet<ResidueId>,
    context: &OptimizationContext<PlacementConfig>,
    el_cache: &ELCache,
) -> Result<(), EngineError> {
    let iterations = context.config.optimization.final_refinement_iterations;
    if iterations == 0 {
        return Ok(());
    }

    context.reporter.report(Progress::PhaseStart {
        name: "Final Refinement",
    });
    info!("Starting final refinement (singlet optimization).");

    for i in 0..iterations {
        context.reporter.report(Progress::StatusUpdate {
            text: format!("Pass {}/{}", i + 1, iterations),
        });

        context.reporter.report(Progress::TaskStart {
            total: active_residues.len() as u64,
        });

        let mut changed_in_cycle = false;
        let mut residues_to_process: Vec<_> = active_residues.iter().cloned().collect();
        residues_to_process.shuffle(&mut thread_rng());

        for res_id in residues_to_process {
            let res_type = state
                .working_state
                .system
                .residue(res_id)
                .unwrap()
                .residue_type
                .unwrap();
            let rotamers = context.rotamer_library.get_rotamers_for(res_type).unwrap();

            let current_rot_idx = state.working_state.rotamers[&res_id];
            let mut best_idx = current_rot_idx;
            let mut best_score = state.current_optimization_score;

            let original_system = state.working_state.system.clone();
            let original_rotamers = state.working_state.rotamers.clone();

            for idx in 0..rotamers.len() {
                if idx == current_rot_idx {
                    continue;
                }

                update_rotamers_in_state(state, res_id, idx, context)?;
                let score =
                    calculate_optimization_score(state, active_residues, context, el_cache)?;

                if score < best_score {
                    best_score = score;
                    best_idx = idx;
                }

                state.working_state.system = original_system.clone();
                state.working_state.rotamers = original_rotamers.clone();
            }

            if best_idx != current_rot_idx {
                changed_in_cycle = true;
                update_rotamers_in_state(state, res_id, best_idx, context)?;
                state.current_optimization_score = best_score;
                state.submit_current_solution();
            }

            context
                .reporter
                .report(Progress::TaskIncrement { amount: 1 });
        }

        context.reporter.report(Progress::TaskFinish);

        if !changed_in_cycle {
            info!(
                iteration = i + 1,
                "Refinement converged as no changes occurred in this pass."
            );
            context
                .reporter
                .report(Progress::Message(format!("Converged after pass {}", i + 1)));
            break;
        }
    }

    context.reporter.report(Progress::PhaseFinish);
    Ok(())
}

fn finalize_results(
    state: OptimizationState,
    initial_state: InitialState,
    energy_offset: EnergyTerm,
    num_solutions: usize,
) -> PlacementResult {
    let mut solutions = state.into_sorted_solutions();

    let should_include_initial = if solutions.len() < num_solutions {
        true
    } else if let Some(worst_solution) = solutions.last() {
        initial_state.optimization_score < worst_solution.optimization_score
    } else {
        true
    };

    if should_include_initial {
        let initial_as_solution = Solution {
            total_energy: initial_state.total_energy,
            optimization_score: initial_state.optimization_score,
            state: SolutionState {
                system: initial_state.system.clone(),
                rotamers: HashMap::new(),
            },
        };
        solutions.push(initial_as_solution);
    }

    solutions.sort_by(|a, b| {
        a.optimization_score
            .partial_cmp(&b.optimization_score)
            .unwrap()
    });
    solutions.dedup_by(|a, b| (a.optimization_score - b.optimization_score).abs() < 1e-6);
    solutions.truncate(num_solutions);

    for sol in &mut solutions {
        sol.total_energy = sol.optimization_score + energy_offset.total();
    }

    PlacementResult {
        initial_state,
        solutions,
    }
}

fn update_rotamers_in_state(
    state: &mut OptimizationState,
    res_id: ResidueId,
    new_rot_idx: usize,
    context: &OptimizationContext<PlacementConfig>,
) -> Result<(), EngineError> {
    let res_type = state
        .working_state
        .system
        .residue(res_id)
        .unwrap()
        .residue_type
        .unwrap();
    let rotamer = &context.rotamer_library.get_rotamers_for(res_type).unwrap()[new_rot_idx];
    let res_name = res_type.to_three_letter();
    let topology = context.topology_registry.get(res_name).unwrap();

    place_rotamer_on_system(&mut state.working_state.system, res_id, rotamer, topology)?;
    state.working_state.rotamers.insert(res_id, new_rot_idx);
    Ok(())
}

fn calculate_optimization_score(
    state: &OptimizationState,
    active_residues: &HashSet<ResidueId>,
    context: &OptimizationContext<PlacementConfig>,
    el_cache: &ELCache,
) -> Result<f64, EngineError> {
    let mut el_sum = EnergyTerm::default();
    for (&res_id, &rot_idx) in &state.working_state.rotamers {
        let res_type = state
            .working_state
            .system
            .residue(res_id)
            .unwrap()
            .residue_type
            .unwrap();
        el_sum += *el_cache
            .get(res_id, res_type, rot_idx)
            .unwrap_or(&EnergyTerm::default());
    }

    let interaction = tasks::interaction_energy::run(
        &state.working_state.system,
        context.forcefield,
        active_residues,
    )?;

    Ok((el_sum + interaction).total())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::parameterization::Parameterizer;
    use crate::core::forcefield::params::Forcefield;
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
                    .unwrap();

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
            pub forcefield_path: PathBuf,
            pub delta_path: PathBuf,
            pub rotamer_lib_path: PathBuf,
            pub topology_registry_path: PathBuf,
            pub initial_system: MolecularSystem,
        }

        impl TestEnvironment {
            pub fn new() -> Self {
                let temp_dir = tempdir().unwrap();
                let forcefield_path = create_dummy_forcefield_file(temp_dir.path());
                let delta_path = create_dummy_delta_file(temp_dir.path());
                let rotamer_lib_path = create_dummy_rotamer_lib_file(temp_dir.path());
                let topology_registry_path = create_dummy_topology_registry_file(temp_dir.path());

                let forcefield = Forcefield::load(&forcefield_path, &delta_path).unwrap();
                let topology_registry = TopologyRegistry::load(&topology_registry_path).unwrap();

                let mut system = TestSystemBuilder::new()
                    .add_residue(ResidueType::Alanine, 1)
                    .add_residue(ResidueType::Glycine, 2)
                    .add_residue(ResidueType::Leucine, 3)
                    .build();

                let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 0.8);
                parameterizer.parameterize_system(&mut system).unwrap();

                Self {
                    temp_dir,
                    forcefield_path,
                    delta_path,
                    rotamer_lib_path,
                    topology_registry_path,
                    initial_system: system,
                }
            }

            pub fn create_default_config_builder(&self) -> PlacementConfigBuilder {
                PlacementConfigBuilder::new()
                    .forcefield_path(&self.forcefield_path)
                    .delta_params_path(&self.delta_path)
                    .rotamer_library_path(&self.rotamer_lib_path)
                    .topology_registry_path(&self.topology_registry_path)
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
                    {{ serial = 4, atom_name = "CB", position = [-1.5, -1.5, -1.5], partial_charge = -0.27, force_field_type = "C_33" }},
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
                    {{ serial = 4, atom_name = "CB", position = [-2.0, -2.0, -2.0], partial_charge = -0.18, force_field_type = "C_32" }},
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

        fn create_dummy_topology_registry_file(dir: &Path) -> PathBuf {
            let path = dir.join("test.topology.toml");
            let mut file = File::create(&path).unwrap();
            write!(
                file,
                r#"
                [ALA]
                anchor_atoms = ["N", "CA", "C"]
                sidechain_atoms = ["CB", "HB1", "HB2", "HB3"]
                [GLY]
                anchor_atoms = ["N", "CA", "C"]
                sidechain_atoms = []
                [LEU]
                anchor_atoms = ["N", "CA", "C"]
                sidechain_atoms = ["CB", "HB1", "HB2", "CG", "HG"]
            "#
            )
            .unwrap();
            path
        }
    }

    fn get_sidechain_centroid(system: &MolecularSystem, res_id: ResidueId) -> Point3<f64> {
        let sc_atom_ids: Vec<_> = system
            .residue(res_id)
            .unwrap()
            .atoms()
            .iter()
            .filter(|id| {
                system.atom(**id).unwrap().role == crate::core::models::atom::AtomRole::Sidechain
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
        let ala_res_id = system_to_run
            .find_residue_by_id(system_to_run.find_chain_by_id('A').unwrap(), 1)
            .unwrap();
        let gly_res_id = system_to_run
            .find_residue_by_id(system_to_run.find_chain_by_id('A').unwrap(), 2)
            .unwrap();

        let initial_ala_centroid = get_sidechain_centroid(&system_to_run, ala_res_id);
        let initial_gly_centroid = get_sidechain_centroid(&system_to_run, gly_res_id);

        let result = run(&system_to_run, &config, &reporter);
        assert!(result.is_ok(), "Workflow failed: {:?}", result.err());

        let placement_result = result.unwrap();
        assert!(!placement_result.solutions.is_empty());

        let best_solution = &placement_result.solutions[0];
        let final_system = &best_solution.state.system;

        let final_ala_centroid = get_sidechain_centroid(final_system, ala_res_id);
        let final_gly_centroid = get_sidechain_centroid(final_system, gly_res_id);

        assert!(
            (final_ala_centroid - initial_ala_centroid).norm() > 1e-6,
            "Alanine sidechain should have moved"
        );
        assert!(
            (final_gly_centroid - initial_gly_centroid).norm() < 1e-6,
            "Glycine has no sidechain to move"
        );
        assert!(
            best_solution.total_energy < placement_result.initial_state.total_energy,
            "Final total energy should be lower than initial"
        );
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

        let solutions = run(&env.initial_system, &config, &reporter)
            .unwrap()
            .solutions;

        assert!(
            solutions.len() <= 3 && !solutions.is_empty(),
            "Expected between 1 and 3 solutions"
        );
        assert!(
            solutions
                .windows(2)
                .all(|w| w[0].optimization_score <= w[1].optimization_score),
            "Solutions should be sorted by optimization score"
        );
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

        let chain_a_id = env.initial_system.find_chain_by_id('A').unwrap();
        let ala_res_id = env
            .initial_system
            .find_residue_by_id(chain_a_id, 1)
            .unwrap();
        let leu_res_id = env
            .initial_system
            .find_residue_by_id(chain_a_id, 3)
            .unwrap();

        let initial_ala_centroid = get_sidechain_centroid(&env.initial_system, ala_res_id);
        let initial_leu_centroid = get_sidechain_centroid(&env.initial_system, leu_res_id);

        let solutions = run(&env.initial_system, &config, &reporter)
            .unwrap()
            .solutions;

        let final_system = &solutions[0].state.system;
        let final_ala_centroid = get_sidechain_centroid(final_system, ala_res_id);
        let final_leu_centroid = get_sidechain_centroid(final_system, leu_res_id);

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
        let scorer = crate::core::forcefield::scoring::Scorer::new(&clashing_system, &forcefield);

        let initial_interaction = scorer
            .score_interaction(
                &clashing_system
                    .residue(ala_res_id)
                    .unwrap()
                    .atoms()
                    .iter()
                    .map(|&id| id)
                    .collect::<Vec<_>>(),
                &clashing_system
                    .residue(leu_res_id)
                    .unwrap()
                    .atoms()
                    .iter()
                    .map(|&id| id)
                    .collect::<Vec<_>>(),
            )
            .unwrap();

        assert!(
            initial_interaction.vdw > 10.0,
            "Test setup should create a severe clash. Got VdW energy: {}",
            initial_interaction.vdw
        );

        let result = run(&clashing_system, &config, &reporter).unwrap();

        assert!(
            !result.solutions.is_empty(),
            "Workflow should produce at least one solution"
        );
        assert!(
            result.solutions[0].total_energy < result.initial_state.total_energy,
            "Final total energy ({:.4}) should be much lower than the initial clashing energy ({:.4})",
            result.solutions[0].total_energy,
            result.initial_state.total_energy
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

        let result = run(&env.initial_system, &config, &reporter);

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

        let result = run(&env.initial_system, &config, &reporter).unwrap();

        assert_eq!(result.solutions.len(), 1);
        let final_system = &result.solutions[0].state.system;

        assert!(
            (result.initial_state.total_energy - result.solutions[0].total_energy).abs() < 1e-9
        );

        for (id, initial_atom) in env.initial_system.atoms_iter() {
            let final_atom = final_system.atom(id).unwrap();
            assert!((final_atom.position - initial_atom.position).norm() < 1e-9);
        }
    }
}
