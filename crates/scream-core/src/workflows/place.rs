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

            let mut best_idx = state.working_state.rotamers[&res_id];
            let mut best_score = state.current_optimization_score;

            let original_system = state.working_state.system.clone();
            let original_rotamers = state.working_state.rotamers.clone();

            for idx in 0..rotamers.len() {
                if idx == best_idx {
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

            if best_idx != original_rotamers[&res_id] {
                changed_in_cycle = true;
                update_rotamers_in_state(state, res_id, best_idx, context)?;
                state.current_optimization_score = best_score;
                state.submit_current_solution();
            }
        }
        if !changed_in_cycle {
            info!(iteration = i + 1, "Refinement converged.");
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
