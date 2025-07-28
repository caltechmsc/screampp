use crate::core::forcefield::params::Forcefield;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::core::rotamers::library::RotamerLibrary;
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
use tracing::{info, instrument, warn};

const CLASH_THRESHOLD: f64 = 25.0; // Clash threshold (kcal/mol)
const RMS_THRESHOLD: f64 = 0.1; // RMSD threshold for duplicate rotamers (in Angstroms)

#[instrument(skip_all, name = "placement_workflow")]
pub fn run(
    initial_system: &MolecularSystem,
    config: &PlacementConfig,
    reporter: &ProgressReporter,
) -> Result<Vec<Solution>, EngineError> {
    // --- Phase 0 & 1: Setup & Pre-computation ---
    let (forcefield, rotamer_library, active_residues) = setup(initial_system, config, reporter)?;
    let context = OptimizationContext::new(
        initial_system,
        &forcefield,
        reporter,
        config,
        &rotamer_library,
    );
    let el_cache = precompute_el_energies(&context, reporter)?;

    // --- Phase 2: Initialization ---
    let mut state = initialize_state(
        initial_system,
        &active_residues,
        &context,
        &el_cache,
        reporter,
    )?;

    // --- Phase 3: Clash-Driven Optimization ---
    resolve_clashes_iteratively(&mut state, &active_residues, &context, &el_cache, reporter)?;

    // --- Phase 4: Optional Global Search (Simulated Annealing) ---
    if context.config.optimization.simulated_annealing.is_some() {
        run_simulated_annealing(&mut state, &active_residues, &context, &el_cache, reporter)?;
    }

    // --- Phase 5: Final Refinement (Singlet Optimization) ---
    final_refinement(&mut state, &active_residues, &context, &el_cache, reporter)?;

    // --- Phase 6: Finalization ---
    let sorted_solutions = state.into_sorted_solutions();
    info!(
        "Workflow complete. Returning {} sorted solutions.",
        sorted_solutions.len()
    );
    Ok(sorted_solutions)
}

#[instrument(skip_all, name = "workflow_setup")]
fn setup<'a>(
    initial_system: &MolecularSystem,
    config: &'a PlacementConfig,
    reporter: &ProgressReporter,
) -> Result<(Forcefield, RotamerLibrary, HashSet<ResidueId>), EngineError> {
    reporter.report(Progress::PhaseStart { name: "Setup" });
    info!("Starting workflow setup: loading resources.");

    let forcefield = Forcefield::load(
        &config.forcefield.forcefield_path,
        &config.forcefield.delta_params_path,
    )?;

    let mut rotamer_library = RotamerLibrary::load(
        &config.sampling.rotamer_library_path,
        &config.sampling.placement_registry_path,
        &forcefield,
        config.forcefield.s_factor,
    )?;

    let active_residues = resolve_selection_to_ids(
        initial_system,
        &config.residues_to_optimize,
        &rotamer_library,
    )?;

    if config.optimization.include_input_conformation {
        info!("Including original side-chain conformations in the rotamer library.");
        rotamer_library.include_system_conformations(
            initial_system,
            &active_residues,
            RMS_THRESHOLD,
        );
    }

    reporter.report(Progress::PhaseFinish);
    Ok((forcefield, rotamer_library, active_residues))
}

#[instrument(skip_all, name = "el_energy_precomputation")]
fn precompute_el_energies<'a>(
    context: &OptimizationContext<'a, PlacementConfig>,
    reporter: &ProgressReporter,
) -> Result<ELCache, EngineError> {
    reporter.report(Progress::PhaseStart {
        name: "EL Pre-computation",
    });
    let el_cache = tasks::el_energy::run(context)?;
    reporter.report(Progress::PhaseFinish);
    Ok(el_cache)
}

#[instrument(skip_all, name = "state_initialization")]
fn initialize_state<'a>(
    initial_system: &MolecularSystem,
    active_residues: &HashSet<ResidueId>,
    context: &OptimizationContext<'a, PlacementConfig>,
    el_cache: &ELCache,
    reporter: &ProgressReporter,
) -> Result<OptimizationState, EngineError> {
    reporter.report(Progress::PhaseStart {
        name: "Initialization",
    });
    info!("Initializing system with ground-state rotamers.");

    let mut initial_rotamers = HashMap::new();
    let mut working_system = initial_system.clone();

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
            let p_info = context
                .rotamer_library
                .get_placement_info_for(residue_type)
                .unwrap();
            place_rotamer_on_system(&mut working_system, residue_id, rotamer, p_info)?;
            initial_rotamers.insert(residue_id, ground_state_idx);
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

    reporter.report(Progress::PhaseFinish);
    Ok(state)
}

#[instrument(skip_all, name = "clash_resolution_loop")]
fn resolve_clashes_iteratively<'a>(
    state: &mut OptimizationState,
    active_residues: &HashSet<ResidueId>,
    context: &OptimizationContext<'a, PlacementConfig>,
    el_cache: &ELCache,
    reporter: &ProgressReporter,
) -> Result<(), EngineError> {
    reporter.report(Progress::PhaseStart {
        name: "Clash Resolution",
    });
    info!("Starting iterative clash resolution loop.");

    let convergence_energy_threshold = context.config.optimization.convergence.energy_threshold;
    let convergence_patience_iterations =
        context.config.optimization.convergence.patience_iterations;

    let mut last_total_energy = state.best_solution().map(|s| s.energy).unwrap_or(f64::MAX);
    let mut iterations_without_significant_improvement = 0;

    for iteration in 0..context.config.optimization.max_iterations {
        let clashes = tasks::clash_detection::run(
            &state.working_state.system,
            context.forcefield,
            active_residues,
            CLASH_THRESHOLD,
            reporter,
        )?;

        if clashes.is_empty() {
            info!(
                "System converged after {} clash resolution iterations (no clashes found).",
                iteration
            );
            reporter.report(Progress::Message(format!(
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
        let p_info_a = context
            .rotamer_library
            .get_placement_info_for(res_a_type)
            .unwrap();
        place_rotamer_on_system(
            &mut state.working_state.system,
            res_a_id,
            rotamer_a,
            p_info_a,
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
        let p_info_b = context
            .rotamer_library
            .get_placement_info_for(res_b_type)
            .unwrap();
        place_rotamer_on_system(
            &mut state.working_state.system,
            res_b_id,
            rotamer_b,
            p_info_b,
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
            reporter.report(Progress::Message(format!(
                "Converged after {} iterations (energy stabilized).",
                iteration
            )));
            break;
        }
    }
    reporter.report(Progress::PhaseFinish);
    Ok(())
}

#[instrument(skip_all, name = "final_refinement_loop")]
fn final_refinement<'a>(
    state: &mut OptimizationState,
    active_residues: &HashSet<ResidueId>,
    context: &OptimizationContext<'a, PlacementConfig>,
    el_cache: &ELCache,
    reporter: &ProgressReporter,
) -> Result<(), EngineError> {
    reporter.report(Progress::PhaseStart {
        name: "Final Refinement",
    });
    info!("Starting final refinement (singlet optimization).");

    let final_refinement_iterations = context.config.optimization.final_refinement_iterations;

    for i in 0..final_refinement_iterations {
        let mut changed_in_cycle = false;

        reporter.report(Progress::TaskStart {
            total_steps: active_residues.len() as u64,
        });

        let mut residues_to_process: Vec<ResidueId> = active_residues.iter().cloned().collect();
        residues_to_process.shuffle(&mut thread_rng());

        for &residue_id in &residues_to_process {
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
            let p_info = context
                .rotamer_library
                .get_placement_info_for(residue_type)
                .unwrap();

            let current_rot_idx = *state.working_state.rotamers.get(&residue_id).unwrap();
            let mut best_rotamer_idx = current_rot_idx;
            let mut best_energy = state.current_energy;

            let mut temp_system_for_eval = state.working_state.system.clone();
            let mut temp_rotamers_for_eval = state.working_state.rotamers.clone();

            for (idx, rotamer) in rotamers.iter().enumerate() {
                if idx == current_rot_idx {
                    continue;
                }

                place_rotamer_on_system(&mut temp_system_for_eval, residue_id, rotamer, p_info)?;
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
                    p_info,
                )?;
                temp_rotamers_for_eval.insert(residue_id, current_rot_idx);
            }

            if best_rotamer_idx != current_rot_idx {
                changed_in_cycle = true;
                let best_rotamer = &rotamers[best_rotamer_idx];
                place_rotamer_on_system(
                    &mut state.working_state.system,
                    residue_id,
                    best_rotamer,
                    p_info,
                )?;
                state
                    .working_state
                    .rotamers
                    .insert(residue_id, best_rotamer_idx);
                state.current_energy = best_energy;
                state.submit_current_solution();
            }
            reporter.report(Progress::TaskIncrement);
        }

        reporter.report(Progress::TaskFinish);

        if !changed_in_cycle {
            info!(
                "Refinement converged after {} cycle(s). No further improvements.",
                i + 1
            );
            break;
        }
    }

    state.submit_current_solution();

    reporter.report(Progress::PhaseFinish);
    Ok(())
}

#[instrument(skip_all, name = "simulated_annealing_loop")]
fn run_simulated_annealing<'a>(
    state: &mut OptimizationState,
    active_residues: &HashSet<ResidueId>,
    context: &OptimizationContext<'a, PlacementConfig>,
    el_cache: &ELCache,
    reporter: &ProgressReporter,
) -> Result<(), EngineError> {
    reporter.report(Progress::PhaseStart {
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

    while current_temperature > final_temperature {
        reporter.report(Progress::Message(format!(
            "SA Temp: {:.4}",
            current_temperature
        )));
        reporter.report(Progress::TaskStart {
            total_steps: steps_per_temperature as u64,
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
            let p_info = context
                .rotamer_library
                .get_placement_info_for(residue_type)
                .unwrap();

            if rotamers.len() <= 1 {
                reporter.report(Progress::TaskIncrement);
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
                p_info,
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
            reporter.report(Progress::TaskIncrement);
        }

        reporter.report(Progress::TaskFinish);

        current_temperature *= cooling_rate;
    }
    info!("Simulated Annealing finished. Final temperature reached.");
    reporter.report(Progress::PhaseFinish);
    Ok(())
}
