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
use std::collections::{HashMap, HashSet};
use tracing::{info, instrument, warn};

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

    // --- Phase 4: Final Refinement (Singlet Optimization) ---
    final_refinement(&mut state, &active_residues, &context, &el_cache, reporter)?;

    // --- Phase 5: Finalization ---
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
        rotamer_library.include_system_conformations(initial_system, &active_residues, 0.1);
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
        let res_type = working_system
            .residue(residue_id)
            .and_then(|r| r.res_type)
            .ok_or(EngineError::Internal(format!(
                "Active residue {:?} has no residue type.",
                residue_id
            )))?;

        if let Some((ground_state_idx, _)) = el_cache.find_ground_state_for(residue_id, res_type) {
            let rotamer =
                &context.rotamer_library.get_rotamers_for(res_type).unwrap()[ground_state_idx];
            let p_info = context
                .rotamer_library
                .get_placement_info_for(res_type)
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

// Replace the existing function with this one.
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

    for iteration in 0..context.config.optimization.max_iterations {
        let clashes = tasks::clash_detection::run(
            &state.working_state.system,
            context.forcefield,
            active_residues,
            25.0,
            reporter,
        )?;

        if clashes.is_empty() {
            info!(
                "System converged after {} clash resolution iterations.",
                iteration
            );
            reporter.report(Progress::Message(format!(
                "Converged after {} iterations.",
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
            .res_type
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
            .res_type
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
    }
    reporter.report(Progress::PhaseFinish);
    Ok(())
}

// TODO: Implement simulated annealing
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

    for i in 0..2 {
        let mut changed_in_cycle = false;

        reporter.report(Progress::TaskStart {
            total_steps: active_residues.len() as u64,
        });

        for &residue_id in active_residues {
            let res_type = state
                .working_state
                .system
                .residue(residue_id)
                .unwrap()
                .res_type
                .unwrap();
            let rotamers = context.rotamer_library.get_rotamers_for(res_type).unwrap();
            let p_info = context
                .rotamer_library
                .get_placement_info_for(res_type)
                .unwrap();

            let mut best_rotamer_idx = *state.working_state.rotamers.get(&residue_id).unwrap();
            let mut best_energy = state.current_energy;

            let mut temp_system = state.working_state.system.clone();

            for (idx, rotamer) in rotamers.iter().enumerate() {
                if idx == best_rotamer_idx {
                    continue;
                }

                place_rotamer_on_system(&mut temp_system, residue_id, rotamer, p_info)?;

                let mut temp_rotamers = state.working_state.rotamers.clone();
                temp_rotamers.insert(residue_id, idx);
                let current_energy = tasks::total_energy::run(
                    &temp_system,
                    context.forcefield,
                    active_residues,
                    &temp_rotamers,
                    el_cache,
                )?
                .total();

                if current_energy < best_energy {
                    best_energy = current_energy;
                    best_rotamer_idx = idx;
                }
            }

            let current_rot_idx = *state.working_state.rotamers.get(&residue_id).unwrap();
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
            }
            reporter.report(Progress::TaskIncrement);
        }

        reporter.report(Progress::TaskFinish);

        if !changed_in_cycle {
            info!("Refinement converged after {} cycle(s).", i + 1);
            break;
        }
    }

    state.submit_current_solution();

    reporter.report(Progress::PhaseFinish);
    Ok(())
}
