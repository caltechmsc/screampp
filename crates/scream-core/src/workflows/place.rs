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
