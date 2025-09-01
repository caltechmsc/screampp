use super::defaults::DefaultsConfig;
use super::file::{FileConfig, FileConvergenceConfig, FileSimulatedAnnealingConfig};
use super::models::AppConfig;
use crate::cli::PlaceArgs;
use crate::data::DataManager;
use crate::error::{CliError, Result};
use crate::utils::parser;
use screampp::engine::config as core_config;
use std::path::{Path, PathBuf};

pub fn build_config(args: &PlaceArgs, data_manager: &DataManager) -> Result<AppConfig> {
    let defaults = DefaultsConfig::default();

    let file_config = if let Some(config_path) = &args.config {
        FileConfig::from_file(config_path)?
    } else {
        FileConfig::default()
    };

    let mut file_config = apply_set_values(file_config, &args.set_values)?;

    let ff_file = file_config.forcefield.take().unwrap_or_default();
    let s_factor = args
        .s_factor
        .or(ff_file.s_factor)
        .unwrap_or(defaults.s_factor);

    let sampling_file = file_config.sampling.take().unwrap_or_default();

    let opt_file = file_config.optimization.take().unwrap_or_default();
    let num_solutions = args
        .num_solutions
        .or(opt_file.num_solutions)
        .unwrap_or(defaults.num_solutions);
    let max_iterations = args
        .max_iterations
        .or(opt_file.max_iterations)
        .unwrap_or(defaults.max_iterations);

    let include_input_conformation = match (
        args.include_input_conformation.with_input_conformation,
        args.include_input_conformation.no_input_conformation,
    ) {
        (true, false) => true,
        (false, true) => false,
        _ => opt_file
            .include_input_conformation
            .unwrap_or(defaults.include_input_conformation),
    };

    let final_refinement_iterations = if args.no_refinement {
        0
    } else {
        opt_file
            .final_refinement_iterations
            .unwrap_or(defaults.final_refinement_iterations)
    };

    let forcefield_path = resolve_path_or_logical_name(
        args.forcefield_path.as_deref(),
        ff_file.forcefield_path.as_deref(),
        &defaults.forcefield,
        "forcefield",
        data_manager,
    )?;
    let delta_params_path = resolve_path_or_logical_name(
        args.delta_params_path.as_deref(),
        ff_file.delta_params_path.as_deref(),
        &defaults.delta_params,
        "delta-params",
        data_manager,
    )?;
    let rotamer_library_path = resolve_path_or_logical_name(
        args.rotamer_library.as_deref(),
        sampling_file.rotamer_library.as_deref(),
        &defaults.rotamer_library,
        "rotamer-library",
        data_manager,
    )?;
    let topology_registry_path = resolve_path_or_logical_name(
        args.topology_registry.as_deref(),
        file_config.topology_registry.as_deref(),
        &defaults.topology_registry,
        "topology-registry",
        data_manager,
    )?;

    let convergence_config = merge_convergence(opt_file.convergence, &defaults)?;
    let sa_config = merge_simulated_annealing(args.no_annealing, opt_file.simulated_annealing)?;

    let residues_to_optimize = file_config
        .residues_to_optimize
        .map(Into::into)
        .unwrap_or(core_config::ResidueSelection::All);
    let energy_weights = ff_file.energy_weights.into();

    let core_config = core_config::PlacementConfigBuilder::new()
        .forcefield_path(forcefield_path)
        .delta_params_path(delta_params_path)
        .s_factor(s_factor)
        .energy_weights(energy_weights)
        .rotamer_library_path(rotamer_library_path)
        .topology_registry_path(topology_registry_path)
        .max_iterations(max_iterations)
        .num_solutions(num_solutions)
        .include_input_conformation(include_input_conformation)
        .convergence_config(convergence_config)
        .simulated_annealing_config(sa_config)
        .final_refinement_iterations(final_refinement_iterations)
        .residues_to_optimize(residues_to_optimize)
        .build()
        .map_err(|e| CliError::Config(e.to_string()))?;

    Ok(AppConfig {
        input_path: args.input.clone(),
        output_template: args.output.clone(),
        core_config,
    })
}
