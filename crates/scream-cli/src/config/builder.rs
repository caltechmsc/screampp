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

fn resolve_path_or_logical_name(
    cli_arg: Option<&str>,
    file_arg: Option<&str>,
    default_arg: &str,
    kind: &str,
    data_manager: &DataManager,
) -> Result<PathBuf> {
    let name_or_path = cli_arg.or(file_arg).unwrap_or(default_arg);

    let path = Path::new(name_or_path);
    if path.is_absolute() || name_or_path.contains(['/', '\\']) {
        if !path.exists() {
            return Err(CliError::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("Provided path does not exist: {}", name_or_path),
            )));
        }
        return Ok(path.to_path_buf());
    }

    let parsed_name = parser::parse_logical_name(name_or_path, kind)
        .map_err(|e| CliError::Argument(e.to_string()))?;

    let resolved = data_manager.resolve_logical_name(&parsed_name)?;
    if !resolved.exists() {
        return Err(CliError::Data(format!(
            "Resolved data file does not exist: {:?}.\nHint: Run 'scream data download' to fetch the default data files.",
            resolved
        )));
    }
    Ok(resolved)
}

fn merge_convergence(
    file_val: Option<FileConvergenceConfig>,
    defaults: &DefaultsConfig,
) -> Result<core_config::ConvergenceConfig> {
    let file_val = file_val.unwrap_or_default();
    Ok(core_config::ConvergenceConfig {
        energy_threshold: file_val
            .energy_threshold
            .unwrap_or(defaults.energy_threshold),
        patience_iterations: file_val
            .patience_iterations
            .unwrap_or(defaults.patience_iterations),
    })
}

fn merge_simulated_annealing(
    cli_no_annealing: bool,
    file_val: Option<FileSimulatedAnnealingConfig>,
) -> Result<Option<core_config::SimulatedAnnealingConfig>> {
    if cli_no_annealing {
        return Ok(None);
    }
    if let Some(p) = file_val {
        let config = core_config::SimulatedAnnealingConfig {
            initial_temperature: p.initial_temperature.ok_or_else(|| {
                CliError::Config("`simulated-annealing` requires `initial-temperature`".to_string())
            })?,
            final_temperature: p.final_temperature.ok_or_else(|| {
                CliError::Config("`simulated-annealing` requires `final-temperature`".to_string())
            })?,
            cooling_rate: p.cooling_rate.ok_or_else(|| {
                CliError::Config("`simulated-annealing` requires `cooling-rate`".to_string())
            })?,
            steps_per_temperature: p.steps_per_temperature.ok_or_else(|| {
                CliError::Config(
                    "`simulated-annealing` requires `steps-per-temperature`".to_string(),
                )
            })?,
        };
        return Ok(Some(config));
    }
    Ok(None)
}

fn apply_set_values(mut config: FileConfig, set_values: &[String]) -> Result<FileConfig> {
    if set_values.is_empty() {
        return Ok(config);
    }
    for kv_pair in set_values {
        let parts: Vec<_> = kv_pair.splitn(2, '=').collect();
        if parts.len() != 2 {
            return Err(CliError::Config(format!(
                "Invalid --set format: '{}'. Expected KEY=VALUE.",
                kv_pair
            )));
        }
        let key = parts[0];
        let value_str = parts[1];

        match key {
            "forcefield.s-factor" => {
                config
                    .forcefield
                    .get_or_insert_with(Default::default)
                    .s_factor = Some(value_str.parse().map_err(|_| {
                    CliError::Config(format!("Invalid float value for {}: {}", key, value_str))
                })?);
            }
            "optimization.max-iterations" => {
                config
                    .optimization
                    .get_or_insert_with(Default::default)
                    .max_iterations = Some(value_str.parse().map_err(|_| {
                    CliError::Config(format!("Invalid integer value for {}: {}", key, value_str))
                })?);
            }
            "optimization.num-solutions" => {
                config
                    .optimization
                    .get_or_insert_with(Default::default)
                    .num_solutions = Some(value_str.parse().map_err(|_| {
                    CliError::Config(format!("Invalid integer value for {}: {}", key, value_str))
                })?);
            }
            "optimization.convergence.energy-threshold" => {
                config
                    .optimization
                    .get_or_insert_with(Default::default)
                    .convergence
                    .get_or_insert_with(Default::default)
                    .energy_threshold = Some(value_str.parse().map_err(|_| {
                    CliError::Config(format!("Invalid float value for {}: {}", key, value_str))
                })?);
            }
            "optimization.convergence.patience-iterations" => {
                config
                    .optimization
                    .get_or_insert_with(Default::default)
                    .convergence
                    .get_or_insert_with(Default::default)
                    .patience_iterations = Some(value_str.parse().map_err(|_| {
                    CliError::Config(format!("Invalid integer value for {}: {}", key, value_str))
                })?);
            }
            _ => {
                return Err(CliError::Config(format!(
                    "Unsupported configuration key for --set: '{}'",
                    key
                )));
            }
        }
    }
    Ok(config)
}
