use crate::cli::{IncludeInputConformation, PlaceArgs};
use crate::data::DataManager;
use crate::error::{CliError, Result};
use crate::utils::parser;
use screampp::engine::config as core_config;
use serde::Deserialize;
use std::path::{Path, PathBuf};
use tracing::debug;

#[derive(Deserialize, Debug, Default, Clone)]
#[serde(deny_unknown_fields)]
struct PartialResidueSpecifier {
    #[serde(rename = "chain-id")]
    chain_id: char,
    #[serde(rename = "residue-number")]
    residue_number: isize,
}

impl From<PartialResidueSpecifier> for core_config::ResidueSpecifier {
    fn from(p: PartialResidueSpecifier) -> Self {
        Self {
            chain_id: p.chain_id,
            residue_number: p.residue_number,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
#[serde(rename_all = "kebab-case", tag = "type")]
enum PartialResidueSelection {
    All,
    List {
        include: Vec<PartialResidueSpecifier>,
        exclude: Vec<PartialResidueSpecifier>,
    },
    LigandBindingSite {
        ligand_residue: PartialResidueSpecifier,
        radius_angstroms: f64,
    },
}

impl From<PartialResidueSelection> for core_config::ResidueSelection {
    fn from(p: PartialResidueSelection) -> Self {
        match p {
            PartialResidueSelection::All => core_config::ResidueSelection::All,
            PartialResidueSelection::List { include, exclude } => {
                core_config::ResidueSelection::List {
                    include: include.into_iter().map(Into::into).collect(),
                    exclude: exclude.into_iter().map(Into::into).collect(),
                }
            }
            PartialResidueSelection::LigandBindingSite {
                ligand_residue,
                radius_angstroms,
            } => core_config::ResidueSelection::LigandBindingSite {
                ligand_residue: ligand_residue.into(),
                radius_angstroms,
            },
        }
    }
}

#[derive(Deserialize, Debug, Default)]
#[serde(deny_unknown_fields)]
struct PartialForcefieldConfig {
    #[serde(rename = "forcefield-path")]
    forcefield_path: Option<String>,
    #[serde(rename = "delta-params-path")]
    delta_params_path: Option<String>,
    #[serde(rename = "s-factor")]
    s_factor: Option<f64>,
}

#[derive(Deserialize, Debug, Default)]
#[serde(deny_unknown_fields)]
struct PartialSamplingConfig {
    #[serde(rename = "rotamer-library")]
    rotamer_library: Option<String>,
    #[serde(rename = "placement-registry")]
    placement_registry: Option<String>,
}

#[derive(Deserialize, Debug, Default)]
#[serde(deny_unknown_fields)]
struct PartialConvergenceConfig {
    #[serde(rename = "energy-threshold")]
    energy_threshold: Option<f64>,
    #[serde(rename = "patience-iterations")]
    patience_iterations: Option<usize>,
}

#[derive(Deserialize, Debug, Default, Clone)]
#[serde(deny_unknown_fields)]
struct PartialSimulatedAnnealingConfig {
    #[serde(rename = "initial-temperature")]
    initial_temperature: Option<f64>,
    #[serde(rename = "final-temperature")]
    final_temperature: Option<f64>,
    #[serde(rename = "cooling-rate")]
    cooling_rate: Option<f64>,
    #[serde(rename = "steps-per-temperature")]
    steps_per_temperature: Option<usize>,
}

#[derive(Deserialize, Debug, Default)]
#[serde(deny_unknown_fields)]
struct PartialOptimizationConfig {
    #[serde(rename = "max-iterations")]
    max_iterations: Option<usize>,
    #[serde(rename = "num-solutions")]
    num_solutions: Option<usize>,
    #[serde(rename = "include-input-conformation")]
    include_input_conformation: Option<bool>,
    convergence: Option<PartialConvergenceConfig>,
    #[serde(rename = "simulated-annealing")]
    simulated_annealing: Option<PartialSimulatedAnnealingConfig>,
    #[serde(rename = "final-refinement-iterations")]
    final_refinement_iterations: Option<usize>,
}

#[derive(Deserialize, Debug, Default)]
#[serde(deny_unknown_fields)]
pub struct PartialPlacementConfig {
    forcefield: Option<PartialForcefieldConfig>,
    sampling: Option<PartialSamplingConfig>,
    optimization: Option<PartialOptimizationConfig>,
    #[serde(rename = "residues-to-optimize")]
    residues_to_optimize: Option<PartialResidueSelection>,
}

impl PartialPlacementConfig {
    pub fn from_file(path: &Path) -> Result<Self> {
        debug!("Loading configuration from file: {:?}", path);
        let content = std::fs::read_to_string(path)?;
        toml::from_str(&content).map_err(|e| CliError::FileParsing {
            path: path.to_path_buf(),
            source: e.into(),
        })
    }

    pub fn merge_with_cli(
        mut self,
        args: &PlaceArgs,
        data_manager: &DataManager,
    ) -> Result<core_config::PlacementConfig> {
        self.apply_set_values(&args.set_values)?;

        let ff_config = self.forcefield.take().unwrap_or_default();
        let sampling_config = self.sampling.take().unwrap_or_default();
        let opt_config = self.optimization.take().unwrap_or_default();

        let resolve = |name_or_path_opt: Option<&String>, kind: &str| -> Result<PathBuf> {
            let name_or_path = name_or_path_opt.ok_or_else(|| {
                CliError::Config(format!(
                    "A value for '{}' is required either in the config file or via CLI argument.",
                    kind
                ))
            })?;

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
        };

        let forcefield_path = resolve(
            args.forcefield_path
                .as_ref()
                .or(ff_config.forcefield_path.as_ref()),
            "forcefield",
        )?;
        let delta_params_path = resolve(
            args.delta_params_path
                .as_ref()
                .or(ff_config.delta_params_path.as_ref()),
            "delta-params",
        )?;
        let rotamer_library_path = resolve(
            args.rotamer_library
                .as_ref()
                .or(sampling_config.rotamer_library.as_ref()),
            "rotamer-library",
        )?;
        let placement_registry_path = resolve(
            args.placement_registry
                .as_ref()
                .or(sampling_config.placement_registry.as_ref()),
            "placement-registry",
        )?;
        let s_factor = args
            .s_factor
            .or(ff_config.s_factor)
            .ok_or_else(|| CliError::Config("`forcefield.s-factor` is required.".to_string()))?;

        let final_residue_selection = self.residues_to_optimize.ok_or_else(|| {
            CliError::Config("`residues-to-optimize` section is required.".to_string())
        })?;

        let mut builder = core_config::PlacementConfigBuilder::new()
            .forcefield_path(forcefield_path)
            .delta_params_path(delta_params_path)
            .s_factor(s_factor)
            .rotamer_library_path(rotamer_library_path)
            .placement_registry_path(placement_registry_path)
            .max_iterations(
                args.max_iterations
                    .or(opt_config.max_iterations)
                    .unwrap_or(100),
            )
            .num_solutions(args.num_solutions.or(opt_config.num_solutions).unwrap_or(1))
            .residues_to_optimize(final_residue_selection.into());

        builder = Self::merge_include_conformation(
            builder,
            args.include_input_conformation,
            opt_config.include_input_conformation,
        );
        builder = Self::merge_final_refinement(
            builder,
            args.no_refinement,
            opt_config.final_refinement_iterations,
        );
        builder = Self::merge_simulated_annealing(
            builder,
            args.no_annealing,
            opt_config.simulated_annealing,
        )?;

        let convergence = Self::merge_convergence(opt_config.convergence)?;
        builder = builder.convergence_config(convergence);

        builder.build().map_err(|e| CliError::Config(e.to_string()))
    }

    fn merge_include_conformation(
        builder: core_config::PlacementConfigBuilder,
        cli_flags: IncludeInputConformation,
        file_val: Option<bool>,
    ) -> core_config::PlacementConfigBuilder {
        if cli_flags.with_input_conformation {
            builder.include_input_conformation(true)
        } else if cli_flags.no_input_conformation {
            builder.include_input_conformation(false)
        } else if let Some(val) = file_val {
            builder.include_input_conformation(val)
        } else {
            builder.include_input_conformation(false)
        }
    }

    fn merge_final_refinement(
        builder: core_config::PlacementConfigBuilder,
        cli_no_refinement: bool,
        file_val: Option<usize>,
    ) -> core_config::PlacementConfigBuilder {
        if cli_no_refinement {
            builder.final_refinement_iterations(0)
        } else if let Some(val) = file_val {
            builder.final_refinement_iterations(val)
        } else {
            builder.final_refinement_iterations(2)
        }
    }

    fn merge_convergence(
        partial: Option<PartialConvergenceConfig>,
    ) -> Result<core_config::ConvergenceConfig> {
        let partial = partial.unwrap_or_default();
        Ok(core_config::ConvergenceConfig {
            energy_threshold: partial.energy_threshold.unwrap_or(0.01),
            patience_iterations: partial.patience_iterations.unwrap_or(5),
        })
    }

    fn merge_simulated_annealing(
        builder: core_config::PlacementConfigBuilder,
        cli_no_annealing: bool,
        partial: Option<PartialSimulatedAnnealingConfig>,
    ) -> Result<core_config::PlacementConfigBuilder> {
        if cli_no_annealing {
            return Ok(builder.simulated_annealing_config(None));
        }
        if let Some(p) = partial {
            let config = core_config::SimulatedAnnealingConfig {
                initial_temperature: p.initial_temperature.ok_or_else(|| {
                    CliError::Config(
                        "`simulated-annealing` requires `initial-temperature`".to_string(),
                    )
                })?,
                final_temperature: p.final_temperature.ok_or_else(|| {
                    CliError::Config(
                        "`simulated-annealing` requires `final-temperature`".to_string(),
                    )
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
            return Ok(builder.simulated_annealing_config(Some(config)));
        }
        Ok(builder.simulated_annealing_config(None))
    }

    fn apply_set_values(&mut self, set_values: &[String]) -> Result<()> {
        if set_values.is_empty() {
            return Ok(());
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
                    self.forcefield
                        .get_or_insert_with(Default::default)
                        .s_factor = Some(value_str.parse().map_err(|_| {
                        CliError::Config(format!("Invalid float value for {}: {}", key, value_str))
                    })?);
                }
                "optimization.max-iterations" => {
                    self.optimization
                        .get_or_insert_with(Default::default)
                        .max_iterations = Some(value_str.parse().map_err(|_| {
                        CliError::Config(format!(
                            "Invalid integer value for {}: {}",
                            key, value_str
                        ))
                    })?);
                }
                "optimization.num-solutions" => {
                    self.optimization
                        .get_or_insert_with(Default::default)
                        .num_solutions = Some(value_str.parse().map_err(|_| {
                        CliError::Config(format!(
                            "Invalid integer value for {}: {}",
                            key, value_str
                        ))
                    })?);
                }
                "optimization.convergence.energy-threshold" => {
                    self.optimization
                        .get_or_insert_with(Default::default)
                        .convergence
                        .get_or_insert_with(Default::default)
                        .energy_threshold = Some(value_str.parse().map_err(|_| {
                        CliError::Config(format!("Invalid float value for {}: {}", key, value_str))
                    })?);
                }
                "optimization.convergence.patience-iterations" => {
                    self.optimization
                        .get_or_insert_with(Default::default)
                        .convergence
                        .get_or_insert_with(Default::default)
                        .patience_iterations = Some(value_str.parse().map_err(|_| {
                        CliError::Config(format!(
                            "Invalid integer value for {}: {}",
                            key, value_str
                        ))
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
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::Cli;
    use clap::Parser;
    use once_cell::sync::Lazy;
    use std::fs;
    use tempfile::{TempDir, tempdir};

    static TEST_DIR: Lazy<TempDir> = Lazy::new(|| tempdir().expect("Failed to create temp dir"));

    fn create_test_data_dir() -> PathBuf {
        let data_path = TEST_DIR.path().join("mock_data");
        let nested_data_path = data_path.join("data");
        fs::create_dir_all(nested_data_path.join("rotamers/amber")).unwrap();
        fs::write(nested_data_path.join("rotamers/amber/rmsd-1.0.toml"), "").unwrap();
        fs::create_dir_all(nested_data_path.join("forcefield")).unwrap();
        fs::write(
            nested_data_path.join("forcefield/dreiding-lj-12-6-0.4.toml"),
            "",
        )
        .unwrap();
        data_path
    }

    fn write_config_file(name: &str, content: &str) -> PathBuf {
        let file_path = TEST_DIR.path().join(name);
        fs::write(&file_path, content).unwrap();
        file_path
    }

    fn get_minimal_cli_args(config_path: &Path) -> Vec<String> {
        vec![
            "scream".to_string(),
            "place".to_string(),
            "-i".to_string(),
            "in.bgf".to_string(),
            "-o".to_string(),
            "out.bgf".to_string(),
            "-c".to_string(),
            config_path.to_str().unwrap().to_string(),
        ]
    }

    #[test]
    fn test_load_from_file_and_merge_with_defaults() {
        let data_path = create_test_data_dir();

        let manager = DataManager::with_custom_path(data_path);

        let delta_csv_path = TEST_DIR.path().join("delta.csv");
        let reg_toml_path = TEST_DIR.path().join("reg.toml");
        fs::write(&delta_csv_path, "").unwrap();
        fs::write(&reg_toml_path, "").unwrap();

        let config_content = format!(
            r#"
        [forcefield]
        s-factor = 0.8
        forcefield-path = "lj-12-6@0.4"
        delta-params-path = "{}"

        [sampling]
        rotamer-library = "amber@rmsd-1.0"
        placement-registry = "{}"

        [residues-to-optimize]
        type = "all"
        "#,
            delta_csv_path.to_str().unwrap(),
            reg_toml_path.to_str().unwrap()
        );

        let config_path = write_config_file("config_defaults.toml", &config_content);
        let args = get_minimal_cli_args(&config_path);
        let cli = Cli::parse_from(args);

        if let crate::cli::Commands::Place(place_args) = cli.command {
            let partial_config = PartialPlacementConfig::from_file(&config_path).unwrap();
            let final_config = partial_config
                .merge_with_cli(&place_args, &manager)
                .unwrap();

            assert_eq!(final_config.forcefield.s_factor, 0.8);
            assert_eq!(
                final_config.sampling.rotamer_library_path,
                manager
                    .get_data_path()
                    .join("data")
                    .join("rotamers/amber/rmsd-1.0.toml")
            );
            assert_eq!(
                final_config.forcefield.forcefield_path,
                manager
                    .get_data_path()
                    .join("data")
                    .join("forcefield/dreiding-lj-12-6-0.4.toml")
            );

            assert_eq!(final_config.optimization.num_solutions, 1);
            assert_eq!(final_config.optimization.max_iterations, 100);
            assert_eq!(final_config.optimization.final_refinement_iterations, 2);
            assert!(!final_config.optimization.include_input_conformation);
        } else {
            panic!("Expected 'place' subcommand");
        }
    }

    #[test]
    fn test_cli_args_override_file_values() {
        let data_path = create_test_data_dir();

        let manager = DataManager::with_custom_path(data_path);

        let delta_csv_path = TEST_DIR.path().join("delta.csv");
        let reg_toml_path = TEST_DIR.path().join("reg.toml");
        fs::write(&delta_csv_path, "").unwrap();
        fs::write(&reg_toml_path, "").unwrap();

        let config_content = format!(
            r#"
        [forcefield]
        s-factor = 0.5 # Will be overridden
        forcefield-path = "lj-12-6@0.4"
        delta-params-path = "{}"

        [sampling]
        rotamer-library = "amber@rmsd-1.0"
        placement-registry = "{}"

        [optimization]
        num-solutions = 5 # Will be overridden

        [residues-to-optimize]
        type = "all"
        "#,
            delta_csv_path.to_str().unwrap(),
            reg_toml_path.to_str().unwrap()
        );

        let config_path = write_config_file("config_override.toml", &config_content);
        let mut args = get_minimal_cli_args(&config_path);
        args.extend_from_slice(&[
            "--s-factor".to_string(),
            "0.9".to_string(),
            "--num-solutions".to_string(),
            "10".to_string(),
            "--with-input-conformation".to_string(),
            "--rotamer-library".to_string(),
            "amber@rmsd-1.0".to_string(),
        ]);
        let cli = Cli::parse_from(args);

        if let crate::cli::Commands::Place(place_args) = cli.command {
            let partial_config = PartialPlacementConfig::from_file(&config_path).unwrap();
            let final_config = partial_config
                .merge_with_cli(&place_args, &manager)
                .unwrap();

            assert_eq!(final_config.forcefield.s_factor, 0.9);
            assert_eq!(final_config.optimization.num_solutions, 10);
            assert!(final_config.optimization.include_input_conformation);
        } else {
            panic!("Expected 'place' subcommand");
        }
    }

    #[test]
    fn test_set_value_overrides_file_and_defaults() {
        let data_path = create_test_data_dir();

        let manager = DataManager::with_custom_path(data_path);

        let delta_csv_path = TEST_DIR.path().join("delta.csv");
        let reg_toml_path = TEST_DIR.path().join("reg.toml");
        fs::write(&delta_csv_path, "").unwrap();
        fs::write(&reg_toml_path, "").unwrap();

        let config_content = format!(
            r#"
        [forcefield]
        s-factor = 0.8
        forcefield-path = "lj-12-6@0.4"
        delta-params-path = "{}"

        [sampling]
        rotamer-library = "amber@rmsd-1.0"
        placement-registry = "{}"

        [optimization]
        num-solutions = 5 # Will be overridden by --set

        [residues-to-optimize]
        type = "all"
        "#,
            delta_csv_path.to_str().unwrap(),
            reg_toml_path.to_str().unwrap()
        );

        let config_path = write_config_file("config_set.toml", &config_content);
        let mut args = get_minimal_cli_args(&config_path);
        args.extend_from_slice(&[
            "-S".to_string(),
            "optimization.num-solutions=20".to_string(),
            "-S".to_string(),
            "optimization.convergence.energy-threshold=0.005".to_string(),
        ]);
        let cli = Cli::parse_from(args);

        if let crate::cli::Commands::Place(place_args) = cli.command {
            let partial_config = PartialPlacementConfig::from_file(&config_path).unwrap();
            let final_config = partial_config
                .merge_with_cli(&place_args, &manager)
                .unwrap();

            assert_eq!(final_config.optimization.num_solutions, 20);
            assert_eq!(
                final_config.optimization.convergence.energy_threshold,
                0.005
            );
            assert_eq!(final_config.optimization.convergence.patience_iterations, 5);
        } else {
            panic!("Expected 'place' subcommand");
        }
    }

    #[test]
    fn test_missing_required_field_returns_error() {
        let data_path = create_test_data_dir();

        let manager = DataManager::with_custom_path(data_path);

        let reg_toml_path = TEST_DIR.path().join("reg.toml");
        fs::write(&reg_toml_path, "").unwrap();

        let config_content = format!(
            r#"
        [sampling]
        rotamer-library = "amber@rmsd-1.0"
        placement-registry = "{}"
        [residues-to-optimize]
        type = "all"
        "#,
            reg_toml_path.to_str().unwrap()
        );

        let config_path = write_config_file("config_missing.toml", &config_content);
        let args = get_minimal_cli_args(&config_path);
        let cli = Cli::parse_from(args);

        if let crate::cli::Commands::Place(place_args) = cli.command {
            let partial_config = PartialPlacementConfig::from_file(&config_path).unwrap();
            let result = partial_config.merge_with_cli(&place_args, &manager);
            assert!(matches!(result, Err(CliError::Config(_))));
            if let Err(CliError::Config(msg)) = result {
                assert!(msg.contains("forcefield") || msg.contains("s-factor"));
            }
        } else {
            panic!("Expected 'place' subcommand");
        }
    }
}
