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
