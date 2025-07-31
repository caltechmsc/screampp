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
