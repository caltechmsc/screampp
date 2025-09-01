use crate::error::{CliError, Result};
use screampp::core::models::atom::AtomRole;
use screampp::engine::config as core_config;
use serde::Deserialize;
use std::path::Path;
use tracing::debug;

#[derive(Deserialize, Debug, Default, Clone)]
#[serde(deny_unknown_fields)]
pub struct FileResidueSpecifier {
    #[serde(rename = "chain-id")]
    pub chain_id: char,
    #[serde(rename = "residue-number")]
    pub residue_number: isize,
}

impl From<FileResidueSpecifier> for core_config::ResidueSpecifier {
    fn from(p: FileResidueSpecifier) -> Self {
        Self {
            chain_id: p.chain_id,
            residue_number: p.residue_number,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
#[serde(rename_all = "kebab-case", tag = "type")]
pub enum FileResidueSelection {
    All,
    List {
        #[serde(default)]
        include: Vec<FileResidueSpecifier>,
        #[serde(default)]
        exclude: Vec<FileResidueSpecifier>,
    },
    LigandBindingSite {
        #[serde(rename = "ligand-residue")]
        ligand_residue: FileResidueSpecifier,
        #[serde(rename = "radius-angstroms")]
        radius_angstroms: f64,
    },
}

impl From<FileResidueSelection> for core_config::ResidueSelection {
    fn from(p: FileResidueSelection) -> Self {
        match p {
            FileResidueSelection::All => core_config::ResidueSelection::All,
            FileResidueSelection::List { include, exclude } => {
                core_config::ResidueSelection::List {
                    include: include.into_iter().map(Into::into).collect(),
                    exclude: exclude.into_iter().map(Into::into).collect(),
                }
            }
            FileResidueSelection::LigandBindingSite {
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
pub struct FileForcefieldConfig {
    #[serde(rename = "forcefield-path")]
    pub forcefield_path: Option<String>,
    #[serde(rename = "delta-params-path")]
    pub delta_params_path: Option<String>,
    #[serde(rename = "s-factor")]
    pub s_factor: Option<f64>,
    #[serde(default, rename = "energy-weights")]
    pub energy_weights: FileEnergyWeights,
}

#[derive(Deserialize, Debug, Clone, Copy)]
#[serde(deny_unknown_fields)]
pub struct FileEnergyComponentWeights {
    #[serde(default = "default_one")]
    pub vdw: f64,
    #[serde(default = "default_one")]
    pub coulomb: f64,
    #[serde(default = "default_one")]
    pub hbond: f64,
}

fn default_one() -> f64 {
    1.0
}

#[derive(Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct FileWeightRule {
    pub groups: [String; 2],
    pub weights: FileEnergyComponentWeights,
}

#[derive(Deserialize, Debug, Clone, Default)]
#[serde(deny_unknown_fields)]
pub struct FileEnergyWeights {
    #[serde(default)]
    pub rules: Vec<FileWeightRule>,
}

impl From<FileEnergyComponentWeights>
    for screampp::core::forcefield::params::EnergyComponentWeights
{
    fn from(w: FileEnergyComponentWeights) -> Self {
        Self {
            vdw: w.vdw,
            coulomb: w.coulomb,
            hbond: w.hbond,
        }
    }
}

impl From<FileEnergyWeights> for screampp::core::forcefield::params::EnergyWeights {
    fn from(w: FileEnergyWeights) -> Self {
        use std::str::FromStr;
        let mut rules: Vec<screampp::core::forcefield::params::WeightRule> = Vec::new();
        for r in w.rules.into_iter() {
            if let (Ok(g1), Ok(g2)) = (
                AtomRole::from_str(&r.groups[0]),
                AtomRole::from_str(&r.groups[1]),
            ) {
                rules.push(screampp::core::forcefield::params::WeightRule {
                    groups: [g1, g2],
                    weights: r.weights.into(),
                });
            }
        }
        screampp::core::forcefield::params::EnergyWeights { rules }
    }
}

#[derive(Deserialize, Debug, Default)]
#[serde(deny_unknown_fields)]
pub struct FileSamplingConfig {
    #[serde(rename = "rotamer-library")]
    pub rotamer_library: Option<String>,
}

#[derive(Deserialize, Debug, Default)]
#[serde(deny_unknown_fields)]
pub struct FileConvergenceConfig {
    #[serde(rename = "energy-threshold")]
    pub energy_threshold: Option<f64>,
    #[serde(rename = "patience-iterations")]
    pub patience_iterations: Option<usize>,
}

#[derive(Deserialize, Debug, Default, Clone)]
#[serde(deny_unknown_fields)]
pub struct FileSimulatedAnnealingConfig {
    #[serde(rename = "initial-temperature")]
    pub initial_temperature: Option<f64>,
    #[serde(rename = "final-temperature")]
    pub final_temperature: Option<f64>,
    #[serde(rename = "cooling-rate")]
    pub cooling_rate: Option<f64>,
    #[serde(rename = "steps-per-temperature")]
    pub steps_per_temperature: Option<usize>,
}

#[derive(Deserialize, Debug, Default)]
#[serde(deny_unknown_fields)]
pub struct FileOptimizationConfig {
    #[serde(rename = "max-iterations")]
    pub max_iterations: Option<usize>,
    #[serde(rename = "num-solutions")]
    pub num_solutions: Option<usize>,
    #[serde(rename = "include-input-conformation")]
    pub include_input_conformation: Option<bool>,
    pub convergence: Option<FileConvergenceConfig>,
    #[serde(rename = "simulated-annealing")]
    pub simulated_annealing: Option<FileSimulatedAnnealingConfig>,
    #[serde(rename = "final-refinement-iterations")]
    pub final_refinement_iterations: Option<usize>,
}

#[derive(Deserialize, Debug, Default)]
#[serde(deny_unknown_fields)]
pub struct FileConfig {
    pub forcefield: Option<FileForcefieldConfig>,
    pub sampling: Option<FileSamplingConfig>,
    pub optimization: Option<FileOptimizationConfig>,
    #[serde(rename = "residues-to-optimize")]
    pub residues_to_optimize: Option<FileResidueSelection>,
    #[serde(rename = "topology-registry-path")]
    pub topology_registry: Option<String>,
}

impl FileConfig {
    pub fn from_file(path: &Path) -> Result<Self> {
        debug!("Loading configuration from file: {:?}", path);
        let content = std::fs::read_to_string(path)?;
        toml::from_str(&content).map_err(|e| CliError::FileParsing {
            path: path.to_path_buf(),
            source: e.into(),
        })
    }
}
