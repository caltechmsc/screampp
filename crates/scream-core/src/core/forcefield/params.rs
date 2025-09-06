use crate::core::models::atom::AtomRole;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use thiserror::Error;

/// Represents the parameters for van der Waals interactions.
///
/// This enum supports two common potential functions used in molecular mechanics:
/// Buckingham exponential-6 and Lennard-Jones 12-6 potentials.
#[derive(Debug, Deserialize, Clone, PartialEq)]
#[serde(untagged)]
pub enum VdwParam {
    /// Parameters for the Buckingham exponential-6 potential.
    ///
    /// This potential combines exponential repulsion with r⁻⁶ dispersion attraction,
    /// providing better long-range behavior than Lennard-Jones for some systems.
    Buckingham {
        /// The van der Waals radius parameter.
        radius: f64,
        /// The potential well depth.
        well_depth: f64,
        /// The exponential decay parameter.
        scale: f64,
    },
    /// Parameters for the Lennard-Jones 12-6 potential.
    ///
    /// This is the classic potential with r⁻¹² repulsion and r⁻⁶ attraction terms.
    LennardJones {
        /// The van der Waals radius parameter.
        radius: f64,
        /// The potential well depth.
        well_depth: f64,
    },
}

/// Parameters for hydrogen bond interactions.
///
/// Hydrogen bonds are directional interactions between donor and acceptor atoms,
/// modeled with a specialized potential function.
#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct HBondParam {
    /// The equilibrium distance for the hydrogen bond.
    pub equilibrium_distance: f64,
    /// The depth of the hydrogen bond potential well.
    pub well_depth: f64,
}

/// Global parameters that apply to the entire force field.
///
/// These parameters define the overall behavior of the force field calculations,
/// such as the dielectric constant and the type of potential function used.
#[derive(Debug, Deserialize, Clone, PartialEq, Default)]
pub struct GlobalParams {
    /// The dielectric constant used in electrostatic calculations.
    pub dielectric_constant: f64,
    /// The name of the potential function to use for van der Waals interactions.
    pub potential_function: String,
}

/// Parameters for non-bonded interactions in the force field.
///
/// This struct contains all the parameters needed for calculating van der Waals,
/// electrostatic, and hydrogen bond interactions between atoms.
#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct NonBondedParams {
    /// Global parameters for the force field.
    pub globals: GlobalParams,
    /// Van der Waals parameters for each atom type.
    pub vdw: HashMap<String, VdwParam>,
    /// Hydrogen bond parameters for donor-acceptor pairs.
    pub hbond: HashMap<String, HBondParam>,

    /// Set of atom types that can act as hydrogen bond donors.
    ///
    /// This field is populated automatically during loading based on the
    /// hydrogen bond parameter keys.
    #[serde(skip)]
    pub hbond_donors: HashSet<String>,
    /// Set of atom types that can act as hydrogen bond acceptors.
    ///
    /// This field is populated automatically during loading based on the
    /// hydrogen bond parameter keys.
    #[serde(skip)]
    pub hbond_acceptors: HashSet<String>,
}

/// Parameters for atom-specific corrections (deltas) in the force field.
///
/// These parameters allow for fine-tuning of atomic properties on a per-residue,
/// per-atom basis, typically used for improving agreement with experimental data.
#[derive(Debug, Deserialize, Clone)]
pub struct DeltaParam {
    /// The residue type for which this delta applies.
    pub residue_type: String,
    /// The atom name within the residue.
    pub atom_name: String,
    /// The mean correction value.
    pub mu: f64,
    /// The standard deviation of the correction.
    pub sigma: f64,
}

/// Weights for different energy components in force field calculations.
///
/// This struct allows scaling the contribution of different interaction types
/// to the total energy, which can be useful for optimization or analysis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyComponentWeights {
    /// Weight for van der Waals interactions.
    pub vdw: f64,
    /// Weight for electrostatic interactions.
    pub coulomb: f64,
    /// Weight for hydrogen bond interactions.
    pub hbond: f64,
}

impl Default for EnergyComponentWeights {
    fn default() -> Self {
        Self {
            vdw: 1.0,
            coulomb: 1.0,
            hbond: 1.0,
        }
    }
}

/// A rule for applying energy weights based on atom roles.
///
/// This struct defines how energy components should be weighted when calculating
/// interactions between atoms of specific roles (e.g., backbone vs. sidechain).
#[derive(Debug, Clone, PartialEq)]
pub struct WeightRule {
    /// The pair of atom roles to which this rule applies.
    pub groups: [AtomRole; 2],
    /// The weights to apply for interactions between these roles.
    pub weights: EnergyComponentWeights,
}

/// Configuration for energy component weights across different atom role pairs.
///
/// This struct contains a collection of rules that define how energy components
/// should be weighted based on the roles of the interacting atoms.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct EnergyWeights {
    /// The list of weighting rules to apply.
    pub rules: Vec<WeightRule>,
}

/// The complete force field parameter set.
///
/// This struct encapsulates all parameters needed for molecular mechanics
/// calculations, including non-bonded interactions, atomic corrections,
/// and energy weighting rules.
#[derive(Debug, Clone)]
pub struct Forcefield {
    /// Parameters for non-bonded interactions.
    pub non_bonded: NonBondedParams,
    /// Atomic correction parameters indexed by (residue_type, atom_name).
    pub deltas: HashMap<(String, String), DeltaParam>,
    /// Energy weights for different atom role pairs.
    pub weight_map: HashMap<(AtomRole, AtomRole), EnergyComponentWeights>,
}

/// Errors that can occur during parameter loading.
///
/// This enum covers various failure modes when loading force field parameters
/// from configuration files.
#[derive(Debug, Error)]
pub enum ParamLoadError {
    /// An I/O error occurred while reading a file.
    #[error("File I/O error for '{path}': {source}")]
    Io {
        /// The path to the file that caused the error.
        path: String,
        /// The underlying I/O error.
        source: std::io::Error,
    },
    /// A CSV parsing error occurred.
    #[error("CSV parsing error for '{path}': {source}")]
    Csv {
        /// The path to the CSV file that caused the error.
        path: String,
        /// The underlying CSV parsing error.
        source: csv::Error,
    },
    /// A TOML parsing error occurred.
    #[error("TOML parsing error for '{path}': {source}")]
    Toml {
        /// The path to the TOML file that caused the error.
        path: String,
        /// The underlying TOML parsing error.
        source: toml::de::Error,
    },
}

impl Forcefield {
    /// Loads a complete force field from configuration files.
    ///
    /// This method reads non-bonded parameters from a TOML file, delta parameters
    /// from a CSV file, and applies energy weighting rules to create a complete
    /// force field parameter set.
    ///
    /// # Arguments
    ///
    /// * `non_bonded_path` - Path to the TOML file containing non-bonded parameters.
    /// * `delta_path` - Path to the CSV file containing delta parameters.
    /// * `energy_weights_config` - Configuration for energy component weights.
    ///
    /// # Return
    ///
    /// Returns a `Forcefield` instance with all parameters loaded.
    ///
    /// # Errors
    ///
    /// Returns a `ParamLoadError` if any of the files cannot be read or parsed.
    pub fn load(
        non_bonded_path: &Path,
        delta_path: &Path,
        energy_weights_config: &EnergyWeights,
    ) -> Result<Self, ParamLoadError> {
        let non_bonded = Self::load_non_bonded(non_bonded_path)?;
        let deltas = Self::load_delta_csv(delta_path)?;

        let mut weight_map: HashMap<(AtomRole, AtomRole), EnergyComponentWeights> = HashMap::new();
        for rule in &energy_weights_config.rules {
            let a = rule.groups[0];
            let b = rule.groups[1];
            // Ensure consistent ordering for the map key
            let key = if a <= b { (a, b) } else { (b, a) };
            weight_map.insert(key, rule.weights);
        }

        Ok(Self {
            non_bonded,
            deltas,
            weight_map,
        })
    }

    /// Loads non-bonded parameters from a TOML file.
    ///
    /// This method parses the TOML file and automatically populates the
    /// hydrogen bond donor and acceptor sets based on the parameter keys.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the TOML file containing non-bonded parameters.
    ///
    /// # Return
    ///
    /// Returns the parsed `NonBondedParams`.
    ///
    /// # Errors
    ///
    /// Returns a `ParamLoadError` if the file cannot be read or parsed.
    fn load_non_bonded(path: &Path) -> Result<NonBondedParams, ParamLoadError> {
        let content = std::fs::read_to_string(path).map_err(|e| ParamLoadError::Io {
            path: path.to_string_lossy().to_string(),
            source: e,
        })?;

        let mut params: NonBondedParams =
            toml::from_str(&content).map_err(|e| ParamLoadError::Toml {
                path: path.to_string_lossy().to_string(),
                source: e,
            })?;

        let mut donors = HashSet::new();
        let mut acceptors = HashSet::new();

        for key in params.hbond.keys() {
            let parts: Vec<&str> = key.splitn(2, '-').collect();
            if parts.len() == 2 {
                donors.insert(parts[0].to_string());
                acceptors.insert(parts[1].to_string());
            }
        }

        params.hbond_donors = donors;
        params.hbond_acceptors = acceptors;

        Ok(params)
    }

    /// Loads delta parameters from a CSV file.
    ///
    /// This method reads atomic correction parameters from a CSV file and
    /// indexes them by residue type and atom name for efficient lookup.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the CSV file containing delta parameters.
    ///
    /// # Return
    ///
    /// Returns a map of delta parameters indexed by (residue_type, atom_name).
    ///
    /// # Errors
    ///
    /// Returns a `ParamLoadError` if the file cannot be read or parsed.
    fn load_delta_csv(
        path: &Path,
    ) -> Result<HashMap<(String, String), DeltaParam>, ParamLoadError> {
        let mut reader = csv::Reader::from_path(path).map_err(|e| ParamLoadError::Csv {
            path: path.to_string_lossy().to_string(),
            source: e,
        })?;

        let mut lib_deltas = HashMap::new();
        for result in reader.deserialize::<DeltaParam>() {
            let record = result.map_err(|e| ParamLoadError::Csv {
                path: path.to_string_lossy().to_string(),
                source: e,
            })?;
            lib_deltas.insert(
                (record.residue_type.clone(), record.atom_name.clone()),
                record,
            );
        }
        Ok(lib_deltas)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::atom::AtomRole;
    use std::fs::{self, File};
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn load_non_bonded_succeeds_with_valid_toml() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test.toml");
        let mut file = File::create(&file_path).unwrap();
        writeln!(
            file,
            r#"
            [globals]
            dielectric_constant = 1.0
            potential_function = "buckingham-exp-6"

            [vdw.C]
            radius = 3.5
            well_depth = 0.1
            scale = 12.0

            [vdw.N]
            radius = 3.2
            well_depth = 0.05

            [hbond.N-H]
            equilibrium_distance = 2.7
            well_depth = 5.0
            "#
        )
        .unwrap();

        let params = Forcefield::load_non_bonded(&file_path).unwrap();
        assert_eq!(params.globals.dielectric_constant, 1.0);
        assert_eq!(
            params.vdw.get("C"),
            Some(&VdwParam::Buckingham {
                radius: 3.5,
                well_depth: 0.1,
                scale: 12.0
            })
        );
        assert_eq!(
            params.vdw.get("N"),
            Some(&VdwParam::LennardJones {
                radius: 3.2,
                well_depth: 0.05,
            })
        );
        assert_eq!(
            params.hbond.get("N-H"),
            Some(&HBondParam {
                equilibrium_distance: 2.7,
                well_depth: 5.0
            })
        );
        assert!(params.hbond_donors.contains("N"));
        assert!(params.hbond_acceptors.contains("H"));
        assert_eq!(params.hbond_donors.len(), 1);
        assert_eq!(params.hbond_acceptors.len(), 1);
    }

    #[test]
    fn load_non_bonded_fails_for_missing_file() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("non_existent.toml");
        let result = Forcefield::load_non_bonded(&file_path);
        assert!(matches!(result, Err(ParamLoadError::Io { .. })));
    }

    #[test]
    fn load_non_bonded_fails_for_malformed_toml() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("malformed.toml");
        fs::write(&file_path, "this is not toml").unwrap();
        let result = Forcefield::load_non_bonded(&file_path);
        assert!(matches!(result, Err(ParamLoadError::Toml { .. })));
    }

    #[test]
    fn load_delta_csv_succeeds_with_valid_csv() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test_lib.csv");
        fs::write(
            &file_path,
            "residue_type,atom_name,mu,sigma\nALA,CA,1.0,0.5",
        )
        .unwrap();

        let deltas = Forcefield::load_delta_csv(&file_path).unwrap();
        let param = deltas.get(&("ALA".to_string(), "CA".to_string())).unwrap();
        assert_eq!(param.mu, 1.0);
        assert_eq!(param.sigma, 0.5);
    }

    #[test]
    fn load_delta_csv_fails_for_malformed_csv() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("malformed.csv");
        fs::write(&file_path, "header1,header2\nval1").unwrap();
        let result = Forcefield::load_delta_csv(&file_path);
        assert!(matches!(result, Err(ParamLoadError::Csv { .. })));
    }

    #[test]
    fn load_forcefield_succeeds_with_valid_files() {
        let dir = tempdir().unwrap();

        let non_bonded_path = dir.path().join("non_bonded.toml");
        fs::write(
            &non_bonded_path,
            r#"[globals]
            dielectric_constant = 1.0
            potential_function = "lennard-jones-12-6"
            [vdw.C]
            radius = 1.0
            well_depth = 1.0
            [hbond.N]
            equilibrium_distance = 1.0
            well_depth = 1.0"#,
        )
        .unwrap();

        let delta_path = dir.path().join("delta.csv");
        fs::write(
            &delta_path,
            "residue_type,atom_name,mu,sigma\nALA,CA,1.0,0.5",
        )
        .unwrap();

        let ff =
            Forcefield::load(&non_bonded_path, &delta_path, &EnergyWeights::default()).unwrap();

        assert!(!ff.non_bonded.vdw.is_empty());
        assert!(!ff.deltas.is_empty());
        assert!(ff.weight_map.is_empty());
    }

    #[test]
    fn load_forcefield_fails_if_any_file_is_missing() {
        let dir = tempdir().unwrap();
        let non_bonded_path = dir.path().join("non_bonded.toml");
        let delta_path = dir.path().join("delta.csv");
        fs::write(&non_bonded_path, "").unwrap();

        let result = Forcefield::load(
            &non_bonded_path,
            &delta_path, // delta.csv does not exist
            &EnergyWeights::default(),
        );
        assert!(matches!(result, Err(ParamLoadError::Toml { .. })));
    }

    #[test]
    fn hbond_donor_acceptor_sets_with_multiple_entries() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("hbonds.toml");
        fs::write(
            &file_path,
            r#"
            [globals]
            dielectric_constant = 2.0
            potential_function = "lennard-jones-12-6"

            [vdw.C]
            radius = 1.0
            well_depth = 1.0

            [hbond.N-H]
            equilibrium_distance = 2.0
            well_depth = 3.0

            [hbond.O-HX-Extra]
            equilibrium_distance = 2.5
            well_depth = 4.0
        "#,
        )
        .unwrap();

        let params = Forcefield::load_non_bonded(&file_path).unwrap();

        assert!(params.hbond.contains_key("N-H"));
        assert!(params.hbond.contains_key("O-HX-Extra"));
        assert!(params.hbond_donors.contains("N"));
        assert!(params.hbond_donors.contains("O"));
        assert!(params.hbond_acceptors.contains("H"));
        assert!(params.hbond_acceptors.contains("HX-Extra"));
        assert_eq!(params.hbond_donors.len(), 2);
        assert_eq!(params.hbond_acceptors.len(), 2);
    }

    #[test]
    fn hbond_donor_acceptor_sets_empty_without_hyphen() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("hbonds_no_hyphen.toml");
        fs::write(
            &file_path,
            r#"
            [globals]
            dielectric_constant = 3.0
            potential_function = "lennard-jones-12-6"

            [vdw.C]
            radius = 1.0
            well_depth = 1.0

            [hbond.N_H]
            equilibrium_distance = 2.0
            well_depth = 3.0
        "#,
        )
        .unwrap();

        let params = Forcefield::load_non_bonded(&file_path).unwrap();
        assert!(params.hbond.contains_key("N_H"));
        assert!(params.hbond_donors.is_empty());
        assert!(params.hbond_acceptors.is_empty());
    }

    #[test]
    fn load_forcefield_builds_weight_map_from_rules() {
        let dir = tempdir().unwrap();

        let nb = dir.path().join("nb.toml");
        fs::write(
            &nb,
            r#"[globals]
            dielectric_constant = 1.0
            potential_function = "lennard-jones-12-6"
            [vdw.C]
            radius = 1.0
            well_depth = 1.0
            [hbond.N]
            equilibrium_distance = 1.0
            well_depth = 1.0
        "#,
        )
        .unwrap();
        let delta = dir.path().join("delta.csv");
        fs::write(&delta, "residue_type,atom_name,mu,sigma\nALA,CA,0.0,0.0").unwrap();

        let rules = vec![
            WeightRule {
                groups: [AtomRole::Backbone, AtomRole::Sidechain],
                weights: EnergyComponentWeights {
                    vdw: 0.5,
                    coulomb: 0.2,
                    hbond: 0.1,
                },
            },
            WeightRule {
                groups: [AtomRole::Ligand, AtomRole::Backbone],
                weights: EnergyComponentWeights {
                    vdw: 1.5,
                    coulomb: 0.8,
                    hbond: 0.0,
                },
            },
        ];
        let ff = Forcefield::load(&nb, &delta, &EnergyWeights { rules }).unwrap();

        assert_eq!(ff.weight_map.len(), 2);

        let key_bs = if AtomRole::Backbone <= AtomRole::Sidechain {
            (AtomRole::Backbone, AtomRole::Sidechain)
        } else {
            (AtomRole::Sidechain, AtomRole::Backbone)
        };
        let w_bs = ff.weight_map.get(&key_bs).copied().unwrap();
        assert!((w_bs.vdw - 0.5).abs() < 1e-12);
        assert!((w_bs.coulomb - 0.2).abs() < 1e-12);
        assert!((w_bs.hbond - 0.1).abs() < 1e-12);

        let key_bl = if AtomRole::Backbone <= AtomRole::Ligand {
            (AtomRole::Backbone, AtomRole::Ligand)
        } else {
            (AtomRole::Ligand, AtomRole::Backbone)
        };
        let w_bl = ff.weight_map.get(&key_bl).copied().unwrap();
        assert!((w_bl.vdw - 1.5).abs() < 1e-12);
        assert!((w_bl.coulomb - 0.8).abs() < 1e-12);
        assert!((w_bl.hbond - 0.0).abs() < 1e-12);
    }

    #[test]
    fn load_forcefield_with_empty_rules_has_empty_weight_map() {
        let dir = tempdir().unwrap();
        let nb = dir.path().join("nb2.toml");
        fs::write(
            &nb,
            r#"[globals]
            dielectric_constant = 1.0
            potential_function = "lennard-jones-12-6"
            [vdw.C]
            radius = 1.0
            well_depth = 1.0
            [hbond]
        "#,
        )
        .unwrap();
        let delta = dir.path().join("delta2.csv");
        fs::write(&delta, "residue_type,atom_name,mu,sigma\nALA,CA,0.0,0.0").unwrap();

        let ff = Forcefield::load(&nb, &delta, &EnergyWeights::default()).unwrap();
        assert!(ff.weight_map.is_empty());
    }
}
