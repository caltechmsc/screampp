use crate::core::models::atom::AtomRole;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Deserialize, Clone, PartialEq)]
#[serde(untagged)]
pub enum VdwParam {
    Buckingham {
        radius: f64,
        well_depth: f64,
        scale: f64,
    },
    LennardJones {
        radius: f64,
        well_depth: f64,
    },
}

#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct HBondParam {
    pub equilibrium_distance: f64,
    pub well_depth: f64,
}

#[derive(Debug, Deserialize, Clone, PartialEq, Default)]
pub struct GlobalParams {
    pub dielectric_constant: f64,
    pub potential_function: String,
}

#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct NonBondedParams {
    pub globals: GlobalParams,
    pub vdw: HashMap<String, VdwParam>,
    pub hbond: HashMap<String, HBondParam>,

    #[serde(skip)]
    pub hbond_donors: HashSet<String>,
    #[serde(skip)]
    pub hbond_acceptors: HashSet<String>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct DeltaParam {
    pub residue_type: String,
    pub atom_name: String,
    pub mu: f64,
    pub sigma: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyComponentWeights {
    pub vdw: f64,
    pub coulomb: f64,
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

#[derive(Debug, Clone, PartialEq)]
pub struct WeightRule {
    pub groups: [AtomRole; 2],
    pub weights: EnergyComponentWeights,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct EnergyWeights {
    pub rules: Vec<WeightRule>,
}

#[derive(Debug, Clone)]
pub struct Forcefield {
    pub non_bonded: NonBondedParams,
    pub deltas: HashMap<(String, String), DeltaParam>,
    pub weight_map: HashMap<(AtomRole, AtomRole), EnergyComponentWeights>,
}

#[derive(Debug, Error)]
pub enum ParamLoadError {
    #[error("File I/O error for '{path}': {source}")]
    Io {
        path: String,
        source: std::io::Error,
    },
    #[error("CSV parsing error for '{path}': {source}")]
    Csv { path: String, source: csv::Error },
    #[error("TOML parsing error for '{path}': {source}")]
    Toml {
        path: String,
        source: toml::de::Error,
    },
}

impl Forcefield {
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
            let key = if a <= b { (a, b) } else { (b, a) };
            weight_map.insert(key, rule.weights);
        }

        Ok(Self {
            non_bonded,
            deltas,
            weight_map,
        })
    }

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
