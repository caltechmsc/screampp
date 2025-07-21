use serde::Deserialize;
use std::collections::HashMap;
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

#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct GlobalParams {
    pub dielectric_constant: f64,
    pub potential_function: String,
}

#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct NonBondedParams {
    pub globals: GlobalParams,
    pub vdw: HashMap<String, VdwParam>,
    pub hbond: HashMap<String, HBondParam>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct DeltaParam {
    pub res_type: String,
    pub atom_name: String,
    pub mu: f64,
    pub sigma: f64,
}

#[derive(Debug, Clone)]
pub struct Forcefield {
    pub non_bonded: NonBondedParams,
    pub deltas: HashMap<(String, String), DeltaParam>,
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
    pub fn load(non_bonded_path: &Path, delta_path: &Path) -> Result<Self, ParamLoadError> {
        let non_bonded = Self::load_non_bonded(non_bonded_path)?;
        let deltas = Self::load_delta_csv(delta_path)?;

        Ok(Self { non_bonded, deltas })
    }

    fn load_non_bonded(path: &Path) -> Result<NonBondedParams, ParamLoadError> {
        let content = std::fs::read_to_string(path).map_err(|e| ParamLoadError::Io {
            path: path.to_string_lossy().to_string(),
            source: e,
        })?;
        toml::from_str(&content).map_err(|e| ParamLoadError::Toml {
            path: path.to_string_lossy().to_string(),
            source: e,
        })
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
            lib_deltas.insert((record.res_type.clone(), record.atom_name.clone()), record);
        }
        Ok(lib_deltas)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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

            [hbond.N_H]
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
            params.hbond.get("N_H"),
            Some(&HBondParam {
                equilibrium_distance: 2.7,
                well_depth: 5.0
            })
        );
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
        fs::write(&file_path, "res_type,atom_name,mu,sigma\nALA,CA,1.0,0.5").unwrap();

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
        fs::write(&delta_path, "res_type,atom_name,mu,sigma\nALA,CA,1.0,0.5").unwrap();

        let ff = Forcefield::load(&non_bonded_path, &delta_path).unwrap();

        assert!(!ff.non_bonded.vdw.is_empty());
        assert!(!ff.deltas.is_empty());
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
        );
        assert!(matches!(result, Err(ParamLoadError::Toml { .. })));
    }
}
