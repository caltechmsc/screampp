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
    pub equilibrium_dist: f64,
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
    pub residue_type: String,
    pub atom_name: String,
    pub mu: f64,
    pub sigma: f64,
}

#[derive(Debug, Deserialize, Clone)]
pub struct ChargeParam {
    pub residue_type: String,
    pub atom_name: String,
    pub partial_charge: f64,
}

#[derive(Debug, Deserialize, Clone)]
pub struct TopologyAtomParam {
    pub name: String,
    pub ff_type: String,
}

#[derive(Debug, Deserialize, Clone)]
pub struct TopologyResidueParams {
    pub atoms: Vec<TopologyAtomParam>,
    pub bonds: Vec<[String; 2]>,
}

#[derive(Debug, Clone)]
pub struct Forcefield {
    pub non_bonded: NonBondedParams,
    pub deltas: HashMap<String, HashMap<(String, String), DeltaParam>>,
    pub charges: HashMap<String, HashMap<(String, String), ChargeParam>>,
    pub topology: HashMap<String, TopologyResidueParams>,
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
    #[error("Directory reading error for '{path}': {source}")]
    ReadDir {
        path: String,
        source: std::io::Error,
    },
}

impl Forcefield {
    pub fn load(base_path: &Path) -> Result<Self, ParamLoadError> {
        let non_bonded =
            Self::load_non_bonded(&base_path.join("forcefield/buckingham_exp_6.toml"))?;
        let deltas = Self::load_delta_directory(&base_path.join("delta/"))?;
        let charges = Self::load_charge_directory(&base_path.join("charges/"))?;
        let topology = Self::load_topology(&base_path.join("topology/topology.toml"))?;

        Ok(Self {
            non_bonded,
            deltas,
            charges,
            topology,
        })
    }

    pub fn load_non_bonded(path: &Path) -> Result<NonBondedParams, ParamLoadError> {
        let content = std::fs::read_to_string(path).map_err(|e| ParamLoadError::Io {
            path: path.to_string_lossy().to_string(),
            source: e,
        })?;
        toml::from_str(&content).map_err(|e| ParamLoadError::Toml {
            path: path.to_string_lossy().to_string(),
            source: e,
        })
    }

    pub fn load_delta_directory(
        dir_path: &Path,
    ) -> Result<HashMap<String, HashMap<(String, String), DeltaParam>>, ParamLoadError> {
        let mut all_deltas = HashMap::new();
        for entry in std::fs::read_dir(dir_path).map_err(|e| ParamLoadError::ReadDir {
            path: dir_path.to_string_lossy().to_string(),
            source: e,
        })? {
            let entry = entry.map_err(|e| ParamLoadError::ReadDir {
                path: dir_path.to_string_lossy().to_string(),
                source: e,
            })?;
            let path = entry.path();
            if path.is_file() && path.extension().and_then(|s| s.to_str()) == Some("csv") {
                let library_key = path.file_stem().unwrap().to_string_lossy().to_string();
                let mut reader =
                    csv::Reader::from_path(&path).map_err(|e| ParamLoadError::Csv {
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
                all_deltas.insert(library_key, lib_deltas);
            }
        }
        Ok(all_deltas)
    }

    pub fn load_charge_directory(
        dir_path: &Path,
    ) -> Result<HashMap<String, HashMap<(String, String), ChargeParam>>, ParamLoadError> {
        let mut all_charges = HashMap::new();

        if !dir_path.exists() || !dir_path.is_dir() {
            return Ok(all_charges);
        }

        for entry in std::fs::read_dir(dir_path).map_err(|e| ParamLoadError::ReadDir {
            path: dir_path.to_string_lossy().to_string(),
            source: e,
        })? {
            let entry = entry.map_err(|e| ParamLoadError::ReadDir {
                path: dir_path.to_string_lossy().to_string(),
                source: e,
            })?;
            let path = entry.path();

            if path.is_file() && path.extension().and_then(|s| s.to_str()) == Some("csv") {
                let charge_scheme_key = path.file_stem().unwrap().to_string_lossy().to_string();

                let mut reader =
                    csv::Reader::from_path(&path).map_err(|e| ParamLoadError::Csv {
                        path: path.to_string_lossy().to_string(),
                        source: e,
                    })?;

                let mut scheme_charges = HashMap::new();
                for result in reader.deserialize::<ChargeParam>() {
                    let record = result.map_err(|e| ParamLoadError::Csv {
                        path: path.to_string_lossy().to_string(),
                        source: e,
                    })?;
                    scheme_charges.insert(
                        (record.residue_type.clone(), record.atom_name.clone()),
                        record,
                    );
                }
                all_charges.insert(charge_scheme_key, scheme_charges);
            }
        }
        Ok(all_charges)
    }

    pub fn load_topology(
        path: &Path,
    ) -> Result<HashMap<String, TopologyResidueParams>, ParamLoadError> {
        let content = std::fs::read_to_string(path).map_err(|e| ParamLoadError::Io {
            path: path.to_string_lossy().to_string(),
            source: e,
        })?;
        let topo_params: HashMap<String, TopologyResidueParams> = toml::from_str(&content)
            .map_err(|e| ParamLoadError::Toml {
                path: path.to_string_lossy().to_string(),
                source: e,
            })?;
        Ok(topo_params)
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
            potential_function = "buckingham_exp_6"

            [vdw.C]
            radius = 3.5
            well_depth = 0.1
            scale = 12.0

            [hbond.N]
            equilibrium_dist = 2.7
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
            params.hbond.get("N"),
            Some(&HBondParam {
                equilibrium_dist: 2.7,
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
    fn load_delta_directory_succeeds_with_valid_csv() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test_lib.csv");
        fs::write(
            &file_path,
            "residue_type,atom_name,mu,sigma\nALA,CA,1.0,0.5",
        )
            .unwrap();

        let deltas = Forcefield::load_delta_directory(dir.path()).unwrap();
        assert!(deltas.contains_key("test_lib"));
        let lib = deltas.get("test_lib").unwrap();
        let param = lib.get(&("ALA".to_string(), "CA".to_string())).unwrap();
        assert_eq!(param.mu, 1.0);
        assert_eq!(param.sigma, 0.5);
    }

    #[test]
    fn load_delta_directory_fails_for_malformed_csv() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("malformed.csv");
        fs::write(&file_path, "header1,header2\nval1").unwrap();
        let result = Forcefield::load_delta_directory(dir.path());
        assert!(matches!(result, Err(ParamLoadError::Csv { .. })));
    }

    #[test]
    fn load_charge_directory_succeeds_with_valid_csv() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test_scheme.csv");
        fs::write(
            &file_path,
            "residue_type,atom_name,partial_charge\nGLY,N,-0.5",
        )
            .unwrap();

        let charges = Forcefield::load_charge_directory(dir.path()).unwrap();
        assert!(charges.contains_key("test_scheme"));
        let scheme = charges.get("test_scheme").unwrap();
        let param = scheme.get(&("GLY".to_string(), "N".to_string())).unwrap();
        assert_eq!(param.partial_charge, -0.5);
    }

    #[test]
    fn load_charge_directory_returns_empty_map_for_non_existent_path() {
        let dir = tempdir().unwrap();
        let non_existent_path = dir.path().join("non_existent_dir");
        let charges = Forcefield::load_charge_directory(&non_existent_path).unwrap();
        assert!(charges.is_empty());
    }

    #[test]
    fn load_topology_succeeds_with_valid_toml() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("topology.toml");
        fs::write(
            &file_path,
            r#"
            [ALA]
            atoms = [ { name = "N", ff_type = "N" }, { name = "CA", ff_type = "C" } ]
            bonds = [ ["N", "CA"] ]
            "#,
        )
            .unwrap();

        let topology = Forcefield::load_topology(&file_path).unwrap();
        assert!(topology.contains_key("ALA"));
        let ala_topo = topology.get("ALA").unwrap();
        assert_eq!(ala_topo.atoms.len(), 2);
        assert_eq!(ala_topo.bonds.len(), 1);
        assert_eq!(ala_topo.bonds[0], ["N", "CA"]);
    }

    #[test]
    fn load_forcefield_succeeds_with_valid_data_directory() {
        let base_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("data");
        let ff = Forcefield::load(&base_path).unwrap();

        assert!(!ff.non_bonded.vdw.is_empty());
        assert!(!ff.deltas.is_empty());
        assert!(!ff.charges.is_empty());
        assert!(!ff.topology.is_empty());
        assert!(ff.topology.contains_key("ALA"));
    }
}
