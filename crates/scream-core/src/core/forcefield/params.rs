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
        unimplemented!("Implement charge directory loading similarly to delta loading.")
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
