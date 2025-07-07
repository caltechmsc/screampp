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
    pub res_type: String,
    pub atom_name: String,
    pub mu: f64,
    pub sigma: f64,
}

#[derive(Debug, Deserialize, Clone)]
pub struct ChargeParam {
    pub res_type: String,
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
    pub deltas: HashMap<(String, String), DeltaParam>,
    pub charges: HashMap<(String, String), ChargeParam>,
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
}

impl Forcefield {
    pub fn load(
        non_bonded_path: &Path,
        delta_path: &Path,
        charge_path: &Path,
        topology_path: &Path,
    ) -> Result<Self, ParamLoadError> {
        let non_bonded = Self::load_non_bonded(non_bonded_path)?;
        let deltas = Self::load_delta_csv(delta_path)?;
        let charges = Self::load_charge_csv(charge_path)?;
        let topology = Self::load_topology(topology_path)?;

        Ok(Self {
            non_bonded,
            deltas,
            charges,
            topology,
        })
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

    fn load_charge_csv(
        path: &Path,
    ) -> Result<HashMap<(String, String), ChargeParam>, ParamLoadError> {
        let mut reader = csv::Reader::from_path(path).map_err(|e| ParamLoadError::Csv {
            path: path.to_string_lossy().to_string(),
            source: e,
        })?;

        let mut scheme_charges = HashMap::new();
        for result in reader.deserialize::<ChargeParam>() {
            let record = result.map_err(|e| ParamLoadError::Csv {
                path: path.to_string_lossy().to_string(),
                source: e,
            })?;
            scheme_charges.insert((record.res_type.clone(), record.atom_name.clone()), record);
        }
        Ok(scheme_charges)
    }

    fn load_topology(
        path: &Path,
    ) -> Result<HashMap<String, TopologyResidueParams>, ParamLoadError> {
        let content = std::fs::read_to_string(path).map_err(|e| ParamLoadError::Io {
            path: path.to_string_lossy().to_string(),
            source: e,
        })?;
        toml::from_str(&content).map_err(|e| ParamLoadError::Toml {
            path: path.to_string_lossy().to_string(),
            source: e,
        })
    }
}
