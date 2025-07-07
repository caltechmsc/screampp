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

#[derive(Debug, Deserialize, Clone)]
pub struct DeltaParam {
    pub residue_type: String,
    pub atom_name: String,
    pub mu: f64,
    pub sigma: f64,
}

#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct Globals {
    pub dielectric_constant: f64,
    pub potential_function: String,
}

#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct ForcefieldParams {
    pub globals: Globals,
    pub vdw: HashMap<String, VdwParam>,
    pub hbond: HashMap<String, HBondParam>,
}

#[derive(Debug, Error)]
pub enum ParamLoadError {
    #[error("File I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("TOML deserialization error: {0}")]
    Toml(#[from] toml::de::Error),
}

impl ForcefieldParams {
    pub fn load_from_toml(path: &Path) -> Result<Self, ParamLoadError> {
        let content = std::fs::read_to_string(path)?;
        let params: ForcefieldParams = toml::from_str(&content)?;
        Ok(params)
    }
}
