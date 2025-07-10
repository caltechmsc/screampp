use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct PlacementInfo {
    pub anchor_atoms: Vec<String>,
    pub sidechain_atoms: Vec<String>,
    pub exact_match_atoms: Vec<String>,
    pub connection_points: Vec<String>,
}

pub type PlacementRegistry = HashMap<String, PlacementInfo>;

#[derive(Debug, Error)]
pub enum PlacementLoadError {
    #[error("File I/O error for '{path}': {source}")]
    Io {
        path: String,
        source: std::io::Error,
    },
    #[error("TOML parsing error for '{path}': {source}")]
    Toml {
        path: String,
        source: toml::de::Error,
    },
}

pub fn load_placement_registry(path: &Path) -> Result<PlacementRegistry, PlacementLoadError> {
    let content = std::fs::read_to_string(path).map_err(|e| PlacementLoadError::Io {
        path: path.to_string_lossy().to_string(),
        source: e,
    })?;
    toml::from_str(&content).map_err(|e| PlacementLoadError::Toml {
        path: path.to_string_lossy().to_string(),
        source: e,
    })
}
