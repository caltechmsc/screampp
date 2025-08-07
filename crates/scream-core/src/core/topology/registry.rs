use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Deserialize, Clone, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ResidueTopology {
    pub anchor_atoms: Vec<String>,
    pub sidechain_atoms: Vec<String>,
}

#[derive(Debug, Clone, Default)]
pub struct TopologyRegistry {
    registry: HashMap<String, ResidueTopology>,
}

impl TopologyRegistry {
    pub fn load(path: &Path) -> Result<Self, TopologyLoadError> {
        let content = std::fs::read_to_string(path).map_err(|e| TopologyLoadError::Io {
            path: path.to_string_lossy().to_string(),
            source: e,
        })?;
        let registry: HashMap<String, ResidueTopology> =
            toml::from_str(&content).map_err(|e| TopologyLoadError::Toml {
                path: path.to_string_lossy().to_string(),
                source: e,
            })?;
        Ok(Self { registry })
    }

    pub fn get(&self, residue_name: &str) -> Option<&ResidueTopology> {
        self.registry.get(residue_name)
    }
}

#[derive(Debug, Error)]
pub enum TopologyLoadError {
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
