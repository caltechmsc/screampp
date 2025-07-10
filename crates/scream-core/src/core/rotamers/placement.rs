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
