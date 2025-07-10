use crate::core::{
    forcefield::parameterization::ParameterizationError, forcefield::scoring::ScoringError,
    io::bgf::BgfError, rotamers::library::LibraryLoadError,
};
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ScreamError {
    #[error("Configuration Error: {0}")]
    Config(String),

    #[error("File I/O Error for path '{path}': {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("BGF Parsing Error: {0}")]
    BfgParse(#[from] BgfError),

    #[error("TOML Parsing Error for path '{path}': {source}")]
    TomlParse {
        path: PathBuf,
        #[source]
        source: toml::de::Error,
    },

    #[error("Rotamer Library Error: {0}")]
    RotamerLibrary(#[from] LibraryLoadError),

    #[error("System Parameterization Error: {0}")]
    Parameterization(#[from] ParameterizationError),

    #[error("Scoring Error: {0}")]
    Scoring(#[from] ScoringError),

    #[error("Placement Algorithm Error: {0}")]
    Placement(String),
}

impl From<(PathBuf, std::io::Error)> for ScreamError {
    fn from(value: (PathBuf, std::io::Error)) -> Self {
        ScreamError::Io {
            path: value.0,
            source: value.1,
        }
    }
}
