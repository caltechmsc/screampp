use crate::core::{
    forcefield::{
        library::LibraryLoadError, parameterization::ParameterizationError,
        placement::PlacementLoadError,
    },
    io::bgf::BgfError,
};
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ScreamError {
    #[error("Configuration Error: {0}")]
    Config(String),

    #[error("File I/O error for path '{path}': {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("Failed to load BGF file: {0}")]
    BfgLoad(#[from] BgfError),

    #[error("Failed to load forcefield parameters: {0}")]
    ForcefieldLoad(String),

    #[error("Failed to load placement registry: {0}")]
    PlacementLoad(#[from] PlacementLoadError),

    #[error("Failed to load rotamer library: {0}")]
    RotamerLibraryLoad(#[from] LibraryLoadError),

    #[error("System or Rotamer Parameterization Error: {0}")]
    Parameterization(#[from] ParameterizationError),

    #[error("Placement Algorithm Error: {0}")]
    PlacementAlgorithm(String),
}
