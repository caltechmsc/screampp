use super::config::ResidueSpecifier;
use crate::core::forcefield::params::ParamLoadError;
use crate::core::forcefield::scoring::ScoringError;
use crate::core::models::ids::ResidueId;
use crate::core::rotamers::library::LibraryLoadError;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum EngineError {
    #[error("Initialization failed: {0}")]
    Initialization(String),

    #[error("Residue not found in system: {spec:?}")]
    ResidueNotFound { spec: ResidueSpecifier },

    #[error("Energy scoring failed: {source}")]
    Scoring {
        #[from]
        source: ScoringError,
    },

    #[error("Topology not found for residue: {residue_name}")]
    TopologyNotFound{
        residue_name: String,
    },

    #[error("Rotamer library error for residue {residue_type:?}: {message}")]
    RotamerLibrary {
        residue_type: String,
        message: String,
    },

    #[error("Failed to place rotamer on residue {residue_id:?}: {message}")]
    Placement {
        residue_id: ResidueId,
        message: String,
    },

    #[error("Optimization phase '{phase}' failed: {reason}")]
    PhaseFailed { phase: &'static str, reason: String },

    #[error("Algorithm failed to converge after {iterations} iterations")]
    Convergence { iterations: usize },

    #[error("Internal logic error: {0}")]
    Internal(String),

    #[error("Failed to load forcefield parameters: {0}")]
    ParamLoad(#[from] ParamLoadError),

    #[error("Failed to load rotamer library: {0}")]
    LibraryLoad(#[from] LibraryLoadError),
}
