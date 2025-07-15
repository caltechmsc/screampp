use thiserror::Error;

use super::config::ResidueSpecifier;
use crate::core::forcefield::scoring::ScoringError;
use crate::core::models::ids::ResidueId;

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
}
