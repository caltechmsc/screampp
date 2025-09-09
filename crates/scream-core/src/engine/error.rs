use super::config::ResidueSpecifier;
use crate::core::forcefield::parameterization::ParameterizationError;
use crate::core::forcefield::params::ParamLoadError;
use crate::core::forcefield::scoring::ScoringError;
use crate::core::models::ids::ResidueId;
use crate::core::rotamers::library::LibraryLoadError;
use crate::core::topology::registry::TopologyLoadError;
use thiserror::Error;

/// Represents errors that can occur during SCREAM++ engine operations.
///
/// This enum encompasses all possible failure modes that may arise during
/// protein side-chain placement, optimization, and related molecular modeling
/// operations. Each variant provides specific context about the type of error
/// and the operation that failed.
#[derive(Debug, Error)]
pub enum EngineError {
    /// Indicates that the engine initialization process failed.
    ///
    /// This error occurs when the SCREAM++ engine cannot be properly set up
    /// due to configuration issues, missing dependencies, or invalid parameters.
    /// Common causes include malformed configuration files or incompatible settings.
    #[error("Initialization failed: {0}")]
    Initialization(String),

    /// Indicates that a specified residue could not be found in the molecular system.
    ///
    /// This error occurs when attempting to operate on a residue that doesn't exist
    /// in the current system, possibly due to incorrect residue identifiers or
    /// system state inconsistencies.
    #[error("Residue not found in system: {spec:?}")]
    ResidueNotFound { spec: ResidueSpecifier },

    /// Indicates that an energy scoring operation failed.
    ///
    /// This error wraps failures from the energy calculation subsystem, which
    /// may occur due to invalid atom coordinates, missing force field parameters,
    /// or numerical instabilities during energy evaluation.
    #[error("Energy scoring failed: {source}")]
    Scoring {
        /// The underlying scoring error that occurred.
        #[from]
        source: ScoringError,
    },

    /// Indicates that topology information is missing for a residue type.
    ///
    /// This error occurs when attempting to process a residue for which no
    /// topology definition exists in the registry, preventing proper atom
    /// classification and connectivity analysis.
    #[error("Topology not found for residue: {residue_name}")]
    TopologyNotFound { residue_name: String },

    /// Indicates an error related to rotamer library operations for a specific residue type.
    ///
    /// This error occurs when issues arise with rotamer data access or processing,
    /// such as missing rotamer definitions, corrupted library files, or
    /// incompatible rotamer formats.
    #[error("Rotamer library error for residue {residue_type:?}: {message}")]
    RotamerLibrary {
        /// The residue type that encountered the library error.
        residue_type: String,
        /// Detailed description of the library error.
        message: String,
    },

    /// Indicates that rotamer placement failed on a specific residue.
    ///
    /// This error occurs when the side-chain placement algorithm cannot successfully
    /// position a rotamer on the target residue, typically due to steric clashes,
    /// geometric constraints, or optimization failures.
    #[error("Failed to place rotamer on residue {residue_id:?}: {message}")]
    Placement {
        /// The ID of the residue where placement failed.
        residue_id: ResidueId,
        /// Detailed description of the placement failure.
        message: String,
    },

    /// Indicates that a specific optimization phase failed during execution.
    ///
    /// This error occurs when one of the algorithmic phases in the side-chain
    /// placement pipeline encounters an unrecoverable error, such as numerical
    /// instabilities or constraint violations.
    #[error("Optimization phase '{phase}' failed: {reason}")]
    PhaseFailed { phase: &'static str, reason: String },

    /// Indicates that the optimization algorithm failed to converge within the allowed iterations.
    ///
    /// This error occurs when iterative optimization methods cannot reach the
    /// convergence criteria within the specified maximum number of iterations,
    /// suggesting the problem may be ill-conditioned or require different parameters.
    #[error("Algorithm failed to converge after {iterations} iterations")]
    Convergence { iterations: usize },

    /// Indicates an internal logic error in the engine.
    ///
    /// This error represents unexpected conditions or programming errors that
    /// should not occur under normal operation, such as inconsistent internal state
    /// or violated invariants.
    #[error("Internal logic error: {0}")]
    Internal(String),

    /// Indicates that force field parameter loading failed.
    ///
    /// This error wraps failures from the parameter loading subsystem, which
    /// may occur due to missing parameter files, malformed data, or incompatible
    /// parameter formats.
    #[error("Failed to load forcefield parameters: {0}")]
    ParamLoad(#[from] ParamLoadError),

    /// Indicates that rotamer library loading failed.
    ///
    /// This error wraps failures from the rotamer library loading process,
    /// including file I/O errors, parsing failures, and data validation issues.
    #[error("Failed to load rotamer library: {0}")]
    LibraryLoad(#[from] LibraryLoadError),

    /// Indicates that topology registry loading failed.
    ///
    /// This error wraps failures from the topology loading subsystem, which
    /// may occur due to missing topology files or malformed topology definitions.
    #[error("Failed to load topology registry: {0}")]
    TopologyLoad(#[from] TopologyLoadError),

    /// Indicates that system parameterization failed.
    ///
    /// This error wraps failures from the parameterization process, which assigns
    /// force field parameters to atoms and may fail due to missing parameters
    /// or incompatible atom types.
    #[error("Failed to parameterize system: {0}")]
    Parameterization(#[from] ParameterizationError),
}
