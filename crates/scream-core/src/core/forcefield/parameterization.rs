use super::params::Forcefield;
use crate::core::models::atom::Atom;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::core::models::system::MolecularSystem;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ParameterizationError {
    #[error("Missing topology for residue type '{0}'")]
    MissingTopology(String),
    #[error("Missing force field type for atom '{atom_name}' in residue '{res_type}'")]
    MissingForcefieldType { res_type: String, atom_name: String },
    #[error("Missing charge parameters for scheme '{0}'")]
    MissingChargeScheme(String),
    #[error("Missing charge for atom '{atom_name}' in residue '{res_type}' for scheme '{scheme}'")]
    MissingCharge {
        scheme: String,
        res_type: String,
        atom_name: String,
    },
    #[error("Missing delta parameters for library '{0}'")]
    MissingDeltaLibrary(String),
    #[error(
        "Missing delta for atom '{atom_name}' in residue '{res_type}' for library '{library_key}'"
    )]
    MissingDelta {
        library_key: String,
        res_type: String,
        atom_name: String,
    },
}
