use crate::core::forcefield::params::Forcefield;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::core::models::topology::BondOrder;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ParameterizationError {
    #[error("Missing topology for residue type: {0}")]
    MissingTopology(String),
    #[error("Missing force field type for atom '{atom_name}' in residue '{res_type}'")]
    MissingForceFieldType { res_type: String, atom_name: String },
    #[error("Missing charge parameter for atom '{atom_name}' in residue '{res_type}'")]
    MissingCharge { res_type: String, atom_name: String },
    #[error("Missing VDW parameter for force field type: {0}")]
    MissingVdwParams(String),
    #[error("Missing delta parameter for atom '{atom_name}' in residue '{res_type}'")]
    MissingDelta { res_type: String, atom_name: String },
    #[error(
        "Inconsistent topology: Atom '{atom_name}' from input file not found in topology for residue '{res_type}'"
    )]
    AtomNotFoundInTopology { res_type: String, atom_name: String },
    #[error("Atom with name '{0}' not found in residue ID {1:?}")]
    AtomNameNotFound(String, ResidueId),
}
