use nalgebra::{Matrix3, Point3, Rotation3, Vector3};
use thiserror::Error;

use super::error::EngineError;
use crate::core::{
    models::{
        atom::Atom,
        ids::{AtomId, ResidueId},
        system::MolecularSystem,
    },
    rotamers::{placement::PlacementInfo, rotamer::Rotamer},
};

#[derive(Debug, Error)]
pub enum PlacementError {
    #[error("Anchor atom '{atom_name}' not found in the target residue {residue_id:?}")]
    AnchorAtomNotFoundInSystem {
        atom_name: String,
        residue_id: ResidueId,
    },
    #[error("Anchor atom '{atom_name}' not found in the rotamer template")]
    AnchorAtomNotFoundInRotamer { atom_name: String },

    #[error("Insufficient anchor atoms: requires at least 3, but found {found}")]
    InsufficientAnchors { found: usize },

    #[error(
        "Side-chain atom '{atom_name}' is defined in PlacementInfo but not found in the rotamer template"
    )]
    SideChainAtomNotFoundInRotamer { atom_name: String },
}
