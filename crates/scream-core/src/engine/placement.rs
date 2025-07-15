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

impl From<PlacementError> for EngineError {
    fn from(e: PlacementError) -> Self {
        EngineError::Placement {
            residue_id: ResidueId::default(), // Placeholder
            message: e.to_string(),
        }
    }
}

pub fn place_rotamer_on_system(
    system: &mut MolecularSystem,
    target_residue_id: ResidueId,
    rotamer: &Rotamer,
    placement_info: &PlacementInfo,
) -> Result<(), PlacementError> {
    // --- 1. Gather and validate anchor points ---
    let (system_anchors, rotamer_anchors) =
        gather_anchor_points(system, target_residue_id, rotamer, placement_info)?;

    // --- 2. Calculate the rigid body transformation ---
    let (rotation, translation) = calculate_transformation(&rotamer_anchors, &system_anchors)?;

    // --- 3. Prepare new side-chain atoms ---
    let new_atoms_to_add = prepare_new_sidechain_atoms(
        rotamer,
        placement_info,
        target_residue_id,
        rotation,
        translation,
    )?;

    // --- 4. Remove old side-chain atoms ---
    let old_atom_ids_to_remove: Vec<AtomId> = {
        let target_residue = system.residue(target_residue_id).unwrap();
        placement_info
            .sidechain_atoms
            .iter()
            .filter_map(|atom_name| target_residue.get_atom_id_by_name(atom_name))
            .collect()
    };

    for atom_id in old_atom_ids_to_remove {
        system.remove_atom(atom_id);
    }

    // --- 5. Add new side-chain atoms ---
    for atom in new_atoms_to_add {
        system.add_atom_to_residue(target_residue_id, atom);
    }

    Ok(())
}

fn gather_anchor_points<'a>(
    system: &'a MolecularSystem,
    target_residue_id: ResidueId,
    rotamer: &'a Rotamer,
    placement_info: &PlacementInfo,
) -> Result<(Vec<Point3<f64>>, Vec<Point3<f64>>), PlacementError> {
    let mut system_points = Vec::with_capacity(placement_info.anchor_atoms.len());
    let mut rotamer_points = Vec::with_capacity(placement_info.anchor_atoms.len());
    let target_residue = system.residue(target_residue_id).unwrap();

    for atom_name in &placement_info.anchor_atoms {
        let system_atom_id = target_residue
            .get_atom_id_by_name(atom_name)
            .ok_or_else(|| PlacementError::AnchorAtomNotFoundInSystem {
                atom_name: atom_name.clone(),
                residue_id: target_residue_id,
            })?;
        let system_atom = system.atom(system_atom_id).unwrap();
        system_points.push(system_atom.position);

        let rotamer_atom = rotamer
            .atoms
            .iter()
            .find(|a| &a.name == atom_name)
            .ok_or_else(|| PlacementError::AnchorAtomNotFoundInRotamer {
                atom_name: atom_name.clone(),
            })?;
        rotamer_points.push(rotamer_atom.position);
    }

    if system_points.len() < 3 {
        return Err(PlacementError::InsufficientAnchors {
            found: system_points.len(),
        });
    }

    Ok((system_points, rotamer_points))
}
