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

fn calculate_transformation(
    from_points: &[Point3<f64>],
    to_points: &[Point3<f64>],
) -> Result<(Rotation3<f64>, Vector3<f64>), PlacementError> {
    // 1. Calculate centroids using fold for robust summation
    let from_centroid_sum: Vector3<f64> = from_points
        .iter()
        .fold(Vector3::zeros(), |acc, p| acc + p.coords);
    let from_centroid = Point3::from(from_centroid_sum / from_points.len() as f64);

    let to_centroid_sum: Vector3<f64> = to_points
        .iter()
        .fold(Vector3::zeros(), |acc, p| acc + p.coords);
    let to_centroid = Point3::from(to_centroid_sum / to_points.len() as f64);

    // 2. Center the points
    let centered_from: Vec<_> = from_points.iter().map(|p| p - from_centroid).collect();
    let centered_to: Vec<_> = to_points.iter().map(|p| p - to_centroid).collect();

    // 3. Build the covariance matrix H
    let h = centered_from
        .iter()
        .zip(centered_to.iter())
        .fold(Matrix3::zeros(), |acc, (f, t)| acc + t * f.transpose());

    // 4. Perform SVD
    let svd = h.svd(true, true);
    let u = svd.u.unwrap();
    let v_t = svd.v_t.unwrap();

    // 5. Calculate the rotation matrix, handling reflections
    let mut d = (u * v_t).determinant();
    if d < 0.0 {
        d = -1.0;
    } else {
        d = 1.0;
    }

    let mut correction = Matrix3::identity();
    correction[(2, 2)] = d;

    let rotation_matrix = u * correction * v_t;
    let rotation = Rotation3::from_matrix(&rotation_matrix);

    // 6. Calculate the translation vector
    let translation = to_centroid.coords - rotation * from_centroid.coords;

    Ok((rotation, translation))
}
