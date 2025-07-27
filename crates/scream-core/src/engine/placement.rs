use super::error::EngineError;
use crate::core::{
    models::{
        atom::Atom,
        ids::{AtomId, ResidueId},
        system::MolecularSystem,
    },
    rotamers::{placement::PlacementInfo, rotamer::Rotamer},
};
use nalgebra::{Matrix3, Point3, Rotation3, Vector3};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum PlacementError {
    #[error("Anchor atom '{atom_name}' not found in the target residue")]
    AnchorAtomNotFoundInSystem { atom_name: String },

    #[error("Anchor atom '{atom_name}' not found in the rotamer template")]
    AnchorAtomNotFoundInRotamer { atom_name: String },

    #[error(
        "Insufficient anchor atoms for stable alignment: requires at least 3, but found {found}"
    )]
    InsufficientAnchors { found: usize },

    #[error(
        "Side-chain atom '{atom_name}' is defined in PlacementInfo but not found in the rotamer template"
    )]
    SideChainAtomNotFoundInRotamer { atom_name: String },
}

pub fn place_rotamer_on_system(
    system: &mut MolecularSystem,
    target_residue_id: ResidueId,
    rotamer: &Rotamer,
    placement_info: &PlacementInfo,
) -> Result<(), EngineError> {
    // --- 1. Gather and validate anchor points ---
    let (system_anchors, rotamer_anchors) =
        gather_anchor_points(system, target_residue_id, rotamer, placement_info).map_err(|e| {
            EngineError::Placement {
                residue_id: target_residue_id,
                message: e.to_string(),
            }
        })?;

    // --- 2. Calculate the rigid body transformation ---
    let (rotation, translation) = calculate_transformation(&rotamer_anchors, &system_anchors)
        .map_err(|e| EngineError::Placement {
            residue_id: target_residue_id,
            message: e.to_string(),
        })?;

    // --- 3. Prepare new side-chain atoms ---
    let new_atoms_to_add = prepare_new_sidechain_atoms(
        rotamer,
        placement_info,
        target_residue_id,
        rotation,
        translation,
    )
    .map_err(|e| EngineError::Placement {
        residue_id: target_residue_id,
        message: e.to_string(),
    })?;

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

fn gather_anchor_points(
    system: &MolecularSystem,
    target_residue_id: ResidueId,
    rotamer: &Rotamer,
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
    let d = (u * v_t.transpose()).determinant();
    let mut correction = Matrix3::identity();
    if d < 0.0 {
        correction[(2, 2)] = -1.0;
    }

    let rotation_matrix = u * correction * v_t;
    let rotation = Rotation3::from_matrix(&rotation_matrix);

    // 6. Calculate the translation vector
    let translation = to_centroid.coords - rotation * from_centroid.coords;

    Ok((rotation, translation))
}

fn prepare_new_sidechain_atoms(
    rotamer: &Rotamer,
    placement_info: &PlacementInfo,
    target_residue_id: ResidueId,
    rotation: Rotation3<f64>,
    translation: Vector3<f64>,
) -> Result<Vec<Atom>, PlacementError> {
    let mut new_atoms = Vec::with_capacity(placement_info.sidechain_atoms.len());

    for atom_name in &placement_info.sidechain_atoms {
        let rotamer_atom = rotamer
            .atoms
            .iter()
            .find(|a| &a.name == atom_name)
            .ok_or_else(|| PlacementError::SideChainAtomNotFoundInRotamer {
                atom_name: atom_name.clone(),
            })?;

        let mut new_atom = rotamer_atom.clone();
        new_atom.residue_id = target_residue_id;
        new_atom.position = rotation * rotamer_atom.position + translation;
        new_atoms.push(new_atom);
    }
    Ok(new_atoms)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::chain::ChainType;
    use nalgebra::{Point3, Vector3};

    fn create_test_system_with_residue() -> (MolecularSystem, ResidueId) {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let residue_id = system.add_residue(chain_id, 1, "ALA", None).unwrap();

        let n = Atom::new("N", residue_id, Point3::new(0.0, 1.0, 0.0));
        let ca = Atom::new("CA", residue_id, Point3::new(0.0, 0.0, 0.0));
        let c = Atom::new("C", residue_id, Point3::new(1.0, 0.0, 0.0));
        system.add_atom_to_residue(residue_id, n);
        system.add_atom_to_residue(residue_id, ca);
        system.add_atom_to_residue(residue_id, c);

        let old_cb = Atom::new("CB", residue_id, Point3::new(-1.0, -1.0, -1.0));
        system.add_atom_to_residue(residue_id, old_cb);

        (system, residue_id)
    }

    fn create_test_rotamer() -> Rotamer {
        Rotamer {
            atoms: vec![
                Atom::new("N", ResidueId::default(), Point3::new(5.0, 5.0, 6.0)),
                Atom::new("CA", ResidueId::default(), Point3::new(5.0, 5.0, 5.0)),
                Atom::new("C", ResidueId::default(), Point3::new(6.0, 5.0, 5.0)),
                Atom::new("CB", ResidueId::default(), Point3::new(4.0, 4.0, 4.0)),
            ],
        }
    }

    fn create_test_placement_info() -> PlacementInfo {
        PlacementInfo {
            anchor_atoms: vec!["N".to_string(), "CA".to_string(), "C".to_string()],
            sidechain_atoms: vec!["CB".to_string()],
            exact_match_atoms: vec![],
            connection_points: vec![],
        }
    }

    #[test]
    fn test_calculate_transformation_identity() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let (rot, trans) = calculate_transformation(&points, &points).unwrap();

        assert!(
            rot.angle().abs() < 1e-9,
            "Rotation angle should be zero for identity"
        );
        assert!(trans.norm() < 1e-9, "Translation should be zero");
    }

    #[test]
    fn test_calculate_transformation_pure_translation() {
        let from_points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let translation_vec = Vector3::new(10.0, 20.0, 30.0);
        let to_points: Vec<Point3<f64>> = from_points.iter().map(|p| p + translation_vec).collect();

        let (rot, trans) = calculate_transformation(&from_points, &to_points).unwrap();

        assert!(
            rot.angle().abs() < 1e-9,
            "Rotation angle should be zero for pure translation"
        );
        assert!(
            (trans - translation_vec).norm() < 1e-9,
            "Translation should match the applied vector"
        );
    }

    #[test]
    fn test_calculate_transformation_pure_rotation() {
        let from_points = vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
        ];

        let rotation = Rotation3::from_axis_angle(&Vector3::z_axis(), std::f64::consts::FRAC_PI_2); // 90 deg around Z
        let to_points: Vec<Point3<f64>> = from_points.iter().map(|p| rotation * p).collect();

        let (rot, trans) = calculate_transformation(&from_points, &to_points).unwrap();

        assert!(
            (rot.matrix() - rotation.matrix()).norm() < 1e-9,
            "Rotation should match the applied rotation"
        );
        assert!(trans.norm() < 1e-9, "Translation should be zero");
    }

    #[test]
    fn test_calculate_transformation_with_reflection_handling() {
        let from_points = vec![
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(-1.0, 1.0, -1.0),
            Point3::new(1.0, -1.0, -1.0),
        ];
        let to_points = vec![
            Point3::new(1.0, 1.0, -1.0),
            Point3::new(-1.0, 1.0, 1.0),
            Point3::new(1.0, -1.0, 1.0),
        ];

        let (rot, _) = calculate_transformation(&from_points, &to_points).unwrap();

        assert!(
            (rot.matrix().determinant() - 1.0).abs() < 1e-9,
            "Determinant should be +1, indicating a proper rotation"
        );
    }

    #[test]
    fn test_place_rotamer_on_system_full_flow() {
        let (mut system, residue_id) = create_test_system_with_residue();
        let rotamer = create_test_rotamer();
        let placement_info = create_test_placement_info();

        let old_cb_id = system
            .residue(residue_id)
            .unwrap()
            .get_atom_id_by_name("CB")
            .unwrap();
        let old_cb_pos = system.atom(old_cb_id).unwrap().position;
        assert!((old_cb_pos - Point3::new(-1.0, -1.0, -1.0)).norm() < 1e-9);
        assert_eq!(system.residue(residue_id).unwrap().atoms().len(), 4);

        let result = place_rotamer_on_system(&mut system, residue_id, &rotamer, &placement_info);
        assert!(result.is_ok());

        let residue = system.residue(residue_id).unwrap();

        assert_eq!(
            residue.atoms().len(),
            4,
            "Atom count should remain 4 (N, CA, C, new CB)"
        );

        let new_cb_id = residue
            .get_atom_id_by_name("CB")
            .expect("New CB atom should exist");
        let new_cb_atom = system.atom(new_cb_id).unwrap();

        let (system_anchors, rotamer_anchors) =
            gather_anchor_points(&system, residue_id, &rotamer, &placement_info).unwrap();
        let (rotation, translation) =
            calculate_transformation(&rotamer_anchors, &system_anchors).unwrap();
        let rotamer_cb_pos = rotamer
            .atoms
            .iter()
            .find(|a| a.name == "CB")
            .unwrap()
            .position;
        let expected_pos = rotation * rotamer_cb_pos + translation;

        let distance = (new_cb_atom.position - expected_pos).norm();
        assert!(
            distance < 1e-9,
            "New CB position is incorrect. Expected: {:?}, Found: {:?}",
            expected_pos,
            new_cb_atom.position
        );

        let rotamer_cb = rotamer.atoms.iter().find(|a| a.name == "CB").unwrap();
        assert_eq!(new_cb_atom.partial_charge, rotamer_cb.partial_charge);
        assert_eq!(new_cb_atom.force_field_type, rotamer_cb.force_field_type);
    }

    #[test]
    fn test_placement_fails_if_system_anchor_missing() {
        let (mut system, residue_id) = create_test_system_with_residue();
        let rotamer = create_test_rotamer();

        let mut bad_placement_info = create_test_placement_info();
        bad_placement_info.anchor_atoms.push("OXT".to_string());

        let result =
            place_rotamer_on_system(&mut system, residue_id, &rotamer, &bad_placement_info);

        assert!(matches!(result, Err(EngineError::Placement { .. })));
    }

    #[test]
    fn test_placement_fails_if_rotamer_anchor_missing() {
        let (mut system, residue_id) = create_test_system_with_residue();
        let mut bad_rotamer = create_test_rotamer();
        bad_rotamer.atoms.retain(|a| a.name != "C");
        let placement_info = create_test_placement_info();

        let result =
            place_rotamer_on_system(&mut system, residue_id, &bad_rotamer, &placement_info);

        assert!(matches!(result, Err(EngineError::Placement { .. })));
    }
}
