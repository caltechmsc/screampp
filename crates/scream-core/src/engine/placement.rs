use super::error::EngineError;
use crate::core::{
    models::{
        ids::{AtomId, ResidueId},
        system::MolecularSystem,
        topology::BondOrder,
    },
    rotamers::{placement::PlacementInfo, rotamer::Rotamer},
};
use nalgebra::{Matrix3, Point3, Rotation3, Vector3};
use std::collections::HashMap;
use thiserror::Error;
use tracing::warn;

#[derive(Debug, Error)]
pub enum PlacementError {
    #[error("Anchor atom '{atom_name}' not found in the target residue in the system")]
    AnchorAtomNotFoundInSystem { atom_name: String },

    #[error("Anchor atom '{atom_name}' not found in the rotamer template")]
    AnchorAtomNotFoundInRotamer { atom_name: String },

    #[error(
        "Could not find any atoms with name '{atom_name}' in the rotamer to build the index map"
    )]
    RotamerAtomNameNotFound { atom_name: String },

    #[error(
        "Insufficient anchor atoms for stable alignment: requires at least 3, but found {found}"
    )]
    InsufficientAnchors { found: usize },

    #[error(
        "Placement logic failed for atom '{atom_name}': not enough instances in residue to fulfill placement requirements"
    )]
    InsufficientAtomsInResidue { atom_name: String },
}

#[inline]
pub fn place_rotamer_on_system(
    system: &mut MolecularSystem,
    target_residue_id: ResidueId,
    rotamer: &Rotamer,
    placement_info: &PlacementInfo,
) -> Result<(), EngineError> {
    let result = || -> Result<(), PlacementError> {
        // --- Phase 1: Preparation ---
        let (rotation, translation) =
            calculate_alignment_transform(system, target_residue_id, rotamer, placement_info)?;

        // --- Phase 2: Side-Chain Replacement ---
        remove_old_sidechain(system, target_residue_id, placement_info)?;
        let index_to_id_map = add_new_sidechain_atoms_and_map(
            system,
            target_residue_id,
            rotamer,
            placement_info,
            rotation,
            translation,
        )?;

        // --- Phase 3: Topology Reconstruction ---
        rebuild_topology(system, rotamer, &index_to_id_map)?;

        Ok(())
    }();

    result.map_err(|e| EngineError::Placement {
        residue_id: target_residue_id,
        message: e.to_string(),
    })
}

fn calculate_alignment_transform(
    system: &MolecularSystem,
    target_residue_id: ResidueId,
    rotamer: &Rotamer,
    placement_info: &PlacementInfo,
) -> Result<(Rotation3<f64>, Vector3<f64>), PlacementError> {
    let target_residue = system.residue(target_residue_id).unwrap();
    let mut system_points = Vec::with_capacity(placement_info.anchor_atoms.len());
    let mut rotamer_points = Vec::with_capacity(placement_info.anchor_atoms.len());

    for atom_name in &placement_info.anchor_atoms {
        let system_atom_id = target_residue
            .get_first_atom_id_by_name(atom_name)
            .ok_or_else(|| PlacementError::AnchorAtomNotFoundInSystem {
                atom_name: atom_name.clone(),
            })?;
        system_points.push(system.atom(system_atom_id).unwrap().position);

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

    calculate_transformation(&rotamer_points, &system_points)
}

fn remove_old_sidechain(
    system: &mut MolecularSystem,
    target_residue_id: ResidueId,
    placement_info: &PlacementInfo,
) -> Result<(), PlacementError> {
    let mut frequency_map = HashMap::new();
    for name in &placement_info.sidechain_atoms {
        *frequency_map.entry(name.as_str()).or_insert(0) += 1;
    }

    let mut ids_to_remove = Vec::new();

    {
        let target_residue = system.residue(target_residue_id).unwrap();
        for (name, count) in frequency_map {
            if let Some(atom_ids) = target_residue.get_atom_ids_by_name(name) {
                if atom_ids.len() < count {
                    warn!(
                        "Residue {:?} has only {} atom(s) named '{}', but placement info requires removing {}. This might indicate a malformed input.",
                        target_residue_id,
                        atom_ids.len(),
                        name,
                        count
                    );
                }
                ids_to_remove.extend(atom_ids.iter().rev().take(count));
            }
        }
    }

    for atom_id in ids_to_remove {
        system.remove_atom(atom_id);
    }

    Ok(())
}

fn add_new_sidechain_atoms_and_map(
    system: &mut MolecularSystem,
    target_residue_id: ResidueId,
    rotamer: &Rotamer,
    placement_info: &PlacementInfo,
    rotation: Rotation3<f64>,
    translation: Vector3<f64>,
) -> Result<HashMap<usize, AtomId>, PlacementError> {
    let mut index_to_id_map = HashMap::new();
    let target_residue = system.residue(target_residue_id).unwrap();

    // 1. Pre-populate the map with anchor atoms already in the system.
    for atom_name in &placement_info.anchor_atoms {
        let rotamer_atom_index = rotamer
            .atoms
            .iter()
            .position(|a| &a.name == atom_name)
            .ok_or_else(|| PlacementError::RotamerAtomNameNotFound {
                atom_name: atom_name.clone(),
            })?;

        let system_atom_id = target_residue.get_atom_id_by_name(atom_name).unwrap();
        index_to_id_map.insert(rotamer_atom_index, system_atom_id);
    }

    // 2. Add new side-chain atoms and populate the map.
    for (index, rotamer_atom) in rotamer.atoms.iter().enumerate() {
        if placement_info.sidechain_atoms.contains(&rotamer_atom.name) {
            let mut new_atom = rotamer_atom.clone();
            new_atom.residue_id = target_residue_id;
            new_atom.position = rotation * rotamer_atom.position + translation;

            let new_atom_id = system
                .add_atom_to_residue(target_residue_id, new_atom)
                .unwrap();
            index_to_id_map.insert(index, new_atom_id);
        }
    }

    Ok(index_to_id_map)
}

fn rebuild_topology(
    system: &mut MolecularSystem,
    rotamer: &Rotamer,
    index_to_id_map: &HashMap<usize, AtomId>,
) -> Result<(), PlacementError> {
    for &(index1, index2) in &rotamer.bonds {
        if let (Some(&id1), Some(&id2)) =
            (index_to_id_map.get(&index1), index_to_id_map.get(&index2))
        {
            system.add_bond(id1, id2, BondOrder::Single);
        }
    }
    Ok(())
}

fn calculate_transformation(
    from_points: &[Point3<f64>],
    to_points: &[Point3<f64>],
) -> Result<(Rotation3<f64>, Vector3<f64>), PlacementError> {
    let from_centroid_sum: Vector3<f64> = from_points.iter().map(|p| p.coords).sum();
    let from_centroid = Point3::from(from_centroid_sum / from_points.len() as f64);
    let to_centroid_sum: Vector3<f64> = to_points.iter().map(|p| p.coords).sum();
    let to_centroid = Point3::from(to_centroid_sum / to_points.len() as f64);

    let centered_from: Vec<_> = from_points.iter().map(|p| p - from_centroid).collect();
    let centered_to: Vec<_> = to_points.iter().map(|p| p - to_centroid).collect();

    let h = centered_from
        .iter()
        .zip(centered_to.iter())
        .fold(Matrix3::zeros(), |acc, (f, t)| acc + t * f.transpose());

    let svd = h.svd(true, true);
    let u = svd.u.unwrap();
    let v_t = svd.v_t.unwrap();

    let d = (u * v_t.transpose()).determinant();
    let mut correction = Matrix3::identity();
    if d < 0.0 {
        correction[(2, 2)] = -1.0;
    }

    let rotation_matrix = u * correction * v_t;
    let rotation = Rotation3::from_matrix(&rotation_matrix);
    let translation = to_centroid.coords - rotation * from_centroid.coords;

    Ok((rotation, translation))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{
        models::{atom::Atom, chain::ChainType, residue::ResidueType, topology::BondOrder},
        rotamers::placement::PlacementInfo,
    };
    use nalgebra::{Point3, Vector3};

    fn create_test_system_with_ala_residue() -> (MolecularSystem, ResidueId) {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let residue_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();

        let n = Atom::new("N", residue_id, Point3::new(0.0, 1.0, 0.0));
        let ca = Atom::new("CA", residue_id, Point3::new(0.0, 0.0, 0.0));
        let c = Atom::new("C", residue_id, Point3::new(1.0, 0.0, 0.0));
        let old_cb = Atom::new("CB", residue_id, Point3::new(-1.0, -1.0, -1.0));

        let n_id = system.add_atom_to_residue(residue_id, n).unwrap();
        let ca_id = system.add_atom_to_residue(residue_id, ca).unwrap();
        let c_id = system.add_atom_to_residue(residue_id, c).unwrap();
        let old_cb_id = system.add_atom_to_residue(residue_id, old_cb).unwrap();

        system.add_bond(n_id, ca_id, BondOrder::Single).unwrap();
        system.add_bond(ca_id, c_id, BondOrder::Single).unwrap();
        system
            .add_bond(ca_id, old_cb_id, BondOrder::Single)
            .unwrap();

        (system, residue_id)
    }

    fn create_leu_rotamer() -> Rotamer {
        let atoms = vec![
            Atom::new("N", ResidueId::default(), Point3::new(5.0, 6.0, 5.0)),
            Atom::new("CA", ResidueId::default(), Point3::new(5.0, 5.0, 5.0)),
            Atom::new("C", ResidueId::default(), Point3::new(6.0, 5.0, 5.0)),
            Atom::new("CB", ResidueId::default(), Point3::new(4.0, 4.0, 5.0)),
            Atom::new("CG", ResidueId::default(), Point3::new(3.0, 4.0, 4.0)),
            Atom::new("CD1", ResidueId::default(), Point3::new(2.0, 3.0, 4.0)),
            Atom::new("CD2", ResidueId::default(), Point3::new(2.0, 5.0, 4.0)),
        ];
        let bonds = vec![(0, 1), (1, 2), (1, 3), (3, 4), (4, 5), (4, 6)];
        Rotamer { atoms, bonds }
    }

    fn leu_placement_info() -> PlacementInfo {
        PlacementInfo {
            anchor_atoms: vec!["N".to_string(), "CA".to_string(), "C".to_string()],
            sidechain_atoms: vec![
                "CB".to_string(),
                "CG".to_string(),
                "CD1".to_string(),
                "CD2".to_string(),
            ],
            connection_points: vec![],
            exact_match_atoms: vec![],
        }
    }

    fn bond_exists(system: &MolecularSystem, res_id: ResidueId, name1: &str, name2: &str) -> bool {
        let res = system.residue(res_id).unwrap();
        let id1 = res.get_atom_id_by_name(name1);
        let id2 = res.get_atom_id_by_name(name2);

        if let (Some(id1), Some(id2)) = (id1, id2) {
            system
                .get_bonded_neighbors(id1)
                .map_or(false, |neighbors| neighbors.contains(&id2))
        } else {
            false
        }
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
    fn test_full_placement_workflow() {
        let (mut system, residue_id) = create_test_system_with_ala_residue();
        let rotamer = create_leu_rotamer();
        let placement_info = leu_placement_info();

        let original_atom_count = system.atoms_iter().count();
        assert_eq!(system.residue(residue_id).unwrap().atoms().len(), 4);
        assert!(
            system
                .residue(residue_id)
                .unwrap()
                .get_atom_id_by_name("CG")
                .is_none()
        );
        assert!(bond_exists(&system, residue_id, "CA", "CB"));

        let result = place_rotamer_on_system(&mut system, residue_id, &rotamer, &placement_info);
        assert!(result.is_ok());

        let residue = system.residue(residue_id).unwrap();

        assert_eq!(
            residue.atoms().len(),
            3 + 4,
            "Should have 3 backbone + 4 new sidechain atoms"
        );
        assert_eq!(
            system.atoms_iter().count(),
            original_atom_count - 1 + 4,
            "System total atoms should reflect replacement"
        );

        assert!(residue.get_atom_id_by_name("CB").is_some());
        assert!(residue.get_atom_id_by_name("CG").is_some());
        assert!(residue.get_atom_id_by_name("CD1").is_some());
        assert!(residue.get_atom_id_by_name("CD2").is_some());

        let new_cb_id = residue.get_atom_id_by_name("CB").unwrap();
        let new_cb_atom = system.atom(new_cb_id).unwrap();
        let ca_pos = system
            .atom(residue.get_atom_id_by_name("CA").unwrap())
            .unwrap()
            .position;
        assert!(
            (new_cb_atom.position - ca_pos).norm() < 2.0,
            "New CB should be close to CA"
        );

        assert!(
            bond_exists(&system, residue_id, "CA", "CB"),
            "CA-CB bond must be created"
        );

        assert!(
            bond_exists(&system, residue_id, "CB", "CG"),
            "CB-CG bond must be created"
        );
        assert!(
            bond_exists(&system, residue_id, "CG", "CD1"),
            "CG-CD1 bond must be created"
        );
        assert!(
            bond_exists(&system, residue_id, "CG", "CD2"),
            "CG-CD2 bond must be created"
        );

        assert!(
            bond_exists(&system, residue_id, "N", "CA"),
            "N-CA bond must be preserved"
        );
        assert!(
            bond_exists(&system, residue_id, "CA", "C"),
            "CA-C bond must be preserved"
        );

        assert!(
            !bond_exists(&system, residue_id, "N", "CB"),
            "N-CB bond should not exist"
        );
        assert!(
            !bond_exists(&system, residue_id, "C", "CB"),
            "C-CB bond should not exist"
        );
    }

    #[test]
    fn place_rotamer_fails_if_system_anchor_is_missing() {
        let (mut system, residue_id) = create_test_system_with_ala_residue();
        let n_id = system
            .residue(residue_id)
            .unwrap()
            .get_atom_id_by_name("N")
            .unwrap();
        system.remove_atom(n_id);

        let rotamer = create_leu_rotamer();
        let placement_info = leu_placement_info();

        let result = place_rotamer_on_system(&mut system, residue_id, &rotamer, &placement_info);

        assert!(matches!(result, Err(EngineError::Placement { .. })));
        if let Err(EngineError::Placement { message, .. }) = result {
            assert!(message.contains("Anchor atom 'N' not found in the target residue"));
        }
    }

    #[test]
    fn place_rotamer_fails_if_rotamer_anchor_is_missing() {
        let (mut system, residue_id) = create_test_system_with_ala_residue();
        let mut rotamer = create_leu_rotamer();
        rotamer.atoms.retain(|a| a.name != "N");

        let placement_info = leu_placement_info();

        let result = place_rotamer_on_system(&mut system, residue_id, &rotamer, &placement_info);

        assert!(matches!(result, Err(EngineError::Placement { .. })));
        if let Err(EngineError::Placement { message, .. }) = result {
            assert!(message.contains("Anchor atom 'N' not found in the rotamer template"));
        }
    }

    #[test]
    fn place_rotamer_fails_if_insufficient_anchors() {
        let (mut system, residue_id) = create_test_system_with_ala_residue();
        let rotamer = create_leu_rotamer();
        let mut placement_info = leu_placement_info();
        placement_info.anchor_atoms = vec!["N".to_string(), "CA".to_string()];

        let result = place_rotamer_on_system(&mut system, residue_id, &rotamer, &placement_info);

        assert!(matches!(result, Err(EngineError::Placement { .. })));
        if let Err(EngineError::Placement { message, .. }) = result {
            assert!(message.contains(
                "Insufficient anchor atoms for stable alignment: requires at least 3, but found 2"
            ));
        }
    }
}
