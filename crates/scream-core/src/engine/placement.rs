use super::error::EngineError;
use crate::core::{
    models::{
        ids::{AtomId, ResidueId},
        system::MolecularSystem,
        topology::BondOrder,
    },
    rotamers::rotamer::Rotamer,
    topology::registry::ResidueTopology,
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
        "Placement logic failed for atom '{atom_name}': not enough instances in residue to fulfill placement requirements based on its topology definition"
    )]
    InsufficientAtomsInResidue { atom_name: String },
}

#[inline]
pub fn place_rotamer_on_system(
    system: &mut MolecularSystem,
    target_residue_id: ResidueId,
    rotamer: &Rotamer,
    topology: &ResidueTopology,
) -> Result<(), EngineError> {
    let result = || -> Result<(), PlacementError> {
        // --- Phase 1: Preparation ---
        let (rotation, translation) =
            calculate_alignment_transform(system, target_residue_id, rotamer, topology)?;

        // --- Phase 2: Side-Chain Replacement ---
        remove_old_sidechain(system, target_residue_id, topology)?;
        let index_to_id_map = add_new_sidechain_atoms_and_map(
            system,
            target_residue_id,
            rotamer,
            topology,
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

    // 1. Create a "pool" of rotamer atoms to be consumed during assignment.
    let mut rotamer_atom_pool: HashMap<&str, Vec<(usize, &crate::core::models::atom::Atom)>> =
        HashMap::new();
    for (index, atom) in rotamer.atoms.iter().enumerate() {
        rotamer_atom_pool
            .entry(&atom.name)
            .or_default()
            .push((index, atom));
    }

    // 2. Pre-populate the map with anchor atoms already present in the system.
    for atom_name in &placement_info.anchor_atoms {
        let (rotamer_atom_index, _) = rotamer_atom_pool
            .get_mut(atom_name.as_str())
            .and_then(|atoms| {
                if atoms.is_empty() {
                    None
                } else {
                    Some(atoms.remove(0))
                }
            })
            .ok_or_else(|| PlacementError::RotamerAtomNameNotFound {
                atom_name: atom_name.clone(),
            })?;

        let system_atom_id = target_residue.get_first_atom_id_by_name(atom_name).unwrap();
        index_to_id_map.insert(rotamer_atom_index, system_atom_id);
    }

    // 3. Add new side-chain atoms, consuming them from the end of the rotamer pool.
    for atom_name in &placement_info.sidechain_atoms {
        let (index, rotamer_atom) = rotamer_atom_pool
            .get_mut(atom_name.as_str())
            .and_then(|atoms| atoms.pop())
            .ok_or_else(|| PlacementError::InsufficientAtomsInResidue {
                atom_name: atom_name.clone(),
            })?;

        let mut new_atom = rotamer_atom.clone();
        new_atom.residue_id = target_residue_id;
        new_atom.position = rotation * rotamer_atom.position + translation;

        let new_atom_id = system
            .add_atom_to_residue(target_residue_id, new_atom)
            .unwrap();
        index_to_id_map.insert(index, new_atom_id);
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

    fn create_system_with_residue(
        res_name: &str,
        res_type: ResidueType,
        atoms_to_add: &[(&str, [f64; 3])],
    ) -> (MolecularSystem, ResidueId) {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let residue_id = system
            .add_residue(chain_id, 1, res_name, Some(res_type))
            .unwrap();

        for (name, pos) in atoms_to_add {
            let atom = Atom::new(name, residue_id, Point3::from(*pos));
            system.add_atom_to_residue(residue_id, atom).unwrap();
        }
        (system, residue_id)
    }

    fn create_test_system_with_ala_residue() -> (MolecularSystem, ResidueId) {
        let (mut system, residue_id) = create_system_with_residue(
            "ALA",
            ResidueType::Alanine,
            &[
                ("N", [0.0, 1.0, 0.0]),
                ("CA", [0.0, 0.0, 0.0]),
                ("C", [1.0, 0.0, 0.0]),
                ("CB", [-1.0, -0.5, 0.0]),
                ("HCB", [-1.5, -0.5, 0.8]),
                ("HCB", [-1.5, -0.5, -0.8]),
                ("HCB", [-1.5, 0.5, 0.0]),
            ],
        );
        let n_id = system
            .residue(residue_id)
            .unwrap()
            .get_first_atom_id_by_name("N")
            .unwrap();
        let ca_id = system
            .residue(residue_id)
            .unwrap()
            .get_first_atom_id_by_name("CA")
            .unwrap();
        let c_id = system
            .residue(residue_id)
            .unwrap()
            .get_first_atom_id_by_name("C")
            .unwrap();
        system.add_bond(n_id, ca_id, BondOrder::Single).unwrap();
        system.add_bond(ca_id, c_id, BondOrder::Single).unwrap();
        (system, residue_id)
    }

    fn create_leu_rotamer() -> Rotamer {
        let atoms = vec![
            Atom::new("N", ResidueId::default(), Point3::new(5.0, 6.0, 5.0)),
            Atom::new("CA", ResidueId::default(), Point3::new(5.0, 5.0, 5.0)),
            Atom::new("C", ResidueId::default(), Point3::new(6.0, 5.0, 5.0)),
            Atom::new("CB", ResidueId::default(), Point3::new(4.0, 4.0, 5.0)),
            Atom::new("CG", ResidueId::default(), Point3::new(3.0, 4.0, 4.0)),
        ];
        let bonds = vec![(0, 1), (1, 2), (1, 3), (3, 4)];
        Rotamer { atoms, bonds }
    }

    fn leu_placement_info() -> PlacementInfo {
        PlacementInfo {
            anchor_atoms: vec!["N".to_string(), "CA".to_string(), "C".to_string()],
            sidechain_atoms: vec!["CB".to_string(), "CG".to_string()],
        }
    }

    fn bond_exists(system: &MolecularSystem, res_id: ResidueId, name1: &str, name2: &str) -> bool {
        let res = system.residue(res_id).unwrap();
        let id1 = res.get_first_atom_id_by_name(name1);
        let id2 = res.get_first_atom_id_by_name(name2);

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

        assert!(rot.angle().abs() < 1e-9);
        assert!((trans - translation_vec).norm() < 1e-9);
    }

    #[test]
    fn test_remove_old_sidechain_handles_duplicates() {
        let (mut system, residue_id) = create_test_system_with_ala_residue();
        let placement_info = PlacementInfo {
            anchor_atoms: vec!["N".to_string(), "CA".to_string(), "C".to_string()],
            sidechain_atoms: vec![
                "CB".to_string(),
                "HCB".to_string(),
                "HCB".to_string(),
                "HCB".to_string(),
            ],
        };

        assert_eq!(system.residue(residue_id).unwrap().atoms().len(), 7);
        assert_eq!(
            system
                .residue(residue_id)
                .unwrap()
                .get_atom_ids_by_name("HCB")
                .unwrap()
                .len(),
            3
        );

        remove_old_sidechain(&mut system, residue_id, &placement_info).unwrap();

        let residue = system.residue(residue_id).unwrap();
        assert_eq!(
            residue.atoms().len(),
            3,
            "Only backbone atoms should remain"
        );
        assert!(
            residue.get_atom_ids_by_name("CB").is_none(),
            "CB should be removed"
        );
        assert!(
            residue.get_atom_ids_by_name("HCB").is_none(),
            "All HCB atoms should be removed"
        );
    }

    #[test]
    fn test_full_placement_workflow_with_new_api() {
        let (mut system, residue_id) = create_test_system_with_ala_residue();
        let original_atom_count = system.atoms_iter().count();
        let ala_sidechain_to_remove = PlacementInfo {
            sidechain_atoms: vec![
                "CB".to_string(),
                "HCB".to_string(),
                "HCB".to_string(),
                "HCB".to_string(),
            ],
            anchor_atoms: vec![],
        };
        remove_old_sidechain(&mut system, residue_id, &ala_sidechain_to_remove).unwrap();

        let rotamer = create_leu_rotamer();
        let placement_info = leu_placement_info();

        let result = place_rotamer_on_system(&mut system, residue_id, &rotamer, &placement_info);
        assert!(result.is_ok(), "Placement failed: {:?}", result.err());

        let residue = system.residue(residue_id).unwrap();

        assert_eq!(
            residue.atoms().len(),
            3 + 2,
            "Should have 3 backbone + 2 new LEU sidechain atoms"
        );
        assert_eq!(
            system.atoms_iter().count(),
            original_atom_count - 4 + 2,
            "System total atoms should reflect replacement"
        );

        assert!(residue.get_first_atom_id_by_name("CB").is_some());
        assert!(residue.get_first_atom_id_by_name("CG").is_some());

        assert!(bond_exists(&system, residue_id, "CA", "CB"));
        assert!(bond_exists(&system, residue_id, "CB", "CG"));
        assert!(bond_exists(&system, residue_id, "N", "CA"));
    }

    #[test]
    fn test_placement_on_glycine_special_case() {
        let (mut system, gly_id) = create_system_with_residue(
            "GLY",
            ResidueType::Glycine,
            &[
                ("N", [0.0, 1.0, 0.0]),
                ("CA", [0.0, 0.0, 0.0]),
                ("C", [1.0, 0.0, 0.0]),
                ("HCA", [0.0, -0.5, 0.8]),
                ("HCA", [0.0, -0.5, -0.8]),
            ],
        );

        let gly_placement = PlacementInfo {
            anchor_atoms: vec!["N".to_string(), "CA".to_string(), "HCA".to_string()],
            sidechain_atoms: vec!["HCA".to_string()],
        };

        let gly_rotamer = {
            let atoms = vec![
                Atom::new("N", ResidueId::default(), Point3::new(0.0, 1.0, 0.0)),
                Atom::new("CA", ResidueId::default(), Point3::new(0.0, 0.0, 0.0)),
                Atom::new("HCA", ResidueId::default(), Point3::new(0.0, -0.5, 0.8)),
                Atom::new("HCA", ResidueId::default(), Point3::new(5.0, 5.0, 5.0)),
                Atom::new("C", ResidueId::default(), Point3::new(1.0, 0.0, 0.0)),
            ];
            let bonds = vec![(1, 0), (1, 2), (1, 3), (1, 4)];
            Rotamer { atoms, bonds }
        };

        let result = place_rotamer_on_system(&mut system, gly_id, &gly_rotamer, &gly_placement);
        assert!(result.is_ok(), "GLY placement failed: {:?}", result.err());

        let residue = system.residue(gly_id).unwrap();
        assert_eq!(residue.atoms().len(), 5, "GLY should still have 5 atoms");

        let hca_ids = residue.get_atom_ids_by_name("HCA").unwrap();
        assert_eq!(hca_ids.len(), 2, "Should still have two HCA atoms");

        let hca1 = system.atom(hca_ids[0]).unwrap();
        let hca2 = system.atom(hca_ids[1]).unwrap();

        let new_pos_atom = if (hca1.position - Point3::new(5.0, 5.0, 5.0)).norm() < 1e-6 {
            hca1
        } else {
            hca2
        };
        let old_pos_atom = if new_pos_atom.position == hca1.position {
            hca2
        } else {
            hca1
        };

        assert!(
            (old_pos_atom.position - Point3::new(0.0, -0.5, 0.8)).norm() < 1e-6,
            "The anchor HCA should not have moved"
        );
    }

    #[test]
    fn place_rotamer_fails_if_system_anchor_is_missing() {
        let (mut system, residue_id) = create_test_system_with_ala_residue();
        let n_id = system
            .residue(residue_id)
            .unwrap()
            .get_first_atom_id_by_name("N")
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
            assert!(message.contains("requires at least 3, but found 2"));
        }
    }
}
