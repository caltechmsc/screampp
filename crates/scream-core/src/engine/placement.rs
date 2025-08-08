use super::error::EngineError;
use crate::core::models::atom::AtomRole;
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
        remove_old_sidechain(system, target_residue_id)?;
        let index_to_id_map = add_new_sidechain_atoms_and_map(
            system,
            target_residue_id,
            rotamer,
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
    topology: &ResidueTopology,
) -> Result<(Rotation3<f64>, Vector3<f64>), PlacementError> {
    let system_backbone_atoms: HashMap<_, _> = system
        .residue(target_residue_id)
        .unwrap()
        .atoms()
        .iter()
        .filter_map(|&atom_id| {
            let atom = system.atom(atom_id)?;
            if atom.role == AtomRole::Backbone {
                Some((atom.name.as_str(), atom_id))
            } else {
                None
            }
        })
        .collect();

    let rotamer_backbone_atoms: HashMap<_, _> = rotamer
        .atoms
        .iter()
        .filter(|atom| atom.role == AtomRole::Backbone)
        .map(|atom| (atom.name.as_str(), atom))
        .collect();

    let mut system_points = Vec::with_capacity(topology.anchor_atoms.len());
    let mut rotamer_points = Vec::with_capacity(topology.anchor_atoms.len());

    for atom_name in &topology.anchor_atoms {
        let system_atom_id = *system_backbone_atoms
            .get(atom_name.as_str())
            .ok_or_else(|| PlacementError::AnchorAtomNotFoundInSystem {
                atom_name: atom_name.clone(),
            })?;
        system_points.push(system.atom(system_atom_id).unwrap().position);

        let rotamer_atom = *rotamer_backbone_atoms
            .get(atom_name.as_str())
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
) -> Result<(), PlacementError> {
    let atom_ids_in_residue = system.residue(target_residue_id).unwrap().atoms().to_vec();

    for atom_id in atom_ids_in_residue {
        let should_remove = if let Some(atom) = system.atom(atom_id) {
            atom.role == AtomRole::Sidechain
        } else {
            false
        };

        if should_remove {
            system.remove_atom(atom_id);
        }
    }

    Ok(())
}

fn add_new_sidechain_atoms_and_map(
    system: &mut MolecularSystem,
    target_residue_id: ResidueId,
    rotamer: &Rotamer,
    rotation: Rotation3<f64>,
    translation: Vector3<f64>,
) -> Result<HashMap<usize, AtomId>, PlacementError> {
    // 1. Create a fast lookup map from backbone atom names in the system to their AtomId.
    let system_backbone_map: HashMap<_, _> = system
        .residue(target_residue_id)
        .unwrap()
        .atoms()
        .iter()
        .filter_map(|&atom_id| {
            let atom = system.atom(atom_id)?;
            if atom.role == AtomRole::Backbone {
                Some((atom.name.clone(), atom_id))
            } else {
                None
            }
        })
        .collect();

    let mut index_to_id_map = HashMap::new();

    // 2. Iterate over all atoms in the rotamer.
    for (rotamer_index, rotamer_atom) in rotamer.atoms.iter().enumerate() {
        match rotamer_atom.role {
            // 2a. If it is a backbone atom, find the corresponding existing atom in the system and add it to the map.
            AtomRole::Backbone => {
                if let Some(&system_atom_id) = system_backbone_map.get(&rotamer_atom.name) {
                    index_to_id_map.insert(rotamer_index, system_atom_id);
                } else {
                    return Err(PlacementError::AnchorAtomNotFoundInSystem {
                        atom_name: rotamer_atom.name.clone(),
                    });
                }
            }
            // 2b. If it is a sidechain atom, apply the transformation, create a new atom, and add it to the system and the map.
            AtomRole::Sidechain => {
                let mut new_atom = rotamer_atom.clone();
                new_atom.residue_id = target_residue_id;
                new_atom.position = rotation * rotamer_atom.position + translation;

                let new_atom_id = system
                    .add_atom_to_residue(target_residue_id, new_atom)
                    .unwrap();

                index_to_id_map.insert(rotamer_index, new_atom_id);
            }
            _ => {
                // Do nothing for Ligand, Water, Other roles in this context.
            }
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
        models::{atom::Atom, chain::ChainType, residue::ResidueType},
        topology::registry::ResidueTopology,
    };

    struct TestSetup {
        system: MolecularSystem,
        residue_id: ResidueId,
        topology: ResidueTopology,
    }

    fn setup_ala_system() -> TestSetup {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let residue_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();

        let atoms_to_add = vec![
            ("N", AtomRole::Backbone, [0.0, 1.4, 0.0]),
            ("CA", AtomRole::Backbone, [0.0, 0.0, 0.0]),
            ("C", AtomRole::Backbone, [1.4, 0.0, 0.0]),
            ("CB", AtomRole::Sidechain, [0.0, -0.7, 1.2]),
            ("HB", AtomRole::Sidechain, [0.5, -1.5, 1.5]),
        ];

        for (name, role, pos) in atoms_to_add {
            let mut atom = Atom::new(name, residue_id, Point3::from(pos));
            atom.role = role;
            system.add_atom_to_residue(residue_id, atom).unwrap();
        }

        let topology = ResidueTopology {
            anchor_atoms: vec!["N".to_string(), "CA".to_string(), "C".to_string()],
            sidechain_atoms: vec!["CB".to_string(), "HB".to_string()],
        };

        TestSetup {
            system,
            residue_id,
            topology,
        }
    }

    fn create_leu_rotamer() -> Rotamer {
        let mut atoms = vec![
            Atom::new("N", ResidueId::default(), Point3::new(5.0, 6.4, 5.0)),
            Atom::new("CA", ResidueId::default(), Point3::new(5.0, 5.0, 5.0)),
            Atom::new("C", ResidueId::default(), Point3::new(6.4, 5.0, 5.0)),
            Atom::new("CB", ResidueId::default(), Point3::new(5.0, 4.3, 3.8)),
            Atom::new("CG", ResidueId::default(), Point3::new(4.3, 3.1, 3.5)),
        ];

        for atom in &mut atoms {
            atom.role = if ["N", "CA", "C"].contains(&atom.name.as_str()) {
                AtomRole::Backbone
            } else {
                AtomRole::Sidechain
            };
        }

        let bonds = vec![(0, 1), (1, 2), (1, 3), (3, 4)];
        Rotamer { atoms, bonds }
    }

    fn bond_exists(system: &MolecularSystem, res_id: ResidueId, name1: &str, name2: &str) -> bool {
        let res = system.residue(res_id).unwrap();
        if let (Some(id1), Some(id2)) = (
            res.get_first_atom_id_by_name(name1),
            res.get_first_atom_id_by_name(name2),
        ) {
            system
                .get_bonded_neighbors(id1)
                .map_or(false, |neighbors| neighbors.contains(&id2))
        } else {
            false
        }
    }

    #[test]
    fn remove_old_sidechain_correctly_removes_only_sidechain_atoms() {
        let TestSetup {
            mut system,
            residue_id,
            ..
        } = setup_ala_system();
        assert_eq!(system.residue(residue_id).unwrap().atoms().len(), 5);

        let result = remove_old_sidechain(&mut system, residue_id);
        assert!(result.is_ok());

        let residue = system.residue(residue_id).unwrap();
        assert_eq!(residue.atoms().len(), 3);
        assert!(residue.get_first_atom_id_by_name("N").is_some());
        assert!(residue.get_first_atom_id_by_name("CA").is_some());
        assert!(residue.get_first_atom_id_by_name("C").is_some());
        assert!(residue.get_first_atom_id_by_name("CB").is_none());
        assert!(residue.get_first_atom_id_by_name("HB").is_none());
    }

    #[test]
    fn calculate_transformation_pure_translation() {
        let from = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let to = vec![
            Point3::new(10.0, 20.0, 30.0),
            Point3::new(11.0, 20.0, 30.0),
            Point3::new(10.0, 21.0, 30.0),
        ];

        let (rot, trans) = calculate_transformation(&from, &to).unwrap();

        assert!(
            (rot.angle()).abs() < 1e-9,
            "Rotation should be near zero for pure translation"
        );
        assert!(
            (trans - Vector3::new(10.0, 20.0, 30.0)).norm() < 1e-9,
            "Translation vector is incorrect"
        );
    }

    #[test]
    fn add_new_sidechain_atoms_maps_backbone_and_adds_sidechain() {
        let TestSetup {
            mut system,
            residue_id,
            ..
        } = setup_ala_system();
        remove_old_sidechain(&mut system, residue_id).unwrap();

        let leu_rotamer = create_leu_rotamer();
        let (rotation, translation) = (Rotation3::identity(), Vector3::zeros());

        let map = add_new_sidechain_atoms_and_map(
            &mut system,
            residue_id,
            &leu_rotamer,
            rotation,
            translation,
        )
        .unwrap();

        assert_eq!(
            map.len(),
            5,
            "Map should contain all 5 atoms from the rotamer"
        );
        assert_eq!(
            system.residue(residue_id).unwrap().atoms().len(),
            5,
            "System residue should now have 5 atoms"
        );

        let cg_atom_id = system
            .residue(residue_id)
            .unwrap()
            .get_first_atom_id_by_name("CG");
        assert!(cg_atom_id.is_some());

        let ca_id_in_system = system
            .residue(residue_id)
            .unwrap()
            .get_first_atom_id_by_name("CA")
            .unwrap();
        let ca_index_in_rotamer = leu_rotamer
            .atoms
            .iter()
            .position(|a| a.name == "CA")
            .unwrap();
        assert_eq!(map.get(&ca_index_in_rotamer), Some(&ca_id_in_system));
    }

    #[test]
    fn rebuild_topology_creates_correct_bonds() {
        let TestSetup {
            mut system,
            residue_id,
            ..
        } = setup_ala_system();
        remove_old_sidechain(&mut system, residue_id).unwrap();

        let leu_rotamer = create_leu_rotamer();
        let (rotation, translation) = (Rotation3::identity(), Vector3::zeros());
        let map = add_new_sidechain_atoms_and_map(
            &mut system,
            residue_id,
            &leu_rotamer,
            rotation,
            translation,
        )
        .unwrap();

        rebuild_topology(&mut system, &leu_rotamer, &map).unwrap();

        assert!(
            bond_exists(&system, residue_id, "CA", "CB"),
            "CA-CB bond missing"
        );
        assert!(
            bond_exists(&system, residue_id, "CB", "CG"),
            "CB-CG bond missing"
        );
        assert!(
            !bond_exists(&system, residue_id, "N", "C"),
            "N-C bond should not exist"
        );
    }

    #[test]
    fn place_rotamer_on_system_full_workflow() {
        let TestSetup {
            mut system,
            residue_id,
            topology,
        } = setup_ala_system();
        let original_atom_count = system.atoms_iter().count();

        let leu_rotamer = create_leu_rotamer();

        let result = place_rotamer_on_system(&mut system, residue_id, &leu_rotamer, &topology);
        assert!(result.is_ok(), "Placement failed: {:?}", result.err());

        let residue = system.residue(residue_id).unwrap();
        assert_eq!(
            residue.atoms().len(),
            5,
            "Residue should have 3 backbone + 2 LEU sidechain atoms"
        );
        assert_eq!(
            system.atoms_iter().count(),
            original_atom_count - 2 + 2,
            "System atom count should be updated"
        );
        assert!(
            residue.get_first_atom_id_by_name("CG").is_some(),
            "CG from LEU should be present"
        );
        assert!(
            residue.get_first_atom_id_by_name("HB").is_none(),
            "HB from ALA should be gone"
        );
        assert!(bond_exists(&system, residue_id, "CA", "CB"));
    }

    #[test]
    fn place_rotamer_fails_if_system_anchor_is_missing() {
        let TestSetup {
            mut system,
            residue_id,
            topology,
        } = setup_ala_system();
        let n_id = system
            .residue(residue_id)
            .unwrap()
            .get_first_atom_id_by_name("N")
            .unwrap();
        system.remove_atom(n_id);

        let leu_rotamer = create_leu_rotamer();
        let result = place_rotamer_on_system(&mut system, residue_id, &leu_rotamer, &topology);

        assert!(matches!(result, Err(EngineError::Placement { .. })));
        if let Err(EngineError::Placement { message, .. }) = result {
            assert!(message.contains("Anchor atom 'N' not found in the target residue"));
        }
    }

    #[test]
    fn place_rotamer_fails_if_rotamer_anchor_is_missing() {
        let TestSetup {
            mut system,
            residue_id,
            topology,
        } = setup_ala_system();
        let mut leu_rotamer = create_leu_rotamer();
        leu_rotamer.atoms.retain(|a| a.name != "CA");

        let result = place_rotamer_on_system(&mut system, residue_id, &leu_rotamer, &topology);

        assert!(matches!(result, Err(EngineError::Placement { .. })));
        if let Err(EngineError::Placement { message, .. }) = result {
            assert!(message.contains("Anchor atom 'CA' not found in the rotamer template"));
        }
    }

    #[test]
    fn place_rotamer_fails_if_insufficient_anchors() {
        let TestSetup {
            mut system,
            residue_id,
            mut topology,
        } = setup_ala_system();
        topology.anchor_atoms = vec!["N".to_string(), "CA".to_string()];

        let leu_rotamer = create_leu_rotamer();
        let result = place_rotamer_on_system(&mut system, residue_id, &leu_rotamer, &topology);

        assert!(matches!(result, Err(EngineError::Placement { .. })));
        if let Err(EngineError::Placement { message, .. }) = result {
            assert!(message.contains("requires at least 3, but found 2"));
        }
    }
}
