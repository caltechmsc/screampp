use crate::core::models::chain::ChainType;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::core::models::system::MolecularSystem;
use crate::core::rotamers::library::RotamerLibrary;
use crate::engine::config::ResidueSelection;
use crate::engine::error::EngineError;
use kiddo::{KdTree, SquaredEuclidean};
use std::collections::{HashMap, HashSet};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Collects all sidechain atoms for a set of active residues.
///
/// This function iterates through the specified active residues and extracts
/// all atoms classified as sidechain atoms according to their role in the molecular system.
/// The result is a mapping from residue IDs to their corresponding sidechain atom IDs.
///
/// # Arguments
///
/// * `system` - The molecular system containing the residues and atoms.
/// * `active_residues` - A set of residue IDs for which to collect sidechain atoms.
///
/// # Return
///
/// Returns a `HashMap` where keys are residue IDs and values are vectors of sidechain atom IDs.
pub fn collect_active_sidechain_atoms(
    system: &MolecularSystem,
    active_residues: &HashSet<ResidueId>,
) -> HashMap<ResidueId, Vec<AtomId>> {
    let mut map = HashMap::with_capacity(active_residues.len());
    for &residue_id in active_residues {
        if let Some(residue) = system.residue(residue_id) {
            let sidechain_ids: Vec<AtomId> = residue
                .atoms()
                .iter()
                .filter_map(|&atom_id| {
                    system.atom(atom_id).and_then(|atom| {
                        if atom.role == crate::core::models::atom::AtomRole::Sidechain {
                            Some(atom_id)
                        } else {
                            None
                        }
                    })
                })
                .collect();
            map.insert(residue_id, sidechain_ids);
        }
    }
    map
}

/// Precomputes the set of environment atoms for optimization.
///
/// Environment atoms are those that are not part of the active residues' sidechains
/// but may interact with them during energy calculations. This includes backbone atoms
/// from active residues and all atoms from inactive residues.
///
/// # Arguments
///
/// * `system` - The molecular system containing all atoms.
/// * `active_residues` - A set of residue IDs that are being optimized.
///
/// # Return
///
/// Returns a vector of atom IDs representing the environment atoms.
pub fn precompute_environment_atoms(
    system: &MolecularSystem,
    active_residues: &HashSet<ResidueId>,
) -> Vec<AtomId> {
    system
        .atoms_iter()
        .filter_map(|(atom_id, atom)| {
            if !active_residues.contains(&atom.residue_id)
                || atom.role == crate::core::models::atom::AtomRole::Backbone
            {
                Some(atom_id)
            } else {
                None
            }
        })
        .collect()
}

/// Resolves a residue selection specification into a set of residue IDs.
///
/// This function interprets various types of residue selection criteria and converts
/// them into concrete residue IDs from the molecular system. It supports selecting
/// all residues, explicit include/exclude lists, and ligand binding site detection.
/// The final result is filtered to only include residues that have available rotamers
/// in the provided library.
///
/// # Arguments
///
/// * `system` - The molecular system containing the residues to select from.
/// * `selection` - The selection criteria specifying which residues to include.
/// * `library` - The rotamer library used to filter residues by rotamer availability.
///
/// # Return
///
/// Returns a set of residue IDs that match the selection criteria and have available rotamers.
///
/// # Errors
///
/// Returns `EngineError::ResidueNotFound` if a specified residue cannot be found in the system.
pub fn resolve_selection_to_ids(
    system: &MolecularSystem,
    selection: &ResidueSelection,
    library: &RotamerLibrary,
) -> Result<HashSet<ResidueId>, EngineError> {
    let mut candidate_ids: HashSet<ResidueId> = HashSet::new();

    match selection {
        ResidueSelection::All => {
            // Select all residues in the system
            candidate_ids = system.residues_iter().map(|(id, _)| id).collect();
        }
        ResidueSelection::List { include, exclude } => {
            // Handle explicit include/exclude lists
            if include.is_empty() && !exclude.is_empty() {
                // If only exclusions specified, start with all residues
                candidate_ids = system.residues_iter().map(|(id, _)| id).collect();
            } else {
                // Add explicitly included residues
                for spec in include {
                    let chain_id = system
                        .find_chain_by_id(spec.chain_id)
                        .ok_or_else(|| EngineError::ResidueNotFound { spec: spec.clone() })?;
                    let residue_id = system
                        .find_residue_by_id(chain_id, spec.residue_number)
                        .ok_or_else(|| EngineError::ResidueNotFound { spec: spec.clone() })?;
                    candidate_ids.insert(residue_id);
                }
            }

            // Remove explicitly excluded residues
            for spec in exclude {
                if let Some(chain_id) = system.find_chain_by_id(spec.chain_id) {
                    if let Some(residue_id) =
                        system.find_residue_by_id(chain_id, spec.residue_number)
                    {
                        candidate_ids.remove(&residue_id);
                    }
                }
            }
        }
        ResidueSelection::LigandBindingSite {
            ligand_residue,
            radius_angstroms,
        } => {
            // Find residues within binding distance of the ligand
            let ligand_chain_id = system
                .find_chain_by_id(ligand_residue.chain_id)
                .ok_or_else(|| EngineError::ResidueNotFound {
                    spec: ligand_residue.clone(),
                })?;
            let ligand_res_id = system
                .find_residue_by_id(ligand_chain_id, ligand_residue.residue_number)
                .ok_or_else(|| EngineError::ResidueNotFound {
                    spec: ligand_residue.clone(),
                })?;

            // Collect heavy atom positions from the ligand
            let mut ligand_heavy_atom_positions: Vec<[f64; 3]> = Vec::new();
            if let Some(ligand_res) = system.residue(ligand_res_id) {
                for atom_id in ligand_res.atoms() {
                    if let Some(atom) = system.atom(*atom_id) {
                        if is_heavy_atom(&atom.name) {
                            ligand_heavy_atom_positions.push([
                                atom.position.x,
                                atom.position.y,
                                atom.position.z,
                            ]);
                        }
                    }
                }
            }

            if ligand_heavy_atom_positions.is_empty() {
                return Ok(HashSet::new());
            }

            // Build spatial index for efficient distance queries
            let kdtree: KdTree<f64, 3> = (&ligand_heavy_atom_positions).into();
            let radius_sq = radius_angstroms * radius_angstroms;

            // Use parallel iteration if available
            #[cfg(not(feature = "parallel"))]
            let iterator = system.residues_iter();

            #[cfg(feature = "parallel")]
            let iterator = system.residues_iter().par_bridge();

            // Find protein residues within the specified radius
            let binding_site_ids: HashSet<ResidueId> = iterator
                .filter_map(|(res_id, residue)| {
                    if res_id == ligand_res_id {
                        return None;
                    }
                    if let Some(chain) = system.chain(residue.chain_id) {
                        if chain.chain_type != ChainType::Protein {
                            return None;
                        }
                    } else {
                        return None;
                    }

                    // Check if any heavy atom of this residue is within radius
                    let is_in_binding_site = residue.atoms().iter().any(|protein_atom_id| {
                        if let Some(protein_atom) = system.atom(*protein_atom_id) {
                            if is_heavy_atom(&protein_atom.name) {
                                let protein_pos = [
                                    protein_atom.position.x,
                                    protein_atom.position.y,
                                    protein_atom.position.z,
                                ];
                                let nearest = kdtree.nearest_one::<SquaredEuclidean>(&protein_pos);
                                return nearest.distance <= radius_sq;
                            }
                        }
                        false
                    });

                    if is_in_binding_site {
                        Some(res_id)
                    } else {
                        None
                    }
                })
                .collect();

            candidate_ids.extend(binding_site_ids);
        }
    };

    // Filter to only include residues with available rotamers
    let final_active_residues = candidate_ids
        .into_iter()
        .filter(|&residue_id| {
            system
                .residue(residue_id)
                .and_then(|res| res.residue_type)
                .map_or(false, |residue_type| {
                    library.get_rotamers_for(residue_type).is_some()
                })
        })
        .collect();

    Ok(final_active_residues)
}

/// Determines whether an atom is a heavy atom (non-hydrogen).
///
/// Heavy atoms are defined as any atom whose name does not start with 'H' or 'D'
/// (deuterium). This classification is used in various molecular calculations
/// where hydrogen atoms are often treated differently due to their small size
/// and different interaction properties.
///
/// # Arguments
///
/// * `atom_name` - The name of the atom to classify.
///
/// # Return
///
/// Returns `true` if the atom is a heavy atom, `false` if it is hydrogen or deuterium.
fn is_heavy_atom(atom_name: &str) -> bool {
    let first_char = atom_name
        .trim()
        .chars()
        .next()
        .map(|c| c.to_ascii_uppercase());
    !matches!(first_char, Some('H') | Some('D'))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{
        forcefield::{
            parameterization::Parameterizer,
            params::{Forcefield, GlobalParams, NonBondedParams, VdwParam},
        },
        models::{atom::Atom, chain::ChainType, residue::ResidueType, system::MolecularSystem},
        rotamers::{library::RotamerLibrary, rotamer::Rotamer},
    };
    use crate::engine::config::{ResidueSelection, ResidueSpecifier};
    use nalgebra::Point3;
    use std::{
        collections::{HashMap, HashSet},
        fs,
        path::Path,
    };
    use tempfile::TempDir;

    struct TestSetup {
        system: MolecularSystem,
        rotamer_library: RotamerLibrary,
        _temp_dir: TempDir,
    }

    fn setup_test_data() -> TestSetup {
        let temp_dir = tempfile::tempdir().expect("Failed to create temp dir");
        let dir_path = temp_dir.path();

        let topology_path = dir_path.join("topology.toml");
        write_file(
            &topology_path,
            r#"
            [ALA]
            anchor_atoms = ["N", "CA", "C"]
            sidechain_atoms = ["CB"]

            [LEU]
            anchor_atoms = ["N", "CA", "C"]
            sidechain_atoms = ["CB"]
            "#,
        );

        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_a_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let res_b_id = system
            .add_residue(chain_id, 2, "LEU", Some(ResidueType::Leucine))
            .unwrap();

        let add_residue_atoms = |system: &mut MolecularSystem, res_id: ResidueId, offset: f64| {
            let backbone_atoms_data = vec![
                ("N", Point3::new(offset, 1.0, 0.0)),
                ("CA", Point3::new(offset, 0.0, 0.0)),
                ("C", Point3::new(offset + 1.0, 0.0, 0.0)),
            ];
            for (name, pos) in backbone_atoms_data {
                let mut atom = Atom::new(name, res_id, pos);
                atom.force_field_type = "BB".to_string();
                system.add_atom_to_residue(res_id, atom).unwrap();
            }
            let mut cb_atom = Atom::new("CB", res_id, Point3::new(offset, -0.5, 1.2));
            cb_atom.force_field_type = "C_SC".to_string();
            system.add_atom_to_residue(res_id, cb_atom).unwrap();
        };

        add_residue_atoms(&mut system, res_a_id, 0.0);
        add_residue_atoms(&mut system, res_b_id, 2.0);

        let mut vdw = HashMap::new();
        vdw.insert(
            "BB".to_string(),
            VdwParam::LennardJones {
                radius: 1.0,
                well_depth: 0.0,
            },
        );
        vdw.insert(
            "C_SC".to_string(),
            VdwParam::LennardJones {
                radius: 3.8,
                well_depth: 0.1,
            },
        );
        let forcefield = Forcefield {
            non_bonded: NonBondedParams {
                globals: GlobalParams {
                    dielectric_constant: 1.0,
                    potential_function: "lennard-jones-12-6".to_string(),
                },
                vdw,
                hbond: HashMap::new(),
                hbond_donors: HashSet::new(),
                hbond_acceptors: HashSet::new(),
            },
            deltas: HashMap::new(),
            weight_map: HashMap::new(),
        };

        let topology_registry =
            crate::core::topology::registry::TopologyRegistry::load(&topology_path).unwrap();
        let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 0.0);

        let create_rotamer = |residue_id, cb_pos: Point3<f64>| -> Rotamer {
            let placeholder_residue_id = ResidueId::default();
            let atoms = vec![
                Atom::new("N", placeholder_residue_id, Point3::new(0.0, 1.0, 0.0)),
                Atom::new("CA", placeholder_residue_id, Point3::new(0.0, 0.0, 0.0)),
                Atom::new("C", placeholder_residue_id, Point3::new(1.0, 0.0, 0.0)),
                Atom::new("CB", placeholder_residue_id, cb_pos),
            ];

            let mut rotamer = Rotamer {
                atoms,
                bonds: vec![(0, 1), (1, 2), (1, 3)],
            };

            rotamer.atoms.iter_mut().for_each(|a| {
                a.force_field_type = if a.name == "CB" {
                    "C_SC".to_string()
                } else {
                    "BB".to_string()
                };
            });

            let res_name = system.residue(residue_id).unwrap().name.as_str();
            let topo = topology_registry.get(res_name).unwrap();

            parameterizer
                .parameterize_rotamer(&mut rotamer, res_name, topo)
                .unwrap();

            rotamer
        };

        let rotamer_a0 = create_rotamer(res_a_id, Point3::new(-0.5, -0.8, 0.0));
        let rotamer_b0 = create_rotamer(res_b_id, Point3::new(0.5, 0.8, 0.0));

        let mut rot_lib_map = HashMap::new();
        rot_lib_map.insert(ResidueType::Alanine, vec![rotamer_a0]);
        rot_lib_map.insert(ResidueType::Leucine, vec![rotamer_b0]);

        let rotamer_library = RotamerLibrary {
            rotamers: rot_lib_map,
        };

        parameterizer.parameterize_system(&mut system).unwrap();

        TestSetup {
            system,
            rotamer_library,
            _temp_dir: temp_dir,
        }
    }

    fn write_file(path: &Path, content: &str) {
        fs::write(path, content).expect("Failed to write temporary file for test");
    }

    #[test]
    fn collect_active_sidechain_atoms_works_correctly() {
        let setup = setup_test_data();
        let active_residues: HashSet<ResidueId> =
            setup.system.residues_iter().map(|(id, _)| id).collect();

        let map = collect_active_sidechain_atoms(&setup.system, &active_residues);

        assert_eq!(map.len(), 2);

        let ala_res_id = setup
            .system
            .residues_iter()
            .find(|(_, res)| res.name == "ALA")
            .unwrap()
            .0;
        assert_eq!(map.get(&ala_res_id).unwrap().len(), 1);

        let leu_res_id = setup
            .system
            .residues_iter()
            .find(|(_, res)| res.name == "LEU")
            .unwrap()
            .0;
        assert_eq!(map.get(&leu_res_id).unwrap().len(), 1);

        for (res_id, atoms) in &map {
            for &atom_id in atoms {
                let atom = setup.system.atom(atom_id).unwrap();
                assert_eq!(atom.role, crate::core::models::atom::AtomRole::Sidechain);
                assert_eq!(atom.residue_id, *res_id);
            }
        }
    }

    #[test]
    fn precompute_environment_atoms_works_correctly() {
        let setup = setup_test_data();
        let active_residues: HashSet<ResidueId> =
            setup.system.residues_iter().map(|(id, _)| id).collect();

        let env_atoms = precompute_environment_atoms(&setup.system, &active_residues);

        let expected_backbone_count = setup
            .system
            .residues_iter()
            .map(|(_, res)| {
                res.atoms()
                    .iter()
                    .filter(|&&atom_id| {
                        setup.system.atom(atom_id).unwrap().role
                            == crate::core::models::atom::AtomRole::Backbone
                    })
                    .count()
            })
            .sum::<usize>();

        assert_eq!(env_atoms.len(), expected_backbone_count);

        for &atom_id in &env_atoms {
            let atom = setup.system.atom(atom_id).unwrap();
            if active_residues.contains(&atom.residue_id) {
                assert_eq!(atom.role, crate::core::models::atom::AtomRole::Backbone);
            }
        }
    }

    #[test]
    fn resolve_selection_to_ids_handles_all_selection() {
        let setup = setup_test_data();
        let selection = ResidueSelection::All;

        let result =
            resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library).unwrap();

        assert_eq!(result.len(), 2);
        for &res_id in &result {
            let residue = setup.system.residue(res_id).unwrap();
            let res_type = residue.residue_type.unwrap();
            assert!(setup.rotamer_library.get_rotamers_for(res_type).is_some());
        }
    }

    #[test]
    fn resolve_selection_to_ids_handles_list_selection() {
        let setup = setup_test_data();
        let residues: Vec<_> = setup.system.residues_iter().collect();
        let first_res = residues[0].1;
        let second_res = residues[1].1;

        let selection = ResidueSelection::List {
            include: vec![
                ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: first_res.residue_number,
                },
                ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: second_res.residue_number,
                },
            ],
            exclude: vec![],
        };

        let result =
            resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library).unwrap();

        assert!(!result.is_empty());
        let first_res_id = residues[0].0;
        let second_res_id = residues[1].0;
        assert!(result.contains(&first_res_id) || result.contains(&second_res_id));
    }

    #[test]
    fn is_heavy_atom_identifies_correctly() {
        assert!(!is_heavy_atom("H"));
        assert!(!is_heavy_atom("H1"));
        assert!(!is_heavy_atom("D"));
        assert!(!is_heavy_atom("D1"));
        assert!(is_heavy_atom("C"));
        assert!(is_heavy_atom("N"));
        assert!(is_heavy_atom("O"));
        assert!(is_heavy_atom("S"));
        assert!(is_heavy_atom("CA"));
        assert!(is_heavy_atom("CB"));
    }

    #[test]
    fn resolve_selection_to_ids_filters_by_rotamer_availability() {
        let setup = setup_test_data();

        let selection = ResidueSelection::All;
        let result =
            resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library).unwrap();

        for &res_id in &result {
            let residue = setup.system.residue(res_id).unwrap();
            let res_type = residue.residue_type.unwrap();
            assert!(setup.rotamer_library.get_rotamers_for(res_type).is_some());
        }
    }

    #[test]
    fn resolve_selection_to_ids_handles_list_selection_with_include_and_exclude() {
        let setup = setup_test_data();
        let residues: Vec<_> = setup.system.residues_iter().collect();
        let first_res = residues[0].1;
        let second_res = residues[1].1;

        let selection = ResidueSelection::List {
            include: vec![
                ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: first_res.residue_number,
                },
                ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: second_res.residue_number,
                },
            ],
            exclude: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: first_res.residue_number,
            }],
        };

        let result =
            resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library).unwrap();

        assert_eq!(result.len(), 1);
        let expected_res_id = residues
            .iter()
            .find(|(_, res)| res.residue_number == second_res.residue_number)
            .unwrap()
            .0;
        assert!(result.contains(&expected_res_id));
    }

    #[test]
    fn resolve_selection_to_ids_handles_list_selection_with_exclude_only() {
        let setup = setup_test_data();
        let residues: Vec<_> = setup.system.residues_iter().collect();
        let first_res = residues[0].1;

        let selection = ResidueSelection::List {
            include: vec![],
            exclude: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: first_res.residue_number,
            }],
        };

        let result =
            resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library).unwrap();

        assert_eq!(result.len(), residues.len() - 1);
        let excluded_res_id = residues[0].0;
        assert!(!result.contains(&excluded_res_id));
    }

    #[test]
    fn resolve_selection_to_ids_fails_for_nonexistent_residue() {
        let setup = setup_test_data();
        let selection = ResidueSelection::List {
            include: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 999,
            }],
            exclude: vec![],
        };

        let result = resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library);
        assert!(result.is_err());
    }

    #[test]
    fn resolve_ligand_binding_site_selection() {
        let setup = setup_test_data();
        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue: ResidueSpecifier {
                chain_id: 'A',
                residue_number: 1,
            },
            radius_angstroms: 6.0,
        };

        let result =
            resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library).unwrap();

        assert!(!result.is_empty());
    }

    #[test]
    fn resolve_ligand_binding_site_selection_with_larger_radius() {
        let setup = setup_test_data();
        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue: ResidueSpecifier {
                chain_id: 'A',
                residue_number: 1,
            },
            radius_angstroms: 11.0,
        };

        let result =
            resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library).unwrap();

        assert!(!result.is_empty());
    }

    #[test]
    fn resolve_ligand_binding_site_selection_empty_when_no_protein_nearby() {
        let setup = setup_test_data();
        let ligand_spec = ResidueSpecifier {
            chain_id: 'B',
            residue_number: 999,
        };
        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue: ligand_spec,
            radius_angstroms: 10.0,
        };

        let result = resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library);
        assert!(result.is_err());
    }
}
