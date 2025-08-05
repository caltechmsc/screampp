use super::placement::PlacementInfo;
use super::rotamer::{Rotamer, RotamerData};
use crate::core::forcefield::parameterization::{ParameterizationError, Parameterizer};
use crate::core::forcefield::params::Forcefield;
use crate::core::models::atom::Atom;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::core::models::residue::ResidueType;
use crate::core::models::system::MolecularSystem;
use nalgebra::Point3;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::str::FromStr;
use thiserror::Error;

type RawRotamerFile = HashMap<String, Vec<RotamerData>>;

#[derive(Debug, Default, Clone)]
pub struct RotamerLibrary {
    pub rotamers: HashMap<ResidueType, Vec<Rotamer>>,
    pub placement_info: HashMap<ResidueType, PlacementInfo>,
}

#[derive(Debug, Error)]
pub enum LibraryLoadError {
    #[error("File I/O error for '{path}': {source}")]
    Io {
        path: String,
        source: std::io::Error,
    },
    #[error("TOML parsing error for '{path}': {source}")]
    Toml {
        path: String,
        source: toml::de::Error,
    },
    #[error("Unknown residue type '{0}' found in library file")]
    UnknownResidueType(String),
    #[error(
        "Missing placement info for residue type '{0}' which is present in the rotamer library"
    )]
    MissingPlacementInfo(String),
    #[error(
        "Parameterization failed for residue '{residue_type}' in rotamer from file '{path}': {source}"
    )]
    Parameterization {
        path: String,
        residue_type: String,
        source: ParameterizationError,
    },
    #[error(
        "Invalid bond definition in rotamer library for residue '{residue_type}': bond references non-existent atom serial '{serial}'"
    )]
    InvalidBondSerial { residue_type: String, serial: usize },
    #[error(
        "Duplicate atom serial '{serial}' found in rotamer definition for residue '{residue_type}'"
    )]
    DuplicateAtomSerial { residue_type: String, serial: usize },
}

impl RotamerLibrary {
    pub fn load(
        rotamer_toml_path: &Path,
        placement_registry_path: &Path,
        forcefield: &Forcefield,
        delta_s_factor: f64,
    ) -> Result<Self, LibraryLoadError> {
        // --- Phase 1: Load auxiliary placement info ---
        let placement_content =
            std::fs::read_to_string(placement_registry_path).map_err(|e| LibraryLoadError::Io {
                path: placement_registry_path.to_string_lossy().to_string(),
                source: e,
            })?;
        let placement_registry =
            toml::from_str::<HashMap<String, PlacementInfo>>(&placement_content).map_err(|e| {
                LibraryLoadError::Toml {
                    path: placement_registry_path.to_string_lossy().to_string(),
                    source: e,
                }
            })?;

        // --- Phase 2: Load raw rotamer data from TOML ---
        let content =
            std::fs::read_to_string(rotamer_toml_path).map_err(|e| LibraryLoadError::Io {
                path: rotamer_toml_path.to_string_lossy().to_string(),
                source: e,
            })?;
        let raw_lib: RawRotamerFile =
            toml::from_str(&content).map_err(|e| LibraryLoadError::Toml {
                path: rotamer_toml_path.to_string_lossy().to_string(),
                source: e,
            })?;

        // --- Phase 3: Process raw data into functional Rotamer objects ---
        let parameterizer = Parameterizer::new(forcefield.clone(), delta_s_factor);
        let mut final_rotamers_map = HashMap::new();
        let mut final_placement_map = HashMap::new();

        for (res_name, raw_rotamer_list) in raw_lib {
            let residue_type = ResidueType::from_str(&res_name)
                .map_err(|_| LibraryLoadError::UnknownResidueType(res_name.clone()))?;

            let placement_data = placement_registry
                .get(&res_name)
                .ok_or_else(|| LibraryLoadError::MissingPlacementInfo(res_name.clone()))?
                .clone();
            final_placement_map.insert(residue_type, placement_data);

            let mut processed_rotamers = Vec::with_capacity(raw_rotamer_list.len());

            for raw_rotamer_data in raw_rotamer_list {
                processed_rotamers.push(Self::process_raw_rotamer(
                    &raw_rotamer_data,
                    &parameterizer,
                    &res_name,
                    rotamer_toml_path,
                )?);
            }

            final_rotamers_map.insert(residue_type, processed_rotamers);
        }

        Ok(Self {
            rotamers: final_rotamers_map,
            placement_info: final_placement_map,
        })
    }

    fn process_raw_rotamer(
        raw_rotamer_data: &RotamerData,
        parameterizer: &Parameterizer,
        res_name: &str,
        path_for_error: &Path,
    ) -> Result<Rotamer, LibraryLoadError> {
        let mut atoms = Vec::with_capacity(raw_rotamer_data.atoms.len());
        let mut serial_to_index_map = HashMap::with_capacity(raw_rotamer_data.atoms.len());
        let placeholder_residue_id = ResidueId::default();

        // 1. Process atoms and build serial-to-index map
        for (index, atom_data) in raw_rotamer_data.atoms.iter().enumerate() {
            if serial_to_index_map
                .insert(atom_data.serial, index)
                .is_some()
            {
                return Err(LibraryLoadError::DuplicateAtomSerial {
                    residue_type: res_name.to_string(),
                    serial: atom_data.serial,
                });
            }

            let mut atom = Atom::new(
                &atom_data.atom_name,
                placeholder_residue_id,
                Point3::from(atom_data.position),
            );
            atom.partial_charge = atom_data.partial_charge;
            atom.force_field_type = atom_data.force_field_type.clone();

            parameterizer
                .parameterize_atom(&mut atom, res_name)
                .map_err(|e| LibraryLoadError::Parameterization {
                    path: path_for_error.to_string_lossy().to_string(),
                    residue_type: res_name.to_string(),
                    source: e,
                })?;

            atoms.push(atom);
        }

        // 2. Process bonds using the map
        let mut bonds = Vec::with_capacity(raw_rotamer_data.bonds.len());
        for bond_serials in &raw_rotamer_data.bonds {
            let index1 = *serial_to_index_map.get(&bond_serials[0]).ok_or_else(|| {
                LibraryLoadError::InvalidBondSerial {
                    residue_type: res_name.to_string(),
                    serial: bond_serials[0],
                }
            })?;
            let index2 = *serial_to_index_map.get(&bond_serials[1]).ok_or_else(|| {
                LibraryLoadError::InvalidBondSerial {
                    residue_type: res_name.to_string(),
                    serial: bond_serials[1],
                }
            })?;
            bonds.push((index1, index2));
        }

        Ok(Rotamer { atoms, bonds })
    }

    pub fn get_rotamers_for(&self, residue_type: ResidueType) -> Option<&Vec<Rotamer>> {
        self.rotamers.get(&residue_type)
    }

    pub fn get_placement_info_for(&self, residue_type: ResidueType) -> Option<&PlacementInfo> {
        self.placement_info.get(&residue_type)
    }

    pub fn include_system_conformations(
        &mut self,
        system: &MolecularSystem,
        active_residues: &HashSet<ResidueId>,
    ) {
        for &residue_id in active_residues {
            if let Some(rotamer) = self.extract_rotamer_from_system(system, residue_id) {
                let residue = system.residue(residue_id).unwrap();
                let residue_type = residue.residue_type.unwrap();
                self.rotamers.entry(residue_type).or_default().push(rotamer);
            }
        }
    }

    fn extract_rotamer_from_system(
        &self,
        system: &MolecularSystem,
        residue_id: ResidueId,
    ) -> Option<Rotamer> {
        let residue = system.residue(residue_id)?;
        let residue_type = residue.residue_type?;
        let placement_info = self.placement_info.get(&residue_type)?;

        // --- 1. Create an "atom pool" from the residue in the system ---
        let mut atom_pool: HashMap<String, Vec<AtomId>> = HashMap::new();
        for atom_id in residue.atoms() {
            if let Some(atom) = system.atom(*atom_id) {
                atom_pool
                    .entry(atom.name.clone())
                    .or_default()
                    .push(*atom_id);
            }
        }

        // --- 2. Sequentially assign roles (anchor/sidechain) and build the new rotamer's atom list ---
        let mut new_rotamer_atoms = Vec::new();
        let mut old_id_to_new_index = HashMap::new();

        for name in &placement_info.anchor_atoms {
            if let Some(ids) = atom_pool.get_mut(name) {
                if !ids.is_empty() {
                    let atom_id = ids.remove(0);
                    if let Some(atom) = system.atom(atom_id) {
                        old_id_to_new_index.insert(atom_id, new_rotamer_atoms.len());
                        new_rotamer_atoms.push(atom.clone());
                    }
                } else {
                    return None;
                }
            }
        }

        for name in &placement_info.sidechain_atoms {
            if let Some(ids) = atom_pool.get_mut(name) {
                if let Some(atom_id) = ids.pop() {
                    if let Some(atom) = system.atom(atom_id) {
                        old_id_to_new_index.insert(atom_id, new_rotamer_atoms.len());
                        new_rotamer_atoms.push(atom.clone());
                    }
                } else {
                    return None;
                }
            }
        }

        // --- 3. Rebuild the topology for the extracted atoms ---
        let mut new_rotamer_bonds = Vec::new();
        for (&old_id_a, &new_idx_a) in &old_id_to_new_index {
            if let Some(neighbors) = system.get_bonded_neighbors(old_id_a) {
                for &old_id_b in neighbors {
                    if let Some(&new_idx_b) = old_id_to_new_index.get(&old_id_b) {
                        if new_idx_a < new_idx_b {
                            new_rotamer_bonds.push((new_idx_a, new_idx_b));
                        }
                    }
                }
            }
        }

        Some(Rotamer {
            atoms: new_rotamer_atoms,
            bonds: new_rotamer_bonds,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::params::{
        DeltaParam, Forcefield, GlobalParams, NonBondedParams, VdwParam,
    };
    use crate::core::models::atom::Atom;
    use crate::core::models::chain::ChainType;
    use crate::core::models::residue::ResidueType;
    use crate::core::models::system::MolecularSystem;
    use crate::core::models::topology::BondOrder;
    use crate::core::rotamers::placement::PlacementInfo;
    use nalgebra::Point3;
    use std::collections::{HashMap, HashSet};
    use std::fs;
    use tempfile::tempdir;

    fn dummy_forcefield() -> Forcefield {
        let mut deltas = HashMap::new();
        deltas.insert(
            ("ALA".to_string(), "CA".to_string()),
            DeltaParam {
                residue_type: "ALA".to_string(),
                atom_name: "CA".to_string(),
                mu: 1.0,
                sigma: 0.0,
            },
        );
        let mut vdw = HashMap::new();
        vdw.insert(
            "C_SP3".to_string(),
            VdwParam::LennardJones {
                radius: 2.0,
                well_depth: 0.1,
            },
        );
        let globals = GlobalParams {
            dielectric_constant: 1.0,
            potential_function: "lennard-jones-12-6".to_string(),
        };
        let non_bonded = NonBondedParams {
            globals,
            vdw,
            hbond: HashMap::new(),
        };
        Forcefield { deltas, non_bonded }
    }

    fn write_file(path: &std::path::Path, content: &str) {
        fs::write(path, content).expect("Failed to write temporary file for test setup");
    }

    #[test]
    fn load_rotamer_library_with_topology_success() {
        let dir = tempdir().unwrap();
        let rotamer_path = dir.path().join("rotamer.toml");
        let placement_path = dir.path().join("placement.toml");

        write_file(
            &rotamer_path,
            r#"
[[ALA]]
atoms = [
    { serial = 1, atom_name = "N", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_SP3" },
    { serial = 2, atom_name = "CA", position = [1.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_SP3" },
]
bonds = [ [1, 2] ]
"#,
        );
        write_file(
            &placement_path,
            r#"ALA = { sidechain_atoms = [], anchor_atoms = ["N", "CA"], connection_points = [], exact_match_atoms = [] }"#,
        );
        let ff = dummy_forcefield();
        let lib = RotamerLibrary::load(&rotamer_path, &placement_path, &ff, 0.0).unwrap();
        let rots = lib.get_rotamers_for(ResidueType::Alanine).unwrap();
        assert_eq!(rots.len(), 1);
        let rotamer1 = &rots[0];
        assert_eq!(rotamer1.atoms.len(), 2);
    }

    #[test]
    fn load_fails_on_duplicate_atom_serial() {
        let dir = tempdir().unwrap();
        let rotamer_path = dir.path().join("rotamer.toml");
        write_file(
            &rotamer_path,
            r#"
[[ALA]]
atoms = [
    { serial = 1, atom_name = "N", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_SP3" },
    { serial = 1, atom_name = "CA", position = [1.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_SP3" },
]
bonds = []
"#,
        );
        let placement_path = dir.path().join("placement.toml");
        write_file(
            &placement_path,
            r#"ALA = { sidechain_atoms = [], anchor_atoms = [], connection_points = [], exact_match_atoms = [] }"#,
        );
        let ff = dummy_forcefield();
        let result = RotamerLibrary::load(&rotamer_path, &placement_path, &ff, 0.0);
        assert!(matches!(
            result,
            Err(LibraryLoadError::DuplicateAtomSerial { .. })
        ));
    }

    #[test]
    fn load_fails_on_invalid_bond_serial() {
        let dir = tempdir().unwrap();
        let rotamer_path = dir.path().join("rotamer.toml");
        write_file(
            &rotamer_path,
            r#"
[[ALA]]
atoms = [ { serial = 1, atom_name = "N", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_SP3" } ]
bonds = [ [1, 99] ]
"#,
        );
        let placement_path = dir.path().join("placement.toml");
        write_file(
            &placement_path,
            r#"ALA = { sidechain_atoms = [], anchor_atoms = [], connection_points = [], exact_match_atoms = [] }"#,
        );
        let ff = dummy_forcefield();
        let result = RotamerLibrary::load(&rotamer_path, &placement_path, &ff, 0.0);
        assert!(matches!(
            result,
            Err(LibraryLoadError::InvalidBondSerial { .. })
        ));
    }

    fn create_library_with_placement_info() -> RotamerLibrary {
        let mut lib = RotamerLibrary::default();
        lib.placement_info.insert(
            ResidueType::Alanine,
            PlacementInfo {
                anchor_atoms: vec!["N".to_string(), "CA".to_string(), "C".to_string()],
                sidechain_atoms: vec![
                    "CB".to_string(),
                    "HCB".to_string(),
                    "HCB".to_string(),
                    "HCB".to_string(),
                ],
                connection_points: vec![],
                exact_match_atoms: vec![],
            },
        );
        lib.placement_info.insert(
            ResidueType::Glycine,
            PlacementInfo {
                anchor_atoms: vec![
                    "N".to_string(),
                    "CA".to_string(),
                    "C".to_string(),
                    "HCA".to_string(),
                ],
                sidechain_atoms: vec!["HCA".to_string()],
                connection_points: vec![],
                exact_match_atoms: vec![],
            },
        );
        lib
    }

    #[test]
    fn extract_rotamer_ala_with_multiple_hcb_correctly() {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();

        let n_id = system
            .add_atom_to_residue(res_id, Atom::new("N", res_id, Point3::origin()))
            .unwrap();
        let ca_id = system
            .add_atom_to_residue(res_id, Atom::new("CA", res_id, Point3::origin()))
            .unwrap();
        let c_id = system
            .add_atom_to_residue(res_id, Atom::new("C", res_id, Point3::origin()))
            .unwrap();
        let cb_id = system
            .add_atom_to_residue(res_id, Atom::new("CB", res_id, Point3::origin()))
            .unwrap();
        let hcb1_id = system
            .add_atom_to_residue(res_id, Atom::new("HCB", res_id, Point3::origin()))
            .unwrap();
        let hcb2_id = system
            .add_atom_to_residue(res_id, Atom::new("HCB", res_id, Point3::origin()))
            .unwrap();
        let hcb3_id = system
            .add_atom_to_residue(res_id, Atom::new("HCB", res_id, Point3::origin()))
            .unwrap();
        system.add_bond(ca_id, cb_id, BondOrder::Single).unwrap();
        system.add_bond(cb_id, hcb1_id, BondOrder::Single).unwrap();

        let mut lib = create_library_with_placement_info();
        let active_residues = [res_id].iter().cloned().collect();

        lib.include_system_conformations(&system, &active_residues);

        let rots = lib.get_rotamers_for(ResidueType::Alanine).unwrap();
        assert_eq!(rots.len(), 1);
        let rotamer = &rots[0];

        assert_eq!(rotamer.atoms.len(), 7);

        let atom_names: Vec<_> = rotamer.atoms.iter().map(|a| a.name.as_str()).collect();
        assert_eq!(&atom_names[0..3], &["N", "CA", "C"]);
        assert_eq!(&atom_names[3..], &["CB", "HCB", "HCB", "HCB"]);

        assert_eq!(rotamer.bonds.len(), 2, "Should have CA-CB and CB-HCB bonds");

        let name_to_idx: HashMap<_, _> = rotamer
            .atoms
            .iter()
            .enumerate()
            .map(|(i, a)| (a.name.clone(), i))
            .collect();
        let bonds_as_set: HashSet<(usize, usize)> = rotamer
            .bonds
            .iter()
            .map(|&(a, b)| (a.min(b), a.max(b)))
            .collect();
        let ca_idx = name_to_idx["CA"];
        let cb_idx = name_to_idx["CB"];

        assert!(bonds_as_set.contains(&(ca_idx, cb_idx)));
    }

    #[test]
    fn extract_rotamer_gly_handles_dual_role_hca_correctly() {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_id = system
            .add_residue(chain_id, 1, "GLY", Some(ResidueType::Glycine))
            .unwrap();

        let n_id = system
            .add_atom_to_residue(res_id, Atom::new("N", res_id, Point3::new(0.0, 1.0, 0.0)))
            .unwrap();
        let ca_id = system
            .add_atom_to_residue(res_id, Atom::new("CA", res_id, Point3::new(0.0, 0.0, 0.0)))
            .unwrap();
        let c_id = system
            .add_atom_to_residue(res_id, Atom::new("C", res_id, Point3::new(1.0, 0.0, 0.0)))
            .unwrap();
        let hca1_id = system
            .add_atom_to_residue(
                res_id,
                Atom::new("HCA", res_id, Point3::new(-0.5, -0.5, 0.0)),
            )
            .unwrap();
        let hca2_id = system
            .add_atom_to_residue(
                res_id,
                Atom::new("HCA", res_id, Point3::new(-0.5, -0.5, 1.0)),
            )
            .unwrap();

        system.add_bond(n_id, ca_id, BondOrder::Single).unwrap();
        system.add_bond(ca_id, c_id, BondOrder::Single).unwrap();
        system.add_bond(ca_id, hca1_id, BondOrder::Single).unwrap();
        system.add_bond(ca_id, hca2_id, BondOrder::Single).unwrap();

        let mut lib = create_library_with_placement_info();
        let active_residues = [res_id].iter().cloned().collect();

        lib.include_system_conformations(&system, &active_residues);

        let rots = lib.get_rotamers_for(ResidueType::Glycine).unwrap();
        assert_eq!(rots.len(), 1);
        let rotamer = &rots[0];

        assert_eq!(rotamer.atoms.len(), 5);

        let atom_ids_in_rotamer: Vec<_> = rotamer.atoms.iter().map(|a| a.residue_id).collect();
        let atom_positions_in_rotamer: Vec<_> = rotamer.atoms.iter().map(|a| a.position).collect();

        assert_eq!(atom_positions_in_rotamer[3], Point3::new(-0.5, -0.5, 0.0));
        assert_eq!(atom_positions_in_rotamer[4], Point3::new(-0.5, -0.5, 1.0));

        let ca_idx = 1;
        let hca_anchor_idx = 3;
        let hca_sidechain_idx = 4;
        let bonds_as_set: HashSet<(usize, usize)> = rotamer
            .bonds
            .iter()
            .map(|&(a, b)| (a.min(b), a.max(b)))
            .collect();

        assert!(
            bonds_as_set.contains(&(ca_idx, hca_anchor_idx)),
            "Bond between CA and anchor HCA should exist"
        );
        assert!(
            bonds_as_set.contains(&(ca_idx, hca_sidechain_idx)),
            "Bond between CA and sidechain HCA should exist"
        );
    }

    #[test]
    fn include_system_conformations_is_idempotent_in_adding() {
        let mut lib = create_library_with_placement_info();

        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        system
            .add_atom_to_residue(res_id, Atom::new("N", res_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(res_id, Atom::new("CA", res_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(res_id, Atom::new("C", res_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(res_id, Atom::new("CB", res_id, Point3::origin()))
            .unwrap();

        let active_residues = [res_id].iter().cloned().collect();

        lib.include_system_conformations(&system, &active_residues);
        assert_eq!(
            lib.get_rotamers_for(ResidueType::Alanine).unwrap().len(),
            1,
            "Should add the rotamer the first time"
        );

        lib.include_system_conformations(&system, &active_residues);
        assert_eq!(
            lib.get_rotamers_for(ResidueType::Alanine).unwrap().len(),
            2,
            "Should add the rotamer again as de-duplication is not implemented"
        );
    }

    #[test]
    fn include_system_conformations_skips_residue_if_placement_info_is_missing() {
        let mut lib = RotamerLibrary::default();
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();

        let active_residues = [res_id].iter().cloned().collect();
        lib.include_system_conformations(&system, &active_residues);

        assert!(
            lib.get_rotamers_for(ResidueType::Alanine).is_none(),
            "Should not add rotamer if placement info is missing"
        );
    }
}
