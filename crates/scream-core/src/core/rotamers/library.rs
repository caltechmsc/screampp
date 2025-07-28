use super::placement::PlacementInfo;
use super::rotamer::{Rotamer, RotamerData};
use crate::core::forcefield::parameterization::{ParameterizationError, Parameterizer};
use crate::core::forcefield::params::Forcefield;
use crate::core::models::atom::Atom;
use crate::core::models::ids::ResidueId;
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
            let residue = match system.residue(residue_id) {
                Some(r) => r,
                None => continue,
            };

            let residue_type = match residue.residue_type {
                Some(rt) => rt,
                None => continue,
            };

            let placement_info = match self.placement_info.get(&residue_type) {
                Some(info) => info,
                None => continue,
            };

            // --- 1. Extract atoms and build a map for topology extraction ---
            let all_rotamer_atom_names: HashSet<_> = placement_info
                .anchor_atoms
                .iter()
                .chain(&placement_info.sidechain_atoms)
                .cloned()
                .collect();

            let mut extracted_atoms = Vec::new();
            let mut global_id_to_local_index = HashMap::new();

            for atom_name in &all_rotamer_atom_names {
                if let Some(atom_id) = residue.get_atom_id_by_name(atom_name) {
                    if let Some(atom) = system.atom(atom_id) {
                        let local_index = extracted_atoms.len();
                        extracted_atoms.push(atom.clone());
                        global_id_to_local_index.insert(atom_id, local_index);
                    }
                }
            }

            // --- 2. Extract topology (bonds) for the extracted atoms ---
            let mut extracted_bonds = Vec::new();
            for (atom_id_a, &local_index_a) in &global_id_to_local_index {
                if let Some(neighbors) = system.get_bonded_neighbors(*atom_id_a) {
                    for &atom_id_b in neighbors {
                        if let Some(&local_index_b) = global_id_to_local_index.get(&atom_id_b) {
                            if local_index_a < local_index_b {
                                extracted_bonds.push((local_index_a, local_index_b));
                            }
                        }
                    }
                }
            }

            let extracted_rotamer = Rotamer {
                atoms: extracted_atoms,
                bonds: extracted_bonds,
            };

            // --- 3. Add the extracted rotamer to the library ---
            self.rotamers
                .entry(residue_type)
                .or_default()
                .push(extracted_rotamer);
        }
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

    fn ala_placement_info() -> PlacementInfo {
        PlacementInfo {
            sidechain_atoms: vec!["CB".to_string()],
            anchor_atoms: vec!["N".to_string(), "CA".to_string(), "C".to_string()],
            connection_points: vec![],
            exact_match_atoms: vec![],
        }
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

[[ALA]]
atoms = [
    { serial = 10, atom_name = "N", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_SP3" },
    { serial = 20, atom_name = "CA", position = [1.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_SP3" },
]
bonds = [ [10, 20] ]
"#,
        );

        write_file(
            &placement_path,
            r#"ALA = { sidechain_atoms = ["CB"], anchor_atoms = ["N", "CA", "C"], connection_points = [], exact_match_atoms = [] }"#,
        );
        let ff = dummy_forcefield();
        let lib = RotamerLibrary::load(&rotamer_path, &placement_path, &ff, 0.0).unwrap();

        let rots = lib.get_rotamers_for(ResidueType::Alanine).unwrap();
        assert_eq!(rots.len(), 2, "Should load two rotamers for ALA");

        let rotamer1 = &rots[0];
        assert_eq!(rotamer1.atoms.len(), 2);
        assert_eq!(rotamer1.bonds.len(), 1);
        assert_eq!(rotamer1.bonds[0], (0, 1));

        let rotamer2 = &rots[1];
        assert_eq!(rotamer2.atoms.len(), 2);
        assert_eq!(rotamer2.bonds.len(), 1);
        assert_eq!(rotamer2.bonds[0], (0, 1));
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

        assert!(
            matches!(result, Err(LibraryLoadError::DuplicateAtomSerial { residue_type, serial }) if residue_type == "ALA" && serial == 1)
        );
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

        assert!(
            matches!(result, Err(LibraryLoadError::InvalidBondSerial { residue_type, serial }) if residue_type == "ALA" && serial == 99)
        );
    }

    fn create_test_system_for_extraction() -> (MolecularSystem, ResidueId) {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();

        let n_atom = Atom::new("N", res_id, Point3::origin());
        let ca_atom = Atom::new("CA", res_id, Point3::new(1.0, 0.0, 0.0));
        let c_atom = Atom::new("C", res_id, Point3::new(1.0, 1.0, 0.0));
        let cb_atom = Atom::new("CB", res_id, Point3::new(2.0, -1.0, 0.0));

        let n_id = system.add_atom_to_residue(res_id, n_atom).unwrap();
        let ca_id = system.add_atom_to_residue(res_id, ca_atom).unwrap();
        let c_id = system.add_atom_to_residue(res_id, c_atom).unwrap();
        let cb_id = system.add_atom_to_residue(res_id, cb_atom).unwrap();

        system.add_bond(n_id, ca_id, BondOrder::Single).unwrap();
        system.add_bond(ca_id, c_id, BondOrder::Single).unwrap();
        system.add_bond(ca_id, cb_id, BondOrder::Single).unwrap();

        (system, res_id)
    }

    #[test]
    fn include_system_conformations_extracts_atoms_and_topology() {
        let mut lib = RotamerLibrary::default();
        let placement_info = ala_placement_info();
        lib.placement_info
            .insert(ResidueType::Alanine, placement_info);
        lib.rotamers.insert(ResidueType::Alanine, vec![]);

        let (system, res_id) = create_test_system_for_extraction();
        let mut active = HashSet::new();
        active.insert(res_id);

        lib.include_system_conformations(&system, &active);

        let rots = lib.get_rotamers_for(ResidueType::Alanine).unwrap();
        assert_eq!(rots.len(), 1);
        let rotamer = &rots[0];

        assert_eq!(rotamer.atoms.len(), 4);
        assert!(rotamer.atoms.iter().any(|a| a.name == "N"));
        assert!(rotamer.atoms.iter().any(|a| a.name == "CA"));
        assert!(rotamer.atoms.iter().any(|a| a.name == "C"));
        assert!(rotamer.atoms.iter().any(|a| a.name == "CB"));

        assert_eq!(rotamer.bonds.len(), 3);
        let name_to_idx: HashMap<_, _> = rotamer
            .atoms
            .iter()
            .enumerate()
            .map(|(idx, atom)| (atom.name.as_str(), idx))
            .collect();
        let expected_bonds: HashSet<(usize, usize)> = [
            (name_to_idx["N"], name_to_idx["CA"]),
            (name_to_idx["CA"], name_to_idx["C"]),
            (name_to_idx["CA"], name_to_idx["CB"]),
        ]
        .iter()
        .map(|&(a, b)| if a < b { (a, b) } else { (b, a) })
        .collect();

        let actual_bonds: HashSet<(usize, usize)> = rotamer
            .bonds
            .iter()
            .map(|&(a, b)| if a < b { (a, b) } else { (b, a) })
            .collect();

        assert_eq!(actual_bonds, expected_bonds);
    }

    #[test]
    fn include_system_conformations_always_adds_extracted_rotamer() {
        let mut lib = RotamerLibrary::default();
        let placement_info = ala_placement_info();
        lib.placement_info
            .insert(ResidueType::Alanine, placement_info.clone());
        lib.rotamers.insert(ResidueType::Alanine, vec![]);

        let (system, res_id) = create_test_system_for_extraction();
        let mut active = HashSet::new();
        active.insert(res_id);

        lib.include_system_conformations(&system, &active);
        let rots = lib.get_rotamers_for(ResidueType::Alanine).unwrap();
        assert_eq!(rots.len(), 1, "Should add the rotamer the first time");

        lib.include_system_conformations(&system, &active);
        let rots_after_second_call = lib.get_rotamers_for(ResidueType::Alanine).unwrap();
        assert_eq!(
            rots_after_second_call.len(),
            2,
            "Should add the rotamer again as de-duplication is disabled"
        );
    }
}
