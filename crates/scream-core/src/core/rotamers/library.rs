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
        rmsd_threshold: f64,
    ) {
        let mut new_rotamers_to_add: HashMap<ResidueType, Rotamer> = HashMap::new();

        for &residue_id in active_residues {
            let residue = match system.residue(residue_id) {
                Some(r) => r,
                None => continue,
            };

            let residue_type = match residue.residue_type {
                Some(rt) => rt,
                None => continue,
            };

            let (existing_rotamers, placement_info) = match (
                self.rotamers.get(&residue_type),
                self.placement_info.get(&residue_type),
            ) {
                (Some(rots), Some(info)) => (rots, info),
                _ => continue,
            };

            // 1. Extract atoms and build a global_id -> local_index map for topology extraction.
            let mut extracted_atoms = Vec::new();
            let mut global_id_to_local_index = HashMap::new();
            let mut sidechain_atom_ids = HashSet::new();

            let all_rotamer_atoms: Vec<_> = placement_info
                .anchor_atoms
                .iter()
                .chain(&placement_info.sidechain_atoms)
                .collect();

            for atom_name in all_rotamer_atoms {
                if let Some(atom_id) = residue.get_atom_id_by_name(atom_name) {
                    if let Some(atom) = system.atom(atom_id) {
                        let local_index = extracted_atoms.len();
                        extracted_atoms.push(atom.clone());
                        global_id_to_local_index.insert(atom_id, local_index);
                        if placement_info.sidechain_atoms.contains(atom_name) {
                            sidechain_atom_ids.insert(atom_id);
                        }
                    }
                }
            }

            // 2. Extract topology (bonds) for the extracted atoms.
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

            // 3. Check for duplicates based on RMSD of side-chain atoms only.
            let is_duplicate = existing_rotamers.iter().any(|existing| {
                let name_to_pos1: HashMap<_, _> = extracted_rotamer
                    .atoms
                    .iter()
                    .filter(|a| {
                        sidechain_atom_ids.contains(
                            &system
                                .residue(a.residue_id)
                                .unwrap()
                                .get_atom_id_by_name(&a.name)
                                .unwrap(),
                        )
                    })
                    .map(|a| (&a.name, a.position))
                    .collect();

                let name_to_pos2: HashMap<_, _> = existing
                    .atoms
                    .iter()
                    .filter(|a| placement_info.sidechain_atoms.contains(&a.name))
                    .map(|a| (&a.name, a.position))
                    .collect();

                let mut coords1 = Vec::new();
                let mut coords2 = Vec::new();

                for (name, pos1) in &name_to_pos1 {
                    if let Some(pos2) = name_to_pos2.get(name) {
                        coords1.push(*pos1);
                        coords2.push(*pos2);
                    }
                }

                if let Some(rmsd) = calculate_rmsd(&coords1, &coords2) {
                    rmsd < rmsd_threshold
                } else {
                    false
                }
            });

            if !is_duplicate {
                new_rotamers_to_add
                    .entry(residue_type)
                    .or_insert(extracted_rotamer);
            }
        }

        for (residue_type, new_rotamer) in new_rotamers_to_add {
            if let Some(rotamers) = self.rotamers.get_mut(&residue_type) {
                rotamers.push(new_rotamer);
            }
        }
    }
}

fn calculate_rmsd(coords1: &[Point3<f64>], coords2: &[Point3<f64>]) -> Option<f64> {
    if coords1.len() != coords2.len() || coords1.is_empty() {
        return None;
    }
    let n = coords1.len() as f64;
    let squared_dist_sum: f64 = coords1
        .iter()
        .zip(coords2.iter())
        .map(|(p1, p2)| (p1 - p2).norm_squared())
        .sum();
    Some((squared_dist_sum / n).sqrt())
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
            VdwParam::Buckingham {
                radius: 2.0,
                well_depth: 0.1,
                scale: 12.0,
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

    fn dummy_placement_info() -> PlacementInfo {
        PlacementInfo {
            sidechain_atoms: vec!["CA".to_string()],
            anchor_atoms: vec![],
            connection_points: vec![],
            exact_match_atoms: vec![],
        }
    }

    fn write_file(path: &std::path::Path, content: &str) {
        fs::write(path, content).expect("Failed to write temporary file for test setup");
    }

    #[test]
    fn load_rotamer_library_success() {
        let dir = tempdir().unwrap();
        let rotamer_path = dir.path().join("rotamer.toml");
        let placement_path = dir.path().join("placement.toml");
        write_file(
            &rotamer_path,
            r#"
ALA = [
    { atoms = [
        { serial = 1, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_SP3" }
    ] }
]
"#,
        );
        write_file(
            &placement_path,
            r#"
ALA = { sidechain_atoms = ["CA"], anchor_atoms = [], connection_points = [], exact_match_atoms = [] }
"#,
        );
        let ff = dummy_forcefield();
        let lib = RotamerLibrary::load(&rotamer_path, &placement_path, &ff, 0.0);
        assert!(lib.is_ok());
        let lib = lib.unwrap();
        assert_eq!(lib.rotamers.len(), 1);
        assert_eq!(lib.placement_info.len(), 1);
        let rots = lib.get_rotamers_for(ResidueType::Alanine).unwrap();
        assert_eq!(rots.len(), 1);
        assert_eq!(rots[0].atoms.len(), 1);
        let ca = &rots[0].atoms[0];
        assert_eq!(ca.name, "CA");
        match ca.vdw_param {
            crate::core::models::atom::CachedVdwParam::Buckingham {
                radius, well_depth, ..
            } => {
                assert_eq!(radius, 2.0);
                assert_eq!(well_depth, 0.1);
            }
            _ => panic!("Expected Buckingham variant for CA atom"),
        }
        let placement = lib.get_placement_info_for(ResidueType::Alanine).unwrap();
        assert_eq!(placement.sidechain_atoms, vec!["CA"]);
    }

    #[test]
    fn load_rotamer_library_missing_placement_info() {
        let dir = tempdir().unwrap();
        let rotamer_path = dir.path().join("rotamer.toml");
        let placement_path = dir.path().join("placement.toml");
        write_file(
            &rotamer_path,
            r#"
ALA = [
    { atoms = [
        { serial = 1, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_SP3" }
    ] }
]
"#,
        );
        write_file(&placement_path, "");
        let ff = dummy_forcefield();
        let lib = RotamerLibrary::load(&rotamer_path, &placement_path, &ff, 0.0);
        assert!(matches!(
            lib,
            Err(LibraryLoadError::Toml { .. }) | Err(LibraryLoadError::MissingPlacementInfo(_))
        ));
    }

    #[test]
    fn load_rotamer_library_returns_io_error_for_nonexistent_file() {
        let ff = dummy_forcefield();
        let rotamer_path = std::path::Path::new("nonexistent_rotamer_file.toml");
        let placement_path = std::path::Path::new("nonexistent_placement_file.toml");
        let lib = RotamerLibrary::load(rotamer_path, placement_path, &ff, 0.0);
        assert!(matches!(lib, Err(LibraryLoadError::Io { .. })));
    }

    #[test]
    fn load_rotamer_library_returns_toml_error_for_malformed_file() {
        let dir = tempdir().unwrap();
        let rotamer_path = dir.path().join("malformed.toml");
        let placement_path = dir.path().join("placement.toml");
        write_file(&rotamer_path, "this is not toml");
        write_file(
            &placement_path,
            r#"ALA = { sidechain_atoms = ["CA"], anchor_atoms = [], connection_points = [], exact_match_atoms = [] }"#,
        );
        let ff = dummy_forcefield();
        let lib = RotamerLibrary::load(&rotamer_path, &placement_path, &ff, 0.0);
        assert!(matches!(lib, Err(LibraryLoadError::Toml { .. })));
    }

    #[test]
    fn load_rotamer_library_returns_error_for_unknown_residue_type() {
        let dir = tempdir().unwrap();
        let rotamer_path = dir.path().join("rotamer.toml");
        let placement_path = dir.path().join("placement.toml");
        write_file(
            &rotamer_path,
            r#"
UNK = [ { atoms = [] } ]
"#,
        );
        write_file(
            &placement_path,
            r#"ALA = { sidechain_atoms = ["CA"], anchor_atoms = [], connection_points = [], exact_match_atoms = [] }"#,
        );
        let ff = dummy_forcefield();
        let lib = RotamerLibrary::load(&rotamer_path, &placement_path, &ff, 0.0);
        assert!(matches!(lib, Err(LibraryLoadError::UnknownResidueType(name)) if name == "UNK"));
    }

    #[test]
    fn load_rotamer_library_returns_error_for_parameterization_failure() {
        let dir = tempdir().unwrap();
        let rotamer_path = dir.path().join("rotamer.toml");
        let placement_path = dir.path().join("placement.toml");
        write_file(
            &rotamer_path,
            r#"
ALA = [ { atoms = [ { serial = 1, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "UnknownType" } ] } ]
"#,
        );
        write_file(
            &placement_path,
            r#"ALA = { sidechain_atoms = ["CA"], anchor_atoms = [], connection_points = [], exact_match_atoms = [] }"#,
        );
        let ff = dummy_forcefield();
        let lib = RotamerLibrary::load(&rotamer_path, &placement_path, &ff, 0.0);
        assert!(matches!(
            lib,
            Err(LibraryLoadError::Parameterization { .. })
        ));
    }

    #[test]
    fn load_rotamer_library_empty_file() {
        let dir = tempdir().unwrap();
        let rotamer_path = dir.path().join("empty.toml");
        let placement_path = dir.path().join("placement.toml");
        write_file(&rotamer_path, "");
        write_file(&placement_path, "");
        let ff = dummy_forcefield();
        let lib = RotamerLibrary::load(&rotamer_path, &placement_path, &ff, 0.0);
        if let Ok(library) = lib {
            assert!(library.rotamers.is_empty());
            assert!(library.placement_info.is_empty());
        }
    }

    #[test]
    fn load_rotamer_library_handles_residue_with_empty_rotamer_list() {
        let dir = tempdir().unwrap();
        let rotamer_path = dir.path().join("rotamer.toml");
        let placement_path = dir.path().join("placement.toml");
        write_file(
            &rotamer_path,
            r#"
ALA = [ { atoms = [] } ]
"#,
        );
        write_file(
            &placement_path,
            r#"ALA = { sidechain_atoms = ["CA"], anchor_atoms = [], connection_points = [], exact_match_atoms = [] }"#,
        );
        let ff = dummy_forcefield();
        let lib = RotamerLibrary::load(&rotamer_path, &placement_path, &ff, 0.0).unwrap();
        assert_eq!(lib.rotamers.len(), 1);
        let rots = lib.get_rotamers_for(ResidueType::Alanine).unwrap();
        assert_eq!(rots.len(), 1);
        assert!(rots[0].atoms.is_empty());
    }

    #[test]
    fn include_system_conformations_adds_unique_original_rotamer() {
        let mut lib = RotamerLibrary::default();
        let placement_info = dummy_placement_info();
        lib.placement_info
            .insert(ResidueType::Alanine, placement_info.clone());
        lib.rotamers.insert(ResidueType::Alanine, vec![]);
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let mut atom = Atom::new("CA", res_id, Point3::new(5.0, 5.0, 5.0));
        atom.force_field_type = "C_SP3".to_string();
        system.add_atom_to_residue(res_id, atom).unwrap();
        let mut active = HashSet::new();
        active.insert(res_id);
        let original_count = lib.get_rotamers_for(ResidueType::Alanine).unwrap().len();
        lib.include_system_conformations(&system, &active, 0.1);
        let new_count = lib.get_rotamers_for(ResidueType::Alanine).unwrap().len();
        assert_eq!(new_count, original_count + 1);
    }

    #[test]
    fn include_system_conformations_skips_duplicate_rotamer() {
        let mut lib = RotamerLibrary::default();
        let placement_info = dummy_placement_info();
        lib.placement_info
            .insert(ResidueType::Alanine, placement_info.clone());
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let mut atom = Atom::new("CA", res_id, Point3::new(1.0, 2.0, 3.0));
        atom.force_field_type = "C_SP3".to_string();
        system.add_atom_to_residue(res_id, atom.clone()).unwrap();
        lib.rotamers.insert(
            ResidueType::Alanine,
            vec![Rotamer {
                atoms: vec![atom.clone()],
            }],
        );
        let mut active = HashSet::new();
        active.insert(res_id);
        let original_count = lib.get_rotamers_for(ResidueType::Alanine).unwrap().len();
        lib.include_system_conformations(&system, &active, 0.1);
        let new_count = lib.get_rotamers_for(ResidueType::Alanine).unwrap().len();
        assert_eq!(new_count, original_count);
    }
}
