use super::placement::PlacementInfo;
use super::rotamer::{Rotamer, RotamerData};
use crate::core::forcefield::parameterization::{ParameterizationError, Parameterizer};
use crate::core::forcefield::params::Forcefield;
use crate::core::models::atom::Atom;
use crate::core::models::ids::ResidueId;
use crate::core::models::residue::ResidueType;
use crate::core::models::system::MolecularSystem;
use crate::core::utils::geometry::calculate_rmsd;
use nalgebra::Point3;
use std::collections::HashMap;
use std::collections::HashSet;
use std::path::Path;
use std::str::FromStr;
use thiserror::Error;

type RawRotamerFile = HashMap<String, Vec<RotamerData>>;

#[derive(Debug, Default)]
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
        "Parameterization failed for residue '{res_type}' in rotamer from file '{path}': {source}"
    )]
    Parameterization {
        path: String,
        res_type: String,
        source: ParameterizationError,
    },
}

impl RotamerLibrary {
    pub fn load(
        rotamer_toml_path: &Path,
        placement_registry: &HashMap<String, PlacementInfo>,
        forcefield: &Forcefield,
        delta_s_factor: f64,
    ) -> Result<Self, LibraryLoadError> {
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

        let parameterizer = Parameterizer::new(forcefield.clone(), delta_s_factor);
        let mut final_rotamers_map = HashMap::new();
        let mut final_placement_map = HashMap::new();

        for (res_name, raw_rotamer_list) in raw_lib {
            let res_type = ResidueType::from_str(&res_name)
                .map_err(|_| LibraryLoadError::UnknownResidueType(res_name.clone()))?;

            let placement_data = placement_registry
                .get(&res_name)
                .ok_or_else(|| LibraryLoadError::MissingPlacementInfo(res_name.clone()))?
                .clone();
            final_placement_map.insert(res_type, placement_data);

            let mut parameterized_rotamers = Vec::with_capacity(raw_rotamer_list.len());
            let placeholder_residue_id = ResidueId::default();

            for raw_rotamer in raw_rotamer_list {
                let mut atoms = Vec::with_capacity(raw_rotamer.atoms.len());
                for atom_data in raw_rotamer.atoms {
                    let mut atom = Atom::new(
                        atom_data.serial,
                        &atom_data.atom_name,
                        placeholder_residue_id,
                        Point3::from(atom_data.position),
                    );
                    atom.partial_charge = atom_data.partial_charge;
                    atom.force_field_type = atom_data.force_field_type;

                    parameterizer
                        .parameterize_atom(&mut atom, &res_name)
                        .map_err(|e| LibraryLoadError::Parameterization {
                            path: rotamer_toml_path.to_string_lossy().to_string(),
                            res_type: res_name.clone(),
                            source: e,
                        })?;

                    atoms.push(atom);
                }
                parameterized_rotamers.push(Rotamer {
                    atoms,
                    empty_lattice_energy: None,
                });
            }

            final_rotamers_map.insert(res_type, parameterized_rotamers);
        }

        Ok(Self {
            rotamers: final_rotamers_map,
            placement_info: final_placement_map,
        })
    }

    pub fn get_rotamers_for(&self, res_type: ResidueType) -> Option<&Vec<Rotamer>> {
        self.rotamers.get(&res_type)
    }

    pub fn get_placement_info_for(&self, res_type: ResidueType) -> Option<&PlacementInfo> {
        self.placement_info.get(&res_type)
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

            let res_type = match residue.res_type {
                Some(rt) => rt,
                None => continue,
            };

            let (existing_rotamers, placement_info) = match (
                self.rotamers.get(&res_type),
                self.placement_info.get(&res_type),
            ) {
                (Some(rots), Some(info)) => (rots, info),
                _ => continue,
            };

            let mut extracted_atoms = Vec::new();
            for atom_name in &placement_info.sidechain_atoms {
                if let Some(atom_id) = residue.get_atom_id_by_name(atom_name) {
                    if let Some(atom) = system.atom(atom_id) {
                        extracted_atoms.push(atom.clone());
                    }
                }
            }

            if extracted_atoms.len() != placement_info.sidechain_atoms.len() {
                continue;
            }

            let extracted_rotamer = Rotamer {
                atoms: extracted_atoms,
                empty_lattice_energy: None,
            };

            let is_duplicate = existing_rotamers.iter().any(|existing| {
                let coords1: Vec<_> = extracted_rotamer.atoms.iter().map(|a| a.position).collect();
                let coords2: Vec<_> = existing.atoms.iter().map(|a| a.position).collect();

                if let Some(rmsd) = calculate_rmsd(&coords1, &coords2) {
                    rmsd < rmsd_threshold
                } else {
                    false
                }
            });

            if !is_duplicate {
                new_rotamers_to_add
                    .entry(res_type)
                    .or_insert(extracted_rotamer);
            }
        }

        for (res_type, new_rotamer) in new_rotamers_to_add {
            if let Some(rotamers) = self.rotamers.get_mut(&res_type) {
                rotamers.push(new_rotamer);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::params::{Forcefield, GlobalParams, NonBondedParams, VdwParam};
    use crate::core::models::atom::Atom;
    use crate::core::models::chain::ChainType;
    use crate::core::models::system::MolecularSystem;
    use crate::core::rotamers::placement::PlacementInfo;
    use nalgebra::Point3;
    use std::collections::HashSet;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_forcefield() -> Forcefield {
        let globals = GlobalParams {
            dielectric_constant: 1.0,
            potential_function: "mixed".to_string(),
        };

        let mut vdw = HashMap::new();
        vdw.insert(
            "CA".to_string(),
            VdwParam::LennardJones {
                radius: 1.9,
                well_depth: 0.1,
            },
        );
        vdw.insert(
            "CB".to_string(),
            VdwParam::LennardJones {
                radius: 1.9,
                well_depth: 0.1,
            },
        );
        vdw.insert(
            "O".to_string(),
            VdwParam::LennardJones {
                radius: 1.7,
                well_depth: 0.2,
            },
        );

        let non_bonded = NonBondedParams {
            globals,
            vdw,
            hbond: HashMap::new(),
        };

        Forcefield {
            non_bonded,
            deltas: HashMap::new(),
        }
    }

    fn create_test_placement_registry() -> HashMap<String, PlacementInfo> {
        let mut registry = HashMap::new();
        registry.insert(
            "ALA".to_string(),
            PlacementInfo {
                anchor_atoms: vec!["N".to_string(), "CA".to_string(), "C".to_string()],
                sidechain_atoms: vec!["CB".to_string()],
                exact_match_atoms: vec![
                    "N".to_string(),
                    "CA".to_string(),
                    "C".to_string(),
                    "O".to_string(),
                ],
                connection_points: vec!["C".to_string(), "N".to_string()],
            },
        );
        registry
    }

    fn create_temp_rotamer_file(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        write!(file, "{}", content).unwrap();
        file
    }

    fn create_test_system_with_ala(cb_pos: Point3<f64>) -> (MolecularSystem, ResidueId) {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();

        let mut cb_atom = Atom::new(10, "CB", res_id, cb_pos);
        cb_atom.force_field_type = "CB".to_string();
        system.add_atom_to_residue(res_id, cb_atom).unwrap();

        (system, res_id)
    }

    #[test]
    fn loads_rotamers_and_placement_info_successfully() {
        let rotamer_content = r#"
[[ALA]]
atoms = [
    { serial = 1, atom_name = "CA", partial_charge = 0.1, position = [0.0, 0.0, 0.0], force_field_type = "CA" },
    { serial = 2, atom_name = "CB", partial_charge = -0.1, position = [1.5, 0.0, 0.0], force_field_type = "CB" },
]
"#;
        let rotamer_file = create_temp_rotamer_file(rotamer_content);

        let ff = create_test_forcefield();
        let placement_registry = create_test_placement_registry();

        let library =
            RotamerLibrary::load(rotamer_file.path(), &placement_registry, &ff, 0.0).unwrap();

        assert_eq!(library.rotamers.len(), 1);
        assert_eq!(library.placement_info.len(), 1);

        let ala_rotamers = library.get_rotamers_for(ResidueType::Alanine).unwrap();
        assert_eq!(ala_rotamers.len(), 1);
        assert_eq!(ala_rotamers[0].atoms.len(), 2);

        let ca_atom = &ala_rotamers[0]
            .atoms
            .iter()
            .find(|a| a.name == "CA")
            .unwrap();
        assert_eq!(ca_atom.vdw_radius, 1.9);
        assert_eq!(ca_atom.vdw_well_depth, 0.1);

        let ala_placement = library
            .get_placement_info_for(ResidueType::Alanine)
            .unwrap();
        assert_eq!(ala_placement.sidechain_atoms, vec!["CB"]);
    }

    #[test]
    fn load_returns_io_error_for_nonexistent_file() {
        let path = Path::new("nonexistent_rotamer_file.toml");
        let ff = create_test_forcefield();
        let placement_registry = create_test_placement_registry();

        let result = RotamerLibrary::load(path, &placement_registry, &ff, 0.0);

        assert!(matches!(result, Err(LibraryLoadError::Io { .. })));
    }

    #[test]
    fn load_returns_toml_error_for_malformed_file() {
        let content = "this is not valid toml";
        let file = create_temp_rotamer_file(content);
        let ff = create_test_forcefield();
        let placement_registry = create_test_placement_registry();

        let result = RotamerLibrary::load(file.path(), &placement_registry, &ff, 0.0);

        assert!(matches!(result, Err(LibraryLoadError::Toml { .. })));
    }

    #[test]
    fn load_returns_error_for_unknown_residue_type() {
        let rotamer_content = r#"
[[UNK]]
atoms = []
"#;
        let rotamer_file = create_temp_rotamer_file(rotamer_content);
        let ff = create_test_forcefield();
        let placement_registry = create_test_placement_registry();

        let result = RotamerLibrary::load(rotamer_file.path(), &placement_registry, &ff, 0.0);

        assert!(matches!(
            result,
            Err(LibraryLoadError::UnknownResidueType(name)) if name == "UNK"
        ));
    }

    #[test]
    fn load_returns_error_for_missing_placement_info() {
        let rotamer_content = r#"
[[ALA]]
atoms = []
"#;
        let rotamer_file = create_temp_rotamer_file(rotamer_content);
        let ff = create_test_forcefield();
        let empty_placement_registry = HashMap::new();

        let result = RotamerLibrary::load(rotamer_file.path(), &empty_placement_registry, &ff, 0.0);

        assert!(matches!(
            result,
            Err(LibraryLoadError::MissingPlacementInfo(name)) if name == "ALA"
        ));
    }

    #[test]
    fn load_returns_error_for_parameterization_failure() {
        let rotamer_content = r#"
[[ALA]]
atoms = [
    { serial = 1, atom_name = "CA", partial_charge = 0.1, position = [0.0, 0.0, 0.0], force_field_type = "UnknownType" }
]
"#;
        let rotamer_file = create_temp_rotamer_file(rotamer_content);
        let ff = create_test_forcefield();
        let placement_registry = create_test_placement_registry();

        let result = RotamerLibrary::load(rotamer_file.path(), &placement_registry, &ff, 0.0);

        assert!(matches!(
            result,
            Err(LibraryLoadError::Parameterization { .. })
        ));
    }

    #[test]
    fn loads_empty_library_from_empty_file() {
        let rotamer_file = create_temp_rotamer_file("");
        let ff = create_test_forcefield();
        let placement_registry = create_test_placement_registry();

        let library =
            RotamerLibrary::load(rotamer_file.path(), &placement_registry, &ff, 0.0).unwrap();

        assert!(library.rotamers.is_empty());
        assert!(library.placement_info.is_empty());
    }

    #[test]
    fn handles_residue_with_empty_rotamer_list() {
        let rotamer_content = r#"
[[ALA]]
atoms = []
"#;
        let rotamer_file = create_temp_rotamer_file(rotamer_content);
        let ff = create_test_forcefield();
        let placement_registry = create_test_placement_registry();

        let library =
            RotamerLibrary::load(rotamer_file.path(), &placement_registry, &ff, 0.0).unwrap();

        assert_eq!(library.rotamers.len(), 1);
        assert_eq!(library.placement_info.len(), 1);
        let ala_rotamers = library.get_rotamers_for(ResidueType::Alanine).unwrap();
        assert_eq!(ala_rotamers.len(), 1);
        assert!(ala_rotamers[0].atoms.is_empty());
    }

    #[test]
    fn include_system_conformations_adds_unique_original_rotamer() {
        let rotamer_content = r#"
[[ALA]]
atoms = [
    { serial = 2, atom_name = "CB", partial_charge = -0.1, position = [1.5, 0.0, 0.0], force_field_type = "CB" }
]
"#;
        let rotamer_file = create_temp_rotamer_file(rotamer_content);
        let ff = create_test_forcefield();
        let placement_registry = create_test_placement_registry();

        let mut library =
            RotamerLibrary::load(rotamer_file.path(), &placement_registry, &ff, 0.0).unwrap();

        let original_count = library
            .get_rotamers_for(ResidueType::Alanine)
            .unwrap()
            .len();
        assert_eq!(original_count, 1);

        let (system, res_id) = create_test_system_with_ala(Point3::new(5.0, 5.0, 5.0));

        let mut active_residues = HashSet::new();
        active_residues.insert(res_id);

        library.include_system_conformations(&system, &active_residues, 0.1);

        let new_count = library
            .get_rotamers_for(ResidueType::Alanine)
            .unwrap()
            .len();
        assert_eq!(
            new_count,
            original_count + 1,
            "A new unique rotamer should have been added"
        );
    }

    #[test]
    fn include_system_conformations_skips_duplicate_rotamer() {
        let rotamer_content = r#"
[[ALA]]
atoms = [
    { serial = 2, atom_name = "CB", partial_charge = -0.1, position = [1.5, 0.0, 0.0], force_field_type = "CB" }
]
"#;
        let rotamer_file = create_temp_rotamer_file(rotamer_content);
        let ff = create_test_forcefield();
        let placement_registry = create_test_placement_registry();

        let mut library =
            RotamerLibrary::load(rotamer_file.path(), &placement_registry, &ff, 0.0).unwrap();

        let original_count = library
            .get_rotamers_for(ResidueType::Alanine)
            .unwrap()
            .len();
        assert_eq!(original_count, 1);

        // Create a system where the ALA side-chain is identical to the one in the library
        let (system, res_id) = create_test_system_with_ala(Point3::new(1.5, 0.0, 0.0));

        let mut active_residues = HashSet::new();
        active_residues.insert(res_id);

        library.include_system_conformations(&system, &active_residues, 0.1);

        let new_count = library
            .get_rotamers_for(ResidueType::Alanine)
            .unwrap()
            .len();
        assert_eq!(
            new_count, original_count,
            "Duplicate rotamer should not be added"
        );
    }
}
