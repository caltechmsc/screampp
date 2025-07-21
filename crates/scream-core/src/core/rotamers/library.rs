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
        placement_registry_path: &Path,
        forcefield: &Forcefield,
        delta_s_factor: f64,
    ) -> Result<Self, LibraryLoadError> {
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
                parameterized_rotamers.push(Rotamer { atoms });
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

// TODO: Unit tests for RotamerLibrary
