use super::rotamer::{Rotamer, RotamerData};
use crate::core::forcefield::parameterization::{ParameterizationError, Parameterizer};
use crate::core::forcefield::params::Forcefield;
use crate::core::models::atom::Atom;
use crate::core::models::ids::ResidueId;
use crate::core::models::residue::ResidueType;
use nalgebra::Point3;
use std::collections::HashMap;
use std::path::Path;
use std::str::FromStr;
use thiserror::Error;

type RawRotamerFile = HashMap<String, Vec<RotamerData>>;

#[derive(Debug, Default)]
pub struct RotamerLibrary {
    pub rotamers: HashMap<ResidueType, Vec<Rotamer>>,
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
        library_path: &Path,
        forcefield: &Forcefield,
        delta_s_factor: f64,
    ) -> Result<Self, LibraryLoadError> {
        let content = std::fs::read_to_string(library_path).map_err(|e| LibraryLoadError::Io {
            path: library_path.to_string_lossy().to_string(),
            source: e,
        })?;

        let raw_lib: RawRotamerFile =
            toml::from_str(&content).map_err(|e| LibraryLoadError::Toml {
                path: library_path.to_string_lossy().to_string(),
                source: e,
            })?;

        let parameterizer = Parameterizer::new(forcefield.clone(), delta_s_factor);
        let mut final_rotamers_map = HashMap::new();

        for (res_name, raw_rotamer_list) in raw_lib {
            let res_type = ResidueType::from_str(&res_name)
                .map_err(|_| LibraryLoadError::UnknownResidueType(res_name.clone()))?;

            let mut parameterized_rotamers = Vec::new();

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

                    parameterizer
                        .parameterize_atom(&mut atom, &res_name)
                        .map_err(|e| LibraryLoadError::Parameterization {
                            path: library_path.to_string_lossy().to_string(),
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
        })
    }
}
