use super::rotamer::{Rotamer, RotamerData};
use crate::core::{
    forcefield::{
        parameterization::{ParameterizationError, Parameterizer},
        params::Forcefield,
    },
    models::{
        atom::Atom,
        ids::{AtomId, ResidueId},
        residue::ResidueType,
        system::MolecularSystem,
    },
    topology::registry::TopologyRegistry,
};
use nalgebra::Point3;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::str::FromStr;
use thiserror::Error;

type RawRotamerFile = HashMap<String, Vec<RotamerData>>;

#[derive(Debug, Default, Clone)]
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
        "Missing topology definition for residue type '{0}', which is present in the rotamer library"
    )]
    MissingTopology(String),
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
        topology_registry: &TopologyRegistry,
        forcefield: &Forcefield,
        delta_s_factor: f64,
    ) -> Result<Self, LibraryLoadError> {
        // --- Phase 1: Load raw rotamer data from TOML ---
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

        // --- Phase 2: Create a parameterizer for pre-parameterizing rotamers ---
        let parameterizer = Parameterizer::new(forcefield, topology_registry, delta_s_factor);
        let mut final_rotamers_map = HashMap::new();

        // --- Phase 3: Process and parameterize each rotamer ---
        for (res_name, raw_rotamer_list) in raw_lib {
            let residue_type = ResidueType::from_str(&res_name)
                .map_err(|_| LibraryLoadError::UnknownResidueType(res_name.clone()))?;

            let topology = topology_registry
                .get(&res_name)
                .ok_or_else(|| LibraryLoadError::MissingTopology(res_name.clone()))?;

            let mut processed_rotamers = Vec::with_capacity(raw_rotamer_list.len());
            for raw_rotamer_data in raw_rotamer_list {
                let rotamer = Self::process_raw_rotamer(
                    &raw_rotamer_data,
                    &parameterizer,
                    &res_name,
                    topology,
                    rotamer_toml_path,
                )?;
                processed_rotamers.push(rotamer);
            }
            final_rotamers_map.insert(residue_type, processed_rotamers);
        }

        Ok(Self {
            rotamers: final_rotamers_map,
        })
    }

    fn process_raw_rotamer(
        raw_rotamer_data: &RotamerData,
        parameterizer: &Parameterizer,
        res_name: &str,
        topology: &crate::core::topology::registry::ResidueTopology,
        path_for_error: &Path,
    ) -> Result<Rotamer, LibraryLoadError> {
        let mut atoms = Vec::with_capacity(raw_rotamer_data.atoms.len());
        let mut serial_to_index_map = HashMap::with_capacity(raw_rotamer_data.atoms.len());
        let placeholder_residue_id = ResidueId::default();

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
                .parameterize_protein_atom(&mut atom, res_name, topology)
                .map_err(|e| LibraryLoadError::Parameterization {
                    path: path_for_error.to_string_lossy().to_string(),
                    residue_type: res_name.to_string(),
                    source: e,
                })?;

            atoms.push(atom);
        }

        let mut bonds = Vec::new();
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

    pub fn include_system_conformations(
        &mut self,
        system: &MolecularSystem,
        active_residues: &HashSet<ResidueId>,
        topology_registry: &TopologyRegistry,
    ) {
        for &residue_id in active_residues {
            if let Some(rotamer) =
                self.extract_rotamer_from_system(system, residue_id, topology_registry)
            {
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
        topology_registry: &TopologyRegistry,
    ) -> Option<Rotamer> {
        let residue = system.residue(residue_id)?;
        let topology = topology_registry.get(&residue.name)?;

        let mut atom_pool: HashMap<String, Vec<AtomId>> = HashMap::new();
        for &atom_id in residue.atoms() {
            let atom = system.atom(atom_id)?;
            atom_pool
                .entry(atom.name.clone())
                .or_default()
                .push(atom_id);
        }

        let mut new_rotamer_atoms = Vec::new();
        let mut old_id_to_new_index = HashMap::new();

        let all_atom_names: HashSet<_> = topology
            .anchor_atoms
            .iter()
            .chain(&topology.sidechain_atoms)
            .collect();

        for name in all_atom_names {
            if let Some(ids) = atom_pool.get_mut(name) {
                if let Some(atom_id) = ids.pop() {
                    let atom = system.atom(atom_id)?;
                    old_id_to_new_index.insert(atom_id, new_rotamer_atoms.len());
                    new_rotamer_atoms.push(atom.clone());
                }
            }
        }

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
