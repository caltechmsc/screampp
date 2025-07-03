use super::ids::{AtomId, ChainId};
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Residue {
    pub id: isize,                          // Residue sequence number from source file
    pub name: String,                       // Name of the residue (e.g., "ALA", "GLY")
    pub chain_id: ChainId,                  // ID of the parent chain
    pub(crate) atoms: Vec<AtomId>,          // Indices of atoms belonging to this residue
    atom_name_map: HashMap<String, AtomId>, // Map from atom name to its stable ID
}

impl Residue {
    pub(crate) fn new(id: isize, name: &str, chain_id: ChainId) -> Self {
        Self {
            id,
            name: name.to_string(),
            chain_id,
            atoms: Vec::new(),
            atom_name_map: HashMap::new(),
        }
    }

    pub(crate) fn add_atom(&mut self, atom_name: &str, atom_id: AtomId) {
        self.atoms.push(atom_id);
        self.atom_name_map.insert(atom_name.to_string(), atom_id);
    }

    pub(crate) fn remove_atom(&mut self, atom_name: &str, atom_id: AtomId) {
        self.atoms.retain(|&id| id != atom_id);
        self.atom_name_map.remove(atom_name);
    }

    pub fn atoms(&self) -> &[AtomId] {
        &self.atoms
    }

    pub fn get_atom_id_by_name(&self, name: &str) -> Option<AtomId> {
        self.atom_name_map.get(name).copied()
    }
}
