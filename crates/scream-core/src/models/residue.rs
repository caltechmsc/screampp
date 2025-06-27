use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq)]
pub struct Residue {
    pub id: isize,                         // Residue sequence number from source file
    pub name: String,                      // Name of the residue (e.g., "ALA", "GLY")
    pub atom_indices: Vec<usize>,          // Indices of atoms belonging to this residue
    atom_name_map: HashMap<String, usize>, // Map from atom name to its global index
}

impl Residue {
    pub(crate) fn new(id: isize, name: &str) -> Self {
        Self {
            id,
            name: name.to_string(),
            atom_indices: Vec::new(),
            atom_name_map: HashMap::new(),
        }
    }

    pub(crate) fn add_atom(&mut self, atom_name: &str, atom_idx: usize) {
        self.atom_indices.push(atom_idx);
        self.atom_name_map.insert(atom_name.to_string(), atom_idx);
    }

    pub fn get_atom_index_by_name(&self, name: &str) -> Option<usize> {
        self.atom_name_map.get(name).copied()
    }
}
