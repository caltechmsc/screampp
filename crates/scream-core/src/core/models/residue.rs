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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creates_residue_with_correct_id_and_name() {
        let residue = Residue::new(42, "GLY");
        assert_eq!(residue.id, 42);
        assert_eq!(residue.name, "GLY");
        assert!(residue.atom_indices.is_empty());
    }

    #[test]
    fn adds_atom_and_retrieves_index_by_name() {
        let mut residue = Residue::new(1, "ALA");
        residue.add_atom("CA", 5);
        residue.add_atom("CB", 6);

        assert_eq!(residue.atom_indices, vec![5, 6]);
        assert_eq!(residue.get_atom_index_by_name("CA"), Some(5));
        assert_eq!(residue.get_atom_index_by_name("CB"), Some(6));
    }

    #[test]
    fn get_atom_index_by_name_returns_none_for_missing_atom() {
        let mut residue = Residue::new(2, "SER");
        residue.add_atom("OG", 10);

        assert_eq!(residue.get_atom_index_by_name("CA"), None);
        assert_eq!(residue.get_atom_index_by_name("OG"), Some(10));
    }

    #[test]
    fn add_atom_overwrites_existing_atom_name() {
        let mut residue = Residue::new(3, "THR");
        residue.add_atom("CG2", 20);
        residue.add_atom("CG2", 21);

        assert_eq!(residue.atom_indices, vec![20, 21]);
        assert_eq!(residue.get_atom_index_by_name("CG2"), Some(21));
    }

    #[test]
    fn add_atom_with_empty_name_and_retrieve() {
        let mut residue = Residue::new(4, "VAL");
        residue.add_atom("", 30);

        assert_eq!(residue.get_atom_index_by_name(""), Some(30));
    }
}
