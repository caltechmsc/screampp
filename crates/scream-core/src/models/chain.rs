use super::residue::Residue;
use std::collections::HashMap;

pub struct Chain {
    pub id: char,                       // Chain identifier (e.g., 'A', 'B', etc.)
    residues: Vec<Residue>,             // List of residues in the chain
    residue_map: HashMap<usize, usize>, // Map from residue ID to index in the residues vector
}

impl Chain {
    pub fn new(id: char, residues: Vec<Residue>) -> Self {
        let residues_vec = residues;
        let mut residue_map = HashMap::new();
        for (i, residue) in residues_vec.iter().enumerate() {
            residue_map.insert(residue.id, i);
        }
        Self {
            id,
            residues: residues_vec,
            residue_map,
        }
    }

    pub fn add_residue(&mut self, residue: Residue) {
        let index = self.residues.len();
        self.residue_map.insert(residue.id, index);
        self.residues.push(residue);
    }

    pub fn get_residue(&self, index: usize) -> Option<&Residue> {
        self.residues.get(index)
    }

    pub fn get_residue_by_id(&self, id: usize) -> Option<&Residue> {
        self.residue_map
            .get(&id)
            .and_then(|&index| self.residues.get(index))
    }

    pub fn residues(&self) -> &[Residue] {
        &self.residues
    }
}
