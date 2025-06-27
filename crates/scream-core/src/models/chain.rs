use super::residue::Residue;
use std::collections::HashMap;

pub enum ChainType {
    Protein, // Peptides and Proteins
    DNA,     // Deoxyribonucleic Acid
    RNA,     // Ribonucleic Acid
    Other,   // Any other type of chain (e.g., carbohydrates, lipids, etc.)
}

impl ChainType {
    pub fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "protein" => ChainType::Protein,
            "dna" => ChainType::DNA,
            "rna" => ChainType::RNA,
            _ => ChainType::Other,
        }
    }

    pub fn to_str(&self) -> &str {
        match self {
            ChainType::Protein => "Protein",
            ChainType::DNA => "DNA",
            ChainType::RNA => "RNA",
            ChainType::Other => "Other",
        }
    }
}

pub struct Chain {
    pub id: char,                       // Chain identifier (e.g., 'A', 'B', etc.)
    pub chain_type: ChainType,          // Type of the chain (Protein, DNA, RNA, Other)
    residues: Vec<Residue>,             // List of residues in the chain
    residue_map: HashMap<usize, usize>, // Map from residue ID to index in the residues vector
}

impl Chain {
    pub fn new(id: char, chain_type: ChainType, residues: Vec<Residue>) -> Self {
        let residues_vec = residues;
        let mut residue_map = HashMap::new();
        for (i, residue) in residues_vec.iter().enumerate() {
            residue_map.insert(residue.id, i);
        }
        Self {
            id,
            chain_type,
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
