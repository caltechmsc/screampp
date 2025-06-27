use super::residue::Residue;
use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ChainType {
    Protein,
    DNA,
    RNA,
    Other,
}

#[derive(Debug, PartialEq, Eq)]
pub struct ParseChainTypeError;

impl FromStr for ChainType {
    type Err = ParseChainTypeError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "protein" => Ok(ChainType::Protein),
            "dna" => Ok(ChainType::DNA),
            "rna" => Ok(ChainType::RNA),
            _ => Ok(ChainType::Other),
        }
    }
}

impl fmt::Display for ChainType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                ChainType::Protein => "Protein",
                ChainType::DNA => "DNA",
                ChainType::RNA => "RNA",
                ChainType::Other => "Other",
            }
        )
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Chain {
    pub id: char,                                  // Chain identifier (e.g., 'A', 'B')
    pub chain_type: ChainType,                     // Type of the chain
    pub(crate) residues: Vec<Residue>,             // List of residues in the chain
    pub(crate) residue_map: HashMap<isize, usize>, // Map from residue ID to its index in the `residues` vector
}

impl Chain {
    pub(crate) fn new(id: char, chain_type: ChainType) -> Self {
        Self {
            id,
            chain_type,
            residues: Vec::new(),
            residue_map: HashMap::new(),
        }
    }

    pub fn residues(&self) -> &[Residue] {
        &self.residues
    }

    pub fn get_residue(&self, index: usize) -> Option<&Residue> {
        self.residues.get(index)
    }

    pub fn get_residue_by_id(&self, id: isize) -> Option<&Residue> {
        self.residue_map
            .get(&id)
            .and_then(|&index| self.residues.get(index))
    }
}
