use super::ids::ResidueId;
use std::fmt;
use std::str::FromStr;
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ChainType {
    Protein,
    DNA,
    RNA,
    Ligand,
    Water,
    Other,
}

#[derive(Debug, Error)]
#[error("Invalid chain type string")]
pub struct ParseChainTypeError;

impl FromStr for ChainType {
    type Err = ParseChainTypeError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "protein" => Ok(ChainType::Protein),
            "dna" => Ok(ChainType::DNA),
            "rna" => Ok(ChainType::RNA),
            "ligand" => Ok(ChainType::Ligand),
            "water" => Ok(ChainType::Water),
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
                ChainType::Ligand => "Ligand",
                ChainType::Water => "Water",
                ChainType::Other => "Other",
            }
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chain {
    pub id: char,                        // Chain identifier (e.g., 'A', 'B')
    pub chain_type: ChainType,           // Type of the chain
    pub(crate) residues: Vec<ResidueId>, // Ordered list of residue IDs belonging to this chain
}

impl Chain {
    pub(crate) fn new(id: char, chain_type: ChainType) -> Self {
        Self {
            id,
            chain_type,
            residues: Vec::new(),
        }
    }

    pub fn residues(&self) -> &[ResidueId] {
        &self.residues
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_valid_chain_types() {
        assert_eq!(ChainType::from_str("protein").unwrap(), ChainType::Protein);
        assert_eq!(ChainType::from_str("dna").unwrap(), ChainType::DNA);
        assert_eq!(ChainType::from_str("water").unwrap(), ChainType::Water);
    }

    #[test]
    fn parses_chain_type_case_insensitive() {
        assert_eq!(ChainType::from_str("PrOtEiN").unwrap(), ChainType::Protein);
    }

    #[test]
    fn parses_unknown_chain_type_as_other() {
        assert_eq!(
            ChainType::from_str("carbohydrate").unwrap(),
            ChainType::Other
        );
    }

    #[test]
    fn displays_chain_type_correctly() {
        assert_eq!(ChainType::Protein.to_string(), "Protein");
    }

    #[test]
    fn creates_chain_with_correct_id_and_type() {
        let chain = Chain::new('A', ChainType::Protein);
        assert_eq!(chain.id, 'A');
        assert_eq!(chain.chain_type, ChainType::Protein);
        assert!(chain.residues().is_empty());
    }
}
