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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::residue::Residue;

    #[test]
    fn parses_valid_chain_types() {
        assert_eq!(ChainType::from_str("protein").unwrap(), ChainType::Protein);
        assert_eq!(ChainType::from_str("dna").unwrap(), ChainType::DNA);
        assert_eq!(ChainType::from_str("rna").unwrap(), ChainType::RNA);
        assert_eq!(ChainType::from_str("other").unwrap(), ChainType::Other);
    }

    #[test]
    fn parses_chain_type_case_insensitive() {
        assert_eq!(ChainType::from_str("PrOtEiN").unwrap(), ChainType::Protein);
        assert_eq!(ChainType::from_str("DNA").unwrap(), ChainType::DNA);
        assert_eq!(ChainType::from_str("rNa").unwrap(), ChainType::RNA);
    }

    #[test]
    fn parses_unknown_chain_type_as_other() {
        assert_eq!(
            ChainType::from_str("carbohydrate").unwrap(),
            ChainType::Other
        );
        assert_eq!(ChainType::from_str("").unwrap(), ChainType::Other);
        assert_eq!(ChainType::from_str("123").unwrap(), ChainType::Other);
    }

    #[test]
    fn displays_chain_type_correctly() {
        assert_eq!(ChainType::Protein.to_string(), "Protein");
        assert_eq!(ChainType::DNA.to_string(), "DNA");
        assert_eq!(ChainType::RNA.to_string(), "RNA");
        assert_eq!(ChainType::Other.to_string(), "Other");
    }

    #[test]
    fn creates_chain_with_correct_id_and_type() {
        let chain = Chain::new('A', ChainType::Protein);
        assert_eq!(chain.id, 'A');
        assert_eq!(chain.chain_type, ChainType::Protein);
        assert!(chain.residues().is_empty());
    }

    #[test]
    fn get_residue_returns_none_for_empty_chain() {
        let chain = Chain::new('B', ChainType::DNA);
        assert_eq!(chain.get_residue(0), None);
    }

    #[test]
    fn get_residue_by_id_returns_none_for_missing_residue() {
        let chain = Chain::new('C', ChainType::RNA);
        assert_eq!(chain.get_residue_by_id(10), None);
    }

    #[test]
    fn get_residue_and_get_residue_by_id_return_correct_residue() {
        let mut chain = Chain::new('D', ChainType::Protein);
        let residue = Residue::new(5, "GLY");
        chain.residues.push(residue.clone());
        chain.residue_map.insert(5, 0);

        assert_eq!(chain.get_residue(0), Some(&residue));
        assert_eq!(chain.get_residue_by_id(5), Some(&residue));
    }
}
