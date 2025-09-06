use super::ids::ResidueId;
use std::fmt;
use std::str::FromStr;
use thiserror::Error;

/// Represents the type of a molecular chain in a structure.
///
/// This enum categorizes chains based on their molecular composition,
/// which is useful for algorithms that need to distinguish between
/// different types of molecules in simulations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ChainType {
    /// Protein chain.
    Protein,
    /// DNA chain.
    DNA,
    /// RNA chain.
    RNA,
    /// Ligand or small molecule chain.
    Ligand,
    /// Water molecule chain.
    Water,
    /// Other or unspecified chain type.
    Other,
}

/// Error type for failed parsing of chain type strings.
///
/// This error is returned when attempting to parse an invalid
/// string into a `ChainType`. Note that this error is currently
/// not used since unknown strings default to `ChainType::Other`.
#[derive(Debug, Error)]
#[error("Invalid chain type string")]
pub struct ParseChainTypeError;

impl FromStr for ChainType {
    type Err = ParseChainTypeError;

    /// Parses a string into a `ChainType`.
    ///
    /// This implementation converts string representations of chain types
    /// into the corresponding enum variants. It is case-insensitive and
    /// defaults to `ChainType::Other` for unknown strings.
    ///
    /// # Arguments
    ///
    /// * `s` - The string to parse.
    ///
    /// # Return
    ///
    /// Returns the parsed `ChainType`.
    ///
    /// # Errors
    ///
    /// This method does not currently return errors; unknown strings
    /// are mapped to `ChainType::Other`.
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
    /// Formats the `ChainType` as a human-readable string.
    ///
    /// This implementation allows `ChainType` to be displayed as a string
    /// using capitalized names for each variant.
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

/// Represents a molecular chain in a structure.
///
/// This struct encapsulates the properties and residues of a single chain,
/// providing access to its constituent residues in order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chain {
    /// The single-character identifier of the chain (e.g., 'A', 'B').
    pub id: char,
    /// The type of the chain.
    pub chain_type: ChainType,
    /// Ordered list of residue IDs belonging to this chain.
    pub(crate) residues: Vec<ResidueId>,
}

impl Chain {
    /// Creates a new `Chain` with the specified ID and type.
    ///
    /// This constructor initializes a chain with an empty list of residues.
    ///
    /// # Arguments
    ///
    /// * `id` - The single-character identifier for the chain.
    /// * `chain_type` - The type of the chain.
    pub(crate) fn new(id: char, chain_type: ChainType) -> Self {
        Self {
            id,
            chain_type,
            residues: Vec::new(),
        }
    }

    /// Returns a slice of all residue IDs in the chain.
    ///
    /// This provides read-only access to the ordered list of residues
    /// belonging to the chain.
    ///
    /// # Return
    ///
    /// A slice containing all residue IDs in order.
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
