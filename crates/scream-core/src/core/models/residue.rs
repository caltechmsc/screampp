use super::ids::{AtomId, ChainId};
use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ResidueType {
    // --- Aliphatic, Nonpolar ---
    Alanine,    // Alanine (ALA)
    Glycine,    // Glycine (GLY)
    Isoleucine, // Isoleucine (ILE)
    Leucine,    // Leucine (LEU)
    Proline,    // Proline (PRO)
    Valine,     // Valine (VAL)

    // --- Aromatic ---
    Phenylalanine, // Phenylalanine (PHE)
    Tryptophan,    // Tryptophan (TRP)
    Tyrosine,      // Tyrosine (TYR)

    // --- Polar, Uncharged ---
    Asparagine, // Asparagine (ASN)
    Cysteine,   // Cysteine (CYS)
    Glutamine,  // Glutamine (GLN)
    Serine,     // Serine (SER)
    Threonine,  // Threonine (THR)
    Methionine, // Methionine (MET)

    // --- Positively Charged (Basic) ---
    Arginine, // Arginine (ARG)
    Lysine,   // Lysine (LYS)

    // --- Negatively Charged (Acidic) ---
    AsparticAcid, // Aspartic Acid (ASP)
    GlutamicAcid, // Glutamic Acid (GLU)

    // --- Special Case: Histidine and its Variants ---
    Histidine, // Histidine (HIS) - Typically assumed to be the Epsilon-protonated state (HSE)
    HistidineEpsilon, // Epsilon-protonated Histidine (HSE) - An alias for `Histidine`
    HistidineProtonated, // Doubly-protonated Histidine (HSP) - The positively charged variant
}

#[derive(Debug, Error, PartialEq, Eq)]
#[error("Unsupported or unknown three-letter residue code: '{0}'")]
pub struct ParseResidueTypeError(pub String);

impl FromStr for ResidueType {
    type Err = ParseResidueTypeError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        ResidueType::from_str_optional(s)
            .ok_or_else(|| ParseResidueTypeError(s.trim().to_uppercase()))
    }
}

impl fmt::Display for ResidueType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_three_letter())
    }
}

impl ResidueType {
    fn parse_code(s: &str) -> Option<Self> {
        match s.trim().to_uppercase().as_str() {
            "ALA" => Some(ResidueType::Alanine),
            "GLY" => Some(ResidueType::Glycine),
            "ILE" => Some(ResidueType::Isoleucine),
            "LEU" => Some(ResidueType::Leucine),
            "PRO" => Some(ResidueType::Proline),
            "VAL" => Some(ResidueType::Valine),
            "PHE" => Some(ResidueType::Phenylalanine),
            "TRP" => Some(ResidueType::Tryptophan),
            "TYR" => Some(ResidueType::Tyrosine),
            "ASN" => Some(ResidueType::Asparagine),
            "CYS" => Some(ResidueType::Cysteine),
            "GLN" => Some(ResidueType::Glutamine),
            "SER" => Some(ResidueType::Serine),
            "THR" => Some(ResidueType::Threonine),
            "MET" => Some(ResidueType::Methionine),
            "ARG" => Some(ResidueType::Arginine),
            "LYS" => Some(ResidueType::Lysine),
            "ASP" => Some(ResidueType::AsparticAcid),
            "GLU" => Some(ResidueType::GlutamicAcid),
            "HIS" => Some(ResidueType::Histidine),
            "HSE" => Some(ResidueType::HistidineEpsilon),
            "HSP" => Some(ResidueType::HistidineProtonated),
            _ => None,
        }
    }

    pub fn from_str_optional(s: &str) -> Option<Self> {
        Self::parse_code(s)
    }

    pub fn to_three_letter(self) -> &'static str {
        match self {
            ResidueType::Alanine => "ALA",
            ResidueType::Glycine => "GLY",
            ResidueType::Isoleucine => "ILE",
            ResidueType::Leucine => "LEU",
            ResidueType::Proline => "PRO",
            ResidueType::Valine => "VAL",
            ResidueType::Phenylalanine => "PHE",
            ResidueType::Tryptophan => "TRP",
            ResidueType::Tyrosine => "TYR",
            ResidueType::Asparagine => "ASN",
            ResidueType::Cysteine => "CYS",
            ResidueType::Glutamine => "GLN",
            ResidueType::Serine => "SER",
            ResidueType::Threonine => "THR",
            ResidueType::Methionine => "MET",
            ResidueType::Arginine => "ARG",
            ResidueType::Lysine => "LYS",
            ResidueType::AsparticAcid => "ASP",
            ResidueType::GlutamicAcid => "GLU",
            ResidueType::Histidine => "HIS",
            ResidueType::HistidineEpsilon => "HSE",
            ResidueType::HistidineProtonated => "HSP",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Residue {
    pub residue_number: isize, // Residue sequence number from source file
    pub name: String,          // Name of the residue (e.g., "ALA", "GLY")
    pub residue_type: Option<ResidueType>, // Optional residue type (e.g., Alanine, Glycine)
    pub chain_id: ChainId,     // ID of the parent chain
    pub(crate) atoms: Vec<AtomId>, // Indices of atoms belonging to this residue
    atom_name_map: HashMap<String, Vec<AtomId>>, // Map from atom names to their IDs
}

impl Residue {
    pub(crate) fn new(
        residue_number: isize,
        name: &str,
        residue_type: Option<ResidueType>,
        chain_id: ChainId,
    ) -> Self {
        Self {
            residue_number,
            name: name.to_string(),
            residue_type,
            chain_id,
            atoms: Vec::new(),
            atom_name_map: HashMap::new(),
        }
    }

    pub(crate) fn add_atom(&mut self, atom_name: &str, atom_id: AtomId) {
        self.atoms.push(atom_id);
        self.atom_name_map
            .entry(atom_name.to_string())
            .or_default()
            .push(atom_id);
    }

    pub(crate) fn remove_atom(&mut self, atom_name: &str, atom_id_to_remove: AtomId) {
        self.atoms.retain(|&id| id != atom_id_to_remove);

        if let Some(ids) = self.atom_name_map.get_mut(atom_name) {
            ids.retain(|&id| id != atom_id_to_remove);
            if ids.is_empty() {
                self.atom_name_map.remove(atom_name);
            }
        }
    }

    pub fn atoms(&self) -> &[AtomId] {
        &self.atoms
    }

    pub fn get_atom_id_by_name(&self, name: &str) -> Option<AtomId> {
        self.atom_name_map.get(name).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::ids::{AtomId, ChainId};
    use slotmap::KeyData;

    fn dummy_atom_id(n: u64) -> AtomId {
        AtomId::from(KeyData::from_ffi(n))
    }

    fn dummy_chain_id(n: u64) -> ChainId {
        ChainId::from(KeyData::from_ffi(n))
    }

    #[test]
    fn parses_valid_three_letter_residue_codes() {
        assert_eq!(ResidueType::from_str("ALA").unwrap(), ResidueType::Alanine);
        assert_eq!(ResidueType::from_str("gly").unwrap(), ResidueType::Glycine);
        assert_eq!(
            ResidueType::from_str("HSE").unwrap(),
            ResidueType::HistidineEpsilon
        );
        assert_eq!(
            ResidueType::from_str("hsp").unwrap(),
            ResidueType::HistidineProtonated
        );
    }

    #[test]
    fn fails_to_parse_invalid_residue_code() {
        let err = ResidueType::from_str("XXX").unwrap_err();
        assert_eq!(err, ParseResidueTypeError("XXX".to_string()));
    }

    #[test]
    fn from_str_optional_returns_none_for_unknown_code() {
        assert!(ResidueType::from_str_optional("XYZ").is_none());
        assert!(ResidueType::from_str_optional("unknown").is_none());
    }

    #[test]
    fn display_formats_residuetype_as_three_letter_code() {
        assert_eq!(format!("{}", ResidueType::Alanine), "ALA");
        assert_eq!(format!("{}", ResidueType::HistidineProtonated), "HSP");
        assert_eq!(format!("{}", ResidueType::Glycine), "GLY");
    }

    #[test]
    fn to_three_letter_returns_correct_code() {
        assert_eq!(ResidueType::Alanine.to_three_letter(), "ALA");
        assert_eq!(ResidueType::HistidineProtonated.to_three_letter(), "HSP");
    }

    #[test]
    fn residue_adds_and_retrieves_atoms_by_name() {
        let mut residue = Residue::new(1, "ALA", Some(ResidueType::Alanine), dummy_chain_id(0));
        let atom_id = dummy_atom_id(42);
        residue.add_atom("CA", atom_id);
        assert_eq!(residue.atoms(), &[atom_id]);
        assert_eq!(residue.get_atom_id_by_name("CA"), Some(atom_id));
    }

    #[test]
    fn residue_removes_atom_by_name_and_id() {
        let mut residue = Residue::new(2, "GLY", Some(ResidueType::Glycine), dummy_chain_id(1));
        let atom_id1 = dummy_atom_id(1);
        let atom_id2 = dummy_atom_id(2);
        residue.add_atom("N", atom_id1);
        residue.add_atom("CA", atom_id2);
        residue.remove_atom("N", atom_id1);
        assert_eq!(residue.atoms(), &[atom_id2]);
        assert_eq!(residue.get_atom_id_by_name("N"), None);
        assert_eq!(residue.get_atom_id_by_name("CA"), Some(atom_id2));
    }

    #[test]
    fn residue_get_atom_id_by_name_returns_none_for_unknown_atom() {
        let residue = Residue::new(3, "SER", Some(ResidueType::Serine), dummy_chain_id(2));
        assert_eq!(residue.get_atom_id_by_name("CB"), None);
    }

    #[test]
    fn residue_atoms_returns_empty_slice_when_no_atoms() {
        let residue = Residue::new(4, "THR", Some(ResidueType::Threonine), dummy_chain_id(3));
        assert!(residue.atoms().is_empty());
    }
}
