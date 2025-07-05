use super::ids::{AtomId, ChainId};
use std::collections::HashMap;
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
        match s.trim().to_uppercase().as_str() {
            "ALA" => Ok(ResidueType::Alanine),
            "GLY" => Ok(ResidueType::Glycine),
            "ILE" => Ok(ResidueType::Isoleucine),
            "LEU" => Ok(ResidueType::Leucine),
            "PRO" => Ok(ResidueType::Proline),
            "VAL" => Ok(ResidueType::Valine),
            "PHE" => Ok(ResidueType::Phenylalanine),
            "TRP" => Ok(ResidueType::Tryptophan),
            "TYR" => Ok(ResidueType::Tyrosine),
            "ASN" => Ok(ResidueType::Asparagine),
            "CYS" => Ok(ResidueType::Cysteine),
            "GLN" => Ok(ResidueType::Glutamine),
            "SER" => Ok(ResidueType::Serine),
            "THR" => Ok(ResidueType::Threonine),
            "MET" => Ok(ResidueType::Methionine),
            "ARG" => Ok(ResidueType::Arginine),
            "LYS" => Ok(ResidueType::Lysine),
            "ASP" => Ok(ResidueType::AsparticAcid),
            "GLU" => Ok(ResidueType::GlutamicAcid),
            "HIS" => Ok(ResidueType::Histidine),
            "HSE" => Ok(ResidueType::HistidineEpsilon),
            "HSP" => Ok(ResidueType::HistidineProtonated),
            unsupported => Err(ParseResidueTypeError(unsupported.to_string())),
        }
    }
}

impl ResidueType {
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
    pub id: isize,                          // Residue sequence number from source file
    pub name: String,                       // Name of the residue (e.g., "ALA", "GLY")
    pub res_type: ResidueType,              // Type of the residue (e.g., Alanine, Glycine)
    pub chain_id: ChainId,                  // ID of the parent chain
    pub(crate) atoms: Vec<AtomId>,          // Indices of atoms belonging to this residue
    atom_name_map: HashMap<String, AtomId>, // Map from atom name to its stable ID
}

impl Residue {
    pub(crate) fn new(id: isize, name: &str, res_type: ResidueType, chain_id: ChainId) -> Self {
        Self {
            id,
            name: name.to_string(),
            res_type,
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
