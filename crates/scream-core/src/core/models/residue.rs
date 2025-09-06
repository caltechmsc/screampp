use super::atom::AtomRole;
use super::ids::{AtomId, ChainId};
use super::system::MolecularSystem;
use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;
use thiserror::Error;

/// Represents the type of an amino acid residue in a protein structure.
///
/// This enum defines standard amino acid types, including their protonation states
/// for histidine variants, used in molecular modeling and simulations.
/// Each variant corresponds to a three-letter code commonly used in PDB files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ResidueType {
    /// Alanine (ALA).
    Alanine,
    /// Glycine (GLY).
    Glycine,
    /// Isoleucine (ILE).
    Isoleucine,
    /// Leucine (LEU).
    Leucine,
    /// Proline (PRO).
    Proline,
    /// Valine (VAL).
    Valine,
    /// Phenylalanine (PHE).
    Phenylalanine,
    /// Tryptophan (TRP).
    Tryptophan,
    /// Tyrosine (TYR).
    Tyrosine,
    /// Asparagine (ASN).
    Asparagine,
    /// Cysteine (CYS).
    Cysteine,
    /// Glutamine (GLN).
    Glutamine,
    /// Serine (SER).
    Serine,
    /// Threonine (THR).
    Threonine,
    /// Methionine (MET).
    Methionine,
    /// Arginine (ARG).
    Arginine,
    /// Lysine (LYS).
    Lysine,
    /// Aspartic Acid (ASP).
    AsparticAcid,
    /// Glutamic Acid (GLU).
    GlutamicAcid,
    /// Histidine (HIS), typically epsilon-protonated.
    Histidine,
    /// Epsilon-protonated Histidine (HSE).
    HistidineEpsilon,
    /// Doubly-protonated Histidine (HSP).
    HistidineProtonated,
}

/// Error type for failed parsing of residue type strings.
///
/// This error is returned when attempting to parse an invalid or unsupported
/// three-letter residue code into a `ResidueType`.
#[derive(Debug, Error, PartialEq, Eq)]
#[error("Unsupported or unknown three-letter residue code: '{0}'")]
pub struct ParseResidueTypeError(pub String);

impl FromStr for ResidueType {
    type Err = ParseResidueTypeError;

    /// Parses a string into a `ResidueType`.
    ///
    /// This implementation converts a three-letter residue code (case-insensitive)
    /// into the corresponding `ResidueType` variant. It supports standard amino acids
    /// and histidine variants.
    ///
    /// # Arguments
    ///
    /// * `s` - The string to parse, expected to be a three-letter code.
    ///
    /// # Return
    ///
    /// Returns the parsed `ResidueType` if the code is recognized.
    ///
    /// # Errors
    ///
    /// Returns `ParseResidueTypeError` if the code is invalid or unsupported.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        ResidueType::from_str_optional(s)
            .ok_or_else(|| ParseResidueTypeError(s.trim().to_uppercase()))
    }
}

impl fmt::Display for ResidueType {
    /// Formats the `ResidueType` as its three-letter code.
    ///
    /// This implementation allows `ResidueType` to be displayed as a string
    /// using the standard three-letter amino acid codes.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_three_letter())
    }
}

impl ResidueType {
    /// Parses a three-letter residue code into a `ResidueType`.
    ///
    /// This is an internal helper function that performs case-insensitive matching
    /// of standard three-letter codes to `ResidueType` variants.
    ///
    /// # Arguments
    ///
    /// * `s` - The string to parse.
    ///
    /// # Return
    ///
    /// Returns `Some(ResidueType)` if the code is recognized, otherwise `None`.
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

    /// Attempts to parse a string into a `ResidueType` without error handling.
    ///
    /// This method provides a fallible alternative to `FromStr::from_str`,
    /// returning `None` for invalid codes instead of an error.
    ///
    /// # Arguments
    ///
    /// * `s` - The string to parse.
    ///
    /// # Return
    ///
    /// Returns `Some(ResidueType)` if parsing succeeds, otherwise `None`.
    pub fn from_str_optional(s: &str) -> Option<Self> {
        Self::parse_code(s)
    }

    /// Converts the `ResidueType` to its standard three-letter code.
    ///
    /// This method returns the canonical three-letter abbreviation for the amino acid,
    /// as used in PDB files and molecular databases.
    ///
    /// # Return
    ///
    /// A static string slice containing the three-letter code.
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

/// Represents a residue in a molecular structure.
///
/// This struct encapsulates the properties and atoms of a single residue,
/// providing efficient access to backbone and sidechain atoms through caching.
/// It is used in protein modeling and side-chain placement algorithms.
#[derive(Debug, Clone, PartialEq)]
pub struct Residue {
    /// The sequential number of the residue in its chain.
    pub residue_number: isize,
    /// The name of the residue.
    pub name: String,
    /// The type of the residue, if known.
    pub residue_type: Option<ResidueType>,
    /// The ID of the chain this residue belongs to.
    pub chain_id: ChainId,
    /// List of atom IDs belonging to this residue.
    atoms: Vec<AtomId>,
    /// Mapping from atom names to lists of atom IDs for quick lookup.
    atom_name_map: HashMap<String, Vec<AtomId>>,
    /// Cached list of sidechain atom IDs.
    sidechain_atoms_cache: Vec<AtomId>,
    /// Cached list of backbone atom IDs.
    backbone_atoms_cache: Vec<AtomId>,
}

impl Residue {
    /// Creates a new `Residue` with the specified properties.
    ///
    /// This constructor initializes a residue with empty atom lists and caches.
    /// Atoms can be added later using `add_atom`.
    ///
    /// # Arguments
    ///
    /// * `residue_number` - The sequential number of the residue.
    /// * `name` - The name of the residue.
    /// * `residue_type` - The type of the residue, if known.
    /// * `chain_id` - The ID of the chain this residue belongs to.
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
            sidechain_atoms_cache: Vec::new(),
            backbone_atoms_cache: Vec::new(),
        }
    }

    /// Clears the cached lists of backbone and sidechain atoms.
    ///
    /// This method is called whenever the atom list changes to ensure
    /// that cached data remains consistent.
    fn invalidate_caches(&mut self) {
        self.sidechain_atoms_cache.clear();
        self.backbone_atoms_cache.clear();
    }

    /// Adds an atom to the residue.
    ///
    /// This method registers an atom with the given name and ID, updating
    /// the internal lists and invalidating caches as necessary.
    ///
    /// # Arguments
    ///
    /// * `atom_name` - The name of the atom.
    /// * `atom_id` - The ID of the atom to add.
    pub(crate) fn add_atom(&mut self, atom_name: &str, atom_id: AtomId) {
        self.atoms.push(atom_id);
        self.atom_name_map
            .entry(atom_name.to_string())
            .or_default()
            .push(atom_id);
        self.invalidate_caches();
    }

    /// Removes an atom from the residue.
    ///
    /// This method removes the specified atom by name and ID, cleaning up
    /// the internal data structures and invalidating caches.
    ///
    /// # Arguments
    ///
    /// * `atom_name` - The name of the atom to remove.
    /// * `atom_id_to_remove` - The ID of the atom to remove.
    pub(crate) fn remove_atom(&mut self, atom_name: &str, atom_id_to_remove: AtomId) {
        self.atoms.retain(|&id| id != atom_id_to_remove);

        if let Some(ids) = self.atom_name_map.get_mut(atom_name) {
            ids.retain(|&id| id != atom_id_to_remove);
            if ids.is_empty() {
                self.atom_name_map.remove(atom_name);
            }
        }
        self.invalidate_caches();
    }

    /// Returns a slice of all atom IDs in the residue.
    ///
    /// This provides read-only access to the list of atoms belonging to the residue.
    ///
    /// # Return
    ///
    /// A slice containing all atom IDs.
    pub fn atoms(&self) -> &[AtomId] {
        &self.atoms
    }

    /// Returns a slice of sidechain atom IDs, building the cache if necessary.
    ///
    /// This method lazily computes and caches the list of sidechain atoms
    /// by querying the molecular system for atom roles.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system containing the atoms.
    ///
    /// # Return
    ///
    /// A slice containing sidechain atom IDs.
    pub fn sidechain_atoms<'a>(&'a mut self, system: &'a MolecularSystem) -> &'a [AtomId] {
        if self.sidechain_atoms_cache.is_empty() && !self.atoms.is_empty() {
            self.build_caches(system);
        }
        &self.sidechain_atoms_cache
    }

    /// Returns a slice of backbone atom IDs, building the cache if necessary.
    ///
    /// This method lazily computes and caches the list of backbone atoms
    /// by querying the molecular system for atom roles.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system containing the atoms.
    ///
    /// # Return
    ///
    /// A slice containing backbone atom IDs.
    pub fn backbone_atoms<'a>(&'a mut self, system: &'a MolecularSystem) -> &'a [AtomId] {
        if self.backbone_atoms_cache.is_empty() && !self.atoms.is_empty() {
            self.build_caches(system);
        }
        &self.backbone_atoms_cache
    }

    /// Builds the caches for backbone and sidechain atoms.
    ///
    /// This method iterates through all atoms in the residue, determines their roles
    /// using the molecular system, and populates the respective caches.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system for querying atom roles.
    fn build_caches(&mut self, system: &MolecularSystem) {
        self.invalidate_caches();
        for &atom_id in &self.atoms {
            if let Some(atom) = system.atom(atom_id) {
                match atom.role {
                    AtomRole::Sidechain => self.sidechain_atoms_cache.push(atom_id),
                    AtomRole::Backbone => self.backbone_atoms_cache.push(atom_id),
                    _ => {}
                }
            }
        }
    }

    /// Retrieves atom IDs by atom name.
    ///
    /// This method looks up atoms with the specified name in the residue.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the atom to search for.
    ///
    /// # Return
    ///
    /// Returns `Some` slice of atom IDs if found, otherwise `None`.
    pub fn get_atom_ids_by_name(&self, name: &str) -> Option<&[AtomId]> {
        self.atom_name_map.get(name).map(|v| v.as_slice())
    }

    /// Retrieves the first atom ID by atom name.
    ///
    /// This method returns the first atom with the specified name, useful
    /// when there is expected to be only one atom of that name.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the atom to search for.
    ///
    /// # Return
    ///
    /// Returns `Some` atom ID if found, otherwise `None`.
    pub fn get_first_atom_id_by_name(&self, name: &str) -> Option<AtomId> {
        self.get_atom_ids_by_name(name)
            .and_then(|ids| ids.first().copied())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::atom::Atom;
    use crate::core::models::ids::{AtomId, ChainId};
    use crate::core::models::system::MolecularSystem;
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
    }

    #[test]
    fn display_formats_residuetype_as_three_letter_code() {
        assert_eq!(format!("{}", ResidueType::Alanine), "ALA");
    }

    #[test]
    fn to_three_letter_returns_correct_code() {
        assert_eq!(ResidueType::Alanine.to_three_letter(), "ALA");
    }

    #[test]
    fn residue_adds_and_retrieves_single_atom() {
        let mut residue = Residue::new(1, "ALA", Some(ResidueType::Alanine), dummy_chain_id(0));
        let atom_id = dummy_atom_id(42);
        residue.add_atom("CA", atom_id);

        assert_eq!(residue.atoms(), &[atom_id]);
        assert_eq!(
            residue.get_atom_ids_by_name("CA"),
            Some([atom_id].as_slice())
        );
        assert_eq!(residue.get_first_atom_id_by_name("CA"), Some(atom_id));
    }

    #[test]
    fn residue_handles_duplicate_atom_names_correctly() {
        let mut residue = Residue::new(1, "ALA", Some(ResidueType::Alanine), dummy_chain_id(0));
        let hcb1_id = dummy_atom_id(1);
        let hcb2_id = dummy_atom_id(2);
        let hcb3_id = dummy_atom_id(3);

        residue.add_atom("HCB", hcb1_id);
        residue.add_atom("HCB", hcb2_id);
        residue.add_atom("HCB", hcb3_id);

        assert_eq!(residue.atoms(), &[hcb1_id, hcb2_id, hcb3_id]);
        let retrieved_ids = residue.get_atom_ids_by_name("HCB").unwrap();
        assert_eq!(retrieved_ids, &[hcb1_id, hcb2_id, hcb3_id]);
        assert_eq!(residue.get_first_atom_id_by_name("HCB"), Some(hcb1_id));
    }

    #[test]
    fn residue_removes_one_of_many_duplicates() {
        let mut residue = Residue::new(2, "GLY", Some(ResidueType::Glycine), dummy_chain_id(1));
        let atom_id1 = dummy_atom_id(1);
        let atom_id2 = dummy_atom_id(2);
        residue.add_atom("HCA", atom_id1);
        residue.add_atom("HCA", atom_id2);

        residue.remove_atom("HCA", atom_id1);

        assert_eq!(residue.atoms(), &[atom_id2]);
        let remaining_ids = residue.get_atom_ids_by_name("HCA").unwrap();
        assert_eq!(remaining_ids, &[atom_id2]);
    }

    #[test]
    fn residue_removes_last_duplicate_and_cleans_up_map() {
        let mut residue = Residue::new(2, "GLY", Some(ResidueType::Glycine), dummy_chain_id(1));
        let atom_id1 = dummy_atom_id(1);
        residue.add_atom("CA", atom_id1);

        assert!(residue.atom_name_map.contains_key("CA"));
        residue.remove_atom("CA", atom_id1);

        assert!(residue.atoms().is_empty());
        assert!(residue.get_atom_ids_by_name("CA").is_none());
        assert!(
            !residue.atom_name_map.contains_key("CA"),
            "Map entry should be removed when its vec is empty"
        );
    }

    #[test]
    fn get_atom_ids_by_name_returns_none_for_unknown_atom() {
        let residue = Residue::new(3, "SER", Some(ResidueType::Serine), dummy_chain_id(2));
        assert!(residue.get_atom_ids_by_name("CB").is_none());
    }

    #[test]
    fn get_first_atom_id_by_name_returns_none_for_unknown_atom() {
        let residue = Residue::new(3, "SER", Some(ResidueType::Serine), dummy_chain_id(2));
        assert!(residue.get_first_atom_id_by_name("CB").is_none());
    }

    #[test]
    fn residue_atoms_returns_empty_slice_when_no_atoms() {
        let residue = Residue::new(4, "THR", Some(ResidueType::Threonine), dummy_chain_id(3));
        assert!(residue.atoms().is_empty());
    }

    #[test]
    fn caches_and_retrieves_sidechain_and_backbone_atoms() {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', crate::core::models::chain::ChainType::Protein);
        let residue_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let mut residue = Residue::new(1, "ALA", Some(ResidueType::Alanine), chain_id);

        let mut n_atom = Atom::new("N", residue_id, Default::default());
        n_atom.role = AtomRole::Backbone;
        let mut ca_atom = Atom::new("CA", residue_id, Default::default());
        ca_atom.role = AtomRole::Backbone;
        let mut cb_atom = Atom::new("CB", residue_id, Default::default());
        cb_atom.role = AtomRole::Sidechain;

        let n_id = system.add_atom_to_residue(residue_id, n_atom).unwrap();
        let ca_id = system.add_atom_to_residue(residue_id, ca_atom).unwrap();
        let cb_id = system.add_atom_to_residue(residue_id, cb_atom).unwrap();

        residue.add_atom("N", n_id);
        residue.add_atom("CA", ca_id);
        residue.add_atom("CB", cb_id);

        assert!(residue.backbone_atoms_cache.is_empty());
        assert!(residue.sidechain_atoms_cache.is_empty());

        let backbone_atoms = residue.backbone_atoms(&system);
        assert_eq!(backbone_atoms, &[n_id, ca_id]);

        let sidechain_atoms = residue.sidechain_atoms(&system);
        assert_eq!(sidechain_atoms, &[cb_id]);

        assert!(!residue.backbone_atoms_cache.is_empty());
        assert!(!residue.sidechain_atoms_cache.is_empty());
        assert_eq!(residue.backbone_atoms(&system), &[n_id, ca_id]);
    }

    #[test]
    fn modifying_atom_list_invalidates_caches() {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', crate::core::models::chain::ChainType::Protein);
        let residue_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let mut residue = Residue::new(1, "ALA", Some(ResidueType::Alanine), chain_id);

        let mut n_atom = Atom::new("N", residue_id, Default::default());
        n_atom.role = AtomRole::Backbone;
        let n_id = system.add_atom_to_residue(residue_id, n_atom).unwrap();
        residue.add_atom("N", n_id);

        residue.backbone_atoms(&system);
        assert!(!residue.backbone_atoms_cache.is_empty());

        let mut ca_atom = Atom::new("CA", residue_id, Default::default());
        ca_atom.role = AtomRole::Backbone;
        let ca_id = system.add_atom_to_residue(residue_id, ca_atom).unwrap();
        residue.add_atom("CA", ca_id);

        assert!(residue.backbone_atoms_cache.is_empty());

        let backbone_atoms = residue.backbone_atoms(&system);
        assert_eq!(backbone_atoms, &[n_id, ca_id]);
        assert!(!residue.backbone_atoms_cache.is_empty());

        residue.remove_atom("N", n_id);
        assert!(residue.backbone_atoms_cache.is_empty());

        let backbone_atoms_after_remove = residue.backbone_atoms(&system);
        assert_eq!(backbone_atoms_after_remove, &[ca_id]);
    }
}
