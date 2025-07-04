use super::ids::{AtomId, ChainId};
use std::collections::HashMap;
use std::str::FromStr;
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AminoAcidType {
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Residue {
    pub id: isize,                          // Residue sequence number from source file
    pub name: String,                       // Name of the residue (e.g., "ALA", "GLY")
    pub chain_id: ChainId,                  // ID of the parent chain
    pub(crate) atoms: Vec<AtomId>,          // Indices of atoms belonging to this residue
    atom_name_map: HashMap<String, AtomId>, // Map from atom name to its stable ID
}

impl Residue {
    pub(crate) fn new(id: isize, name: &str, chain_id: ChainId) -> Self {
        Self {
            id,
            name: name.to_string(),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::ids::{AtomId, ChainId};
    use slotmap::KeyData;
    use std::collections::HashSet;

    fn dummy_atom_id(n: u64) -> AtomId {
        AtomId::from(KeyData::from_ffi(n))
    }

    fn dummy_chain_id(n: u64) -> ChainId {
        ChainId::from(KeyData::from_ffi(n))
    }

    #[test]
    fn new_residue_initializes_fields_correctly() {
        let chain_id = dummy_chain_id(1);
        let residue = Residue::new(10, "GLY", chain_id);
        assert_eq!(residue.id, 10);
        assert_eq!(residue.name, "GLY");
        assert_eq!(residue.chain_id, chain_id);
        assert!(residue.atoms().is_empty());
        assert!(residue.get_atom_id_by_name("CA").is_none());
    }

    #[test]
    fn add_atom_adds_atom_and_maps_name() {
        let chain_id = dummy_chain_id(2);
        let mut residue = Residue::new(5, "ALA", chain_id);
        let atom_id = dummy_atom_id(42);
        residue.add_atom("CA", atom_id);
        assert_eq!(residue.atoms(), &[atom_id]);
        assert_eq!(residue.get_atom_id_by_name("CA"), Some(atom_id));
    }

    #[test]
    fn add_atom_allows_multiple_atoms_with_different_names() {
        let chain_id = dummy_chain_id(3);
        let mut residue = Residue::new(7, "SER", chain_id);
        let atom_id1 = dummy_atom_id(1);
        let atom_id2 = dummy_atom_id(2);
        residue.add_atom("CA", atom_id1);
        residue.add_atom("CB", atom_id2);
        let atom_set: HashSet<_> = residue.atoms().iter().copied().collect();
        assert!(atom_set.contains(&atom_id1));
        assert!(atom_set.contains(&atom_id2));
        assert_eq!(residue.get_atom_id_by_name("CA"), Some(atom_id1));
        assert_eq!(residue.get_atom_id_by_name("CB"), Some(atom_id2));
    }

    #[test]
    fn remove_atom_removes_atom_and_name_mapping() {
        let chain_id = dummy_chain_id(4);
        let mut residue = Residue::new(8, "THR", chain_id);
        let atom_id = dummy_atom_id(100);
        residue.add_atom("OG1", atom_id);
        residue.remove_atom("OG1", atom_id);
        assert!(residue.atoms().is_empty());
        assert!(residue.get_atom_id_by_name("OG1").is_none());
    }

    #[test]
    fn remove_atom_does_nothing_if_atom_not_present() {
        let chain_id = dummy_chain_id(5);
        let mut residue = Residue::new(9, "VAL", chain_id);
        let atom_id = dummy_atom_id(200);
        residue.add_atom("CG1", atom_id);
        residue.remove_atom("CG2", dummy_atom_id(201));
        assert_eq!(residue.atoms(), &[atom_id]);
        assert_eq!(residue.get_atom_id_by_name("CG1"), Some(atom_id));
    }

    #[test]
    fn get_atom_id_by_name_returns_none_for_unknown_name() {
        let chain_id = dummy_chain_id(6);
        let mut residue = Residue::new(11, "LEU", chain_id);
        residue.add_atom("CD1", dummy_atom_id(300));
        assert!(residue.get_atom_id_by_name("CD2").is_none());
    }
}
