use super::atom::Atom;
use super::chain::{Chain, ChainType};
use super::ids::{AtomId, ChainId, ResidueId};
use super::residue::{Residue, ResidueType};
use super::topology::{Bond, BondOrder};
use slotmap::{SecondaryMap, SlotMap};
use std::collections::HashMap;

#[derive(Debug, Clone, Default)]
pub struct MolecularSystem {
    // --- Primary Data Stores (Source of Truth) ---
    atoms: SlotMap<AtomId, Atom>,
    residues: SlotMap<ResidueId, Residue>,
    chains: SlotMap<ChainId, Chain>,
    bonds: Vec<Bond>,

    // --- Lookup Maps for Fast Access ---
    residue_id_map: HashMap<(ChainId, isize), ResidueId>,
    chain_id_map: HashMap<char, ChainId>,

    // --- Adjacency information (cache for performance) ---
    bond_adjacency: SecondaryMap<AtomId, Vec<AtomId>>,
}

impl MolecularSystem {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn atom(&self, id: AtomId) -> Option<&Atom> {
        self.atoms.get(id)
    }
    pub fn atom_mut(&mut self, id: AtomId) -> Option<&mut Atom> {
        self.atoms.get_mut(id)
    }
    pub fn atoms_iter(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.atoms.iter()
    }

    pub fn residue(&self, id: ResidueId) -> Option<&Residue> {
        self.residues.get(id)
    }
    pub fn residue_mut(&mut self, id: ResidueId) -> Option<&mut Residue> {
        self.residues.get_mut(id)
    }
    pub fn residues_iter(&self) -> impl Iterator<Item = (ResidueId, &Residue)> {
        self.residues.iter()
    }

    pub fn chain(&self, id: ChainId) -> Option<&Chain> {
        self.chains.get(id)
    }
    pub fn chain_mut(&mut self, id: ChainId) -> Option<&mut Chain> {
        self.chains.get_mut(id)
    }
    pub fn chains_iter(&self) -> impl Iterator<Item = (ChainId, &Chain)> {
        self.chains.iter()
    }

    pub fn bonds(&self) -> &[Bond] {
        &self.bonds
    }

    pub fn find_chain_by_id(&self, id: char) -> Option<ChainId> {
        self.chain_id_map.get(&id).copied()
    }
    pub fn find_residue_by_id(&self, chain_id: ChainId, res_seq: isize) -> Option<ResidueId> {
        self.residue_id_map.get(&(chain_id, res_seq)).copied()
    }

    pub fn add_chain(&mut self, id: char, chain_type: ChainType) -> ChainId {
        *self.chain_id_map.entry(id).or_insert_with(|| {
            let chain = Chain::new(id, chain_type);
            self.chains.insert(chain)
        })
    }

    pub fn add_residue(
        &mut self,
        chain_id: ChainId,
        res_seq: isize,
        name: &str,
        residue_type: Option<ResidueType>,
    ) -> Option<ResidueId> {
        let chain = self.chains.get_mut(chain_id)?;
        let key = (chain_id, res_seq);

        let residue_id = *self.residue_id_map.entry(key).or_insert_with(|| {
            let residue = Residue::new(res_seq, name, residue_type, chain_id);
            self.residues.insert(residue)
        });

        if !chain.residues.contains(&residue_id) {
            chain.residues.push(residue_id);
        }

        Some(residue_id)
    }

    pub fn add_atom_to_residue(&mut self, residue_id: ResidueId, atom: Atom) -> Option<AtomId> {
        if !self.residues.contains_key(residue_id) {
            return None;
        }

        let name = atom.name.clone();

        let atom_id = self.atoms.insert(atom);
        self.bond_adjacency.insert(atom_id, Vec::new());

        let residue = self.residues.get_mut(residue_id).unwrap();
        residue.add_atom(&name, atom_id);

        Some(atom_id)
    }

    pub fn add_bond(&mut self, atom1_id: AtomId, atom2_id: AtomId, order: BondOrder) -> Option<()> {
        if !self.atoms.contains_key(atom1_id) || !self.atoms.contains_key(atom2_id) {
            return None;
        }

        if let Some(neighbors) = self.bond_adjacency.get(atom1_id) {
            if neighbors.contains(&atom2_id) {
                // Bond already exists, operation is successful (idempotent)
                return Some(());
            }
        }

        self.bonds.push(Bond::new(atom1_id, atom2_id, order));
        self.bond_adjacency[atom1_id].push(atom2_id);
        self.bond_adjacency[atom2_id].push(atom1_id);
        Some(())
    }

    pub fn remove_atom(&mut self, atom_id: AtomId) -> Option<Atom> {
        let atom = self.atoms.remove(atom_id)?;

        // 1. Remove from parent residue
        if let Some(residue) = self.residues.get_mut(atom.residue_id) {
            residue.remove_atom(&atom.name, atom_id);
        }

        // 2. Remove all bonds connected to this atom
        let original_bonds = std::mem::take(&mut self.bonds);
        self.bonds = original_bonds
            .into_iter()
            .filter(|bond| !bond.contains(atom_id))
            .collect();

        // 3. Clean up adjacency list
        let neighbors = self.bond_adjacency.remove(atom_id).unwrap_or_default();
        for neighbor_id in neighbors {
            if let Some(adjacency) = self.bond_adjacency.get_mut(neighbor_id) {
                adjacency.retain(|&id| id != atom_id);
            }
        }

        Some(atom)
    }

    pub fn remove_residue(&mut self, residue_id: ResidueId) -> Option<Residue> {
        let residue = self.residues.get(residue_id)?.clone(); // Clone to avoid borrow checker issues

        // 1. Remove all atoms within the residue
        for atom_id in residue.atoms().to_vec() {
            // to_vec to avoid borrow issues
            self.remove_atom(atom_id);
        }

        // 2. Remove from parent chain
        if let Some(chain) = self.chains.get_mut(residue.chain_id) {
            chain.residues.retain(|&id| id != residue_id);
        }

        // 3. Remove from residue maps
        self.residue_id_map
            .remove(&(residue.chain_id, residue.residue_number));

        // 4. Finally, remove the residue itself
        self.residues.remove(residue_id)
    }

    pub fn get_bonded_neighbors(&self, atom_id: AtomId) -> Option<&[AtomId]> {
        self.bond_adjacency.get(atom_id).map(|v| v.as_slice())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::residue::ResidueType;
    use nalgebra::Point3;
    use std::collections::HashMap;

    fn create_test_system() -> (MolecularSystem, HashMap<String, AtomId>) {
        let mut system = MolecularSystem::new();
        let mut key_atom_ids = HashMap::new();

        let chain_a_id = system.add_chain('A', ChainType::Protein);

        let gly_id = system
            .add_residue(chain_a_id, 1, "GLY", Some(ResidueType::Glycine))
            .unwrap();
        let n_atom = Atom::new("N", gly_id, Point3::new(0.0, 0.0, 0.0));
        let ca_atom = Atom::new("CA", gly_id, Point3::new(1.4, 0.0, 0.0));

        let n_id = system.add_atom_to_residue(gly_id, n_atom).unwrap();
        key_atom_ids.insert("GLY_1_N".to_string(), n_id);

        let ca_id = system.add_atom_to_residue(gly_id, ca_atom).unwrap();
        key_atom_ids.insert("GLY_1_CA".to_string(), ca_id);
        system.add_bond(n_id, ca_id, BondOrder::Single).unwrap();

        let ala_id = system
            .add_residue(chain_a_id, 2, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let ala_ca_atom = Atom::new("CA", ala_id, Point3::new(2.0, 1.0, 0.0));

        let ala_ca_id = system.add_atom_to_residue(ala_id, ala_ca_atom).unwrap();
        key_atom_ids.insert("ALA_2_CA".to_string(), ala_ca_id);

        system
            .add_bond(ca_id, ala_ca_id, BondOrder::Single)
            .unwrap();

        (system, key_atom_ids)
    }

    #[test]
    fn system_creation_and_access() {
        let (system, key_atom_ids) = create_test_system();

        assert_eq!(system.atoms_iter().count(), 3);
        assert_eq!(system.residues_iter().count(), 2);
        assert_eq!(system.chains_iter().count(), 1);
        assert_eq!(system.bonds().len(), 2);

        let chain_a_id = system.find_chain_by_id('A').unwrap();
        assert!(system.find_chain_by_id('B').is_none());

        let gly_id = system.find_residue_by_id(chain_a_id, 1).unwrap();
        let ala_id = system.find_residue_by_id(chain_a_id, 2).unwrap();

        assert_eq!(system.residue(gly_id).unwrap().name, "GLY");
        assert_eq!(system.residue(ala_id).unwrap().name, "ALA");

        let atom_n_id = *key_atom_ids.get("GLY_1_N").unwrap();
        assert_eq!(system.atom(atom_n_id).unwrap().name, "N");
    }

    #[test]
    fn atom_removal_updates_system_correctly() {
        let (mut system, key_atom_ids) = create_test_system();
        let atom_n_id = *key_atom_ids.get("GLY_1_N").unwrap();
        let atom_ca_id = *key_atom_ids.get("GLY_1_CA").unwrap();

        assert_eq!(system.bonds().len(), 2);
        assert!(system.bonds().iter().any(|b| b.contains(atom_n_id)));
        assert!(system.atom(atom_n_id).is_some());
        let residue_id = system.atom(atom_n_id).unwrap().residue_id;
        assert_eq!(system.residue(residue_id).unwrap().residue_number, 1);

        let removed_atom = system.remove_atom(atom_n_id).unwrap();
        assert_eq!(removed_atom.name, "N");

        assert_eq!(system.atoms_iter().count(), 2);
        assert!(system.atom(atom_n_id).is_none());

        assert_eq!(system.bonds().len(), 1);

        assert!(system.bond_adjacency.get(atom_n_id).is_none());
        assert!(!system.bond_adjacency[atom_ca_id].contains(&atom_n_id));

        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let gly_id = system.find_residue_by_id(chain_a_id, 1).unwrap();
        assert_eq!(system.residue(gly_id).unwrap().atoms().len(), 1);
        assert!(!system.residue(gly_id).unwrap().atoms().contains(&atom_n_id));
    }

    #[test]
    fn residue_removal_updates_system_correctly() {
        let (mut system, key_atom_ids) = create_test_system();
        let atom_ala_ca_id = *key_atom_ids.get("ALA_2_CA").unwrap();

        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let gly_id = system.find_residue_by_id(chain_a_id, 1).unwrap();

        assert_eq!(system.atoms_iter().count(), 3);
        assert_eq!(system.residues_iter().count(), 2);
        assert_eq!(system.bonds().len(), 2);
        assert_eq!(system.chain(chain_a_id).unwrap().residues().len(), 2);

        let removed_residue = system.remove_residue(gly_id).unwrap();
        assert_eq!(removed_residue.name, "GLY");

        assert_eq!(system.residues_iter().count(), 1);
        assert!(system.residue(gly_id).is_none());
        assert!(system.find_residue_by_id(chain_a_id, 1).is_none());

        assert_eq!(system.atoms_iter().count(), 1);
        assert!(system.atom(*key_atom_ids.get("GLY_1_N").unwrap()).is_none());
        assert!(
            system
                .atom(*key_atom_ids.get("GLY_1_CA").unwrap())
                .is_none()
        );
        assert!(system.atom(atom_ala_ca_id).is_some());

        assert!(system.bonds().is_empty());

        assert_eq!(system.chain(chain_a_id).unwrap().residues().len(), 1);
        assert!(
            !system
                .chain(chain_a_id)
                .unwrap()
                .residues()
                .contains(&gly_id)
        );
    }

    #[test]
    fn system_creation_and_access_with_residue_type() {
        let (system, _) = create_test_system();
        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let gly_id = system.find_residue_by_id(chain_a_id, 1).unwrap();
        let ala_id = system.find_residue_by_id(chain_a_id, 2).unwrap();

        assert_eq!(system.residue(gly_id).unwrap().name, "GLY");
        assert_eq!(
            system.residue(gly_id).unwrap().residue_type,
            Some(ResidueType::Glycine)
        );
        assert_eq!(
            system.residue(ala_id).unwrap().residue_type,
            Some(ResidueType::Alanine)
        );
    }

    #[test]
    fn get_bonded_neighbors_returns_correct_neighbors() {
        let (system, _) = create_test_system();

        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let gly_residue = system
            .residue(system.find_residue_by_id(chain_a_id, 1).unwrap())
            .unwrap();
        let ala_residue = system
            .residue(system.find_residue_by_id(chain_a_id, 2).unwrap())
            .unwrap();

        let atom_n_id = gly_residue.get_atom_id_by_name("N").unwrap();
        let atom_gly_ca_id = gly_residue.get_atom_id_by_name("CA").unwrap();
        let atom_ala_ca_id = ala_residue.get_atom_id_by_name("CA").unwrap();

        let neighbors_n = system.get_bonded_neighbors(atom_n_id).unwrap();
        assert_eq!(neighbors_n, &[atom_gly_ca_id]);

        let neighbors_ca = system.get_bonded_neighbors(atom_gly_ca_id).unwrap();
        assert_eq!(neighbors_ca, &[atom_n_id, atom_ala_ca_id]);

        let neighbors_ala_ca = system.get_bonded_neighbors(atom_ala_ca_id).unwrap();
        assert_eq!(neighbors_ala_ca, &[atom_gly_ca_id]);
    }
}
