use super::atom::{Atom, AtomFlags};
use super::chain::{Chain, ChainType};
use super::ids::{AtomId, ChainId, ResidueId};
use super::residue::Residue;
use super::topology::{Bond, BondOrder};
use nalgebra::Point3;
use slotmap::{SecondaryMap, SlotMap};
use std::collections::HashMap;

#[derive(Debug, Clone, Default)]
pub struct MolecularSystem {
    // --- Primary Data Stores (Source of Truth) ---
    atoms: SlotMap<AtomId, Atom>,
    residues: SlotMap<ResidueId, Residue>,
    chains: SlotMap<ChainId, Chain>,
    bonds: Vec<Bond>,

    // --- Lookup Maps for Fast Access by External Identifiers ---
    atom_serial_map: HashMap<usize, AtomId>,
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

    pub fn find_atom_by_serial(&self, serial: usize) -> Option<AtomId> {
        self.atom_serial_map.get(&serial).copied()
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
    ) -> Option<ResidueId> {
        let chain = self.chains.get_mut(chain_id)?;
        let key = (chain_id, res_seq);
        let residue_id = *self.residue_id_map.entry(key).or_insert_with(|| {
            let residue = Residue::new(res_seq, name, chain_id);
            self.residues.insert(residue)
        });

        let chain_residues = &mut self.chains.get_mut(chain_id).unwrap().residues;
        if !chain_residues.contains(&residue_id) {
            chain_residues.push(residue_id);
        }

        Some(residue_id)
    }

    pub fn add_atom_to_residue(&mut self, residue_id: ResidueId, atom: Atom) -> Option<AtomId> {
        if !self.residues.contains_key(residue_id) {
            return None;
        }

        let serial = atom.serial;
        let name = atom.name.clone();

        let atom_id = self.atoms.insert(atom);
        self.bond_adjacency.insert(atom_id, Vec::new());
        self.atom_serial_map.insert(serial, atom_id);

        let residue = self.residues.get_mut(residue_id).unwrap();
        residue.add_atom(&name, atom_id);

        Some(atom_id)
    }

    pub fn add_bond(&mut self, atom1_id: AtomId, atom2_id: AtomId, order: BondOrder) -> Option<()> {
        if self.atoms.contains_key(atom1_id) && self.atoms.contains_key(atom2_id) {
            self.bonds.push(Bond::new(atom1_id, atom2_id, order));
            self.bond_adjacency[atom1_id].push(atom2_id);
            self.bond_adjacency[atom2_id].push(atom1_id);
            Some(())
        } else {
            None
        }
    }

    pub fn remove_atom(&mut self, atom_id: AtomId) -> Option<Atom> {
        let atom = self.atoms.remove(atom_id)?;

        // 1. Remove from parent residue
        let residue = self.residues.get_mut(atom.residue_id).unwrap();
        residue.remove_atom(&atom.name, atom_id);

        // 2. Remove from serial map
        self.atom_serial_map.remove(&atom.serial);

        // 3. Remove all bonds connected to this atom
        let original_bonds = std::mem::take(&mut self.bonds);
        self.bonds = original_bonds
            .into_iter()
            .filter(|bond| !bond.contains(atom_id))
            .collect();

        // 4. Clean up adjacency list
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
        self.residue_id_map.remove(&(residue.chain_id, residue.id));

        // 4. Finally, remove the residue itself
        self.residues.remove(residue_id)
    }
}
