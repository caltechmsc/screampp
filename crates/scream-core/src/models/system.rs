use super::atom::{Atom, Element};
use super::chain::{Chain, ChainType};
use super::residue::Residue;
use super::topology::{Bond, BondOrder};
use nalgebra::Point3;
use std::collections::HashMap;

pub struct MolecularSystem {
    atoms: Vec<Atom>,
    chains: Vec<Chain>,
    bonds: Vec<Bond>,
    atom_serial_map: HashMap<usize, usize>,
    chain_id_map: HashMap<char, usize>,
}

impl MolecularSystem {
    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }
    pub fn chains(&self) -> &[Chain] {
        &self.chains
    }
    pub fn bonds(&self) -> &[Bond] {
        &self.bonds
    }

    pub fn get_atom(&self, index: usize) -> Option<&Atom> {
        self.atoms.get(index)
    }

    pub fn get_atom_by_serial(&self, serial: usize) -> Option<&Atom> {
        self.atom_serial_map
            .get(&serial)
            .and_then(|&idx| self.atoms.get(idx))
    }

    pub fn get_chain(&self, index: usize) -> Option<&Chain> {
        self.chains.get(index)
    }

    pub fn get_chain_by_id(&self, id: char) -> Option<&Chain> {
        self.chain_id_map
            .get(&id)
            .and_then(|&idx| self.chains.get(idx))
    }

    pub fn atoms_in_residue<'a>(&'a self, residue: &'a Residue) -> impl Iterator<Item = &'a Atom> {
        residue
            .atom_indices
            .iter()
            .map(move |&idx| &self.atoms[idx])
    }
}

pub struct MolecularSystemBuilder {
    system: MolecularSystem,
    current_chain_idx: Option<usize>,
    current_residue_idx: Option<usize>,
}

impl Default for MolecularSystemBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl MolecularSystemBuilder {
    pub fn new() -> Self {
        Self {
            system: MolecularSystem {
                atoms: Vec::new(),
                chains: Vec::new(),
                bonds: Vec::new(),
                atom_serial_map: HashMap::new(),
                chain_id_map: HashMap::new(),
            },
            current_chain_idx: None,
            current_residue_idx: None,
        }
    }

    pub fn start_chain(&mut self, id: char, chain_type: ChainType) -> &mut Self {
        let idx = *self.system.chain_id_map.entry(id).or_insert_with(|| {
            let index = self.system.chains.len();
            self.system.chains.push(Chain::new(id, chain_type));
            index
        });
        self.current_chain_idx = Some(idx);
        self.current_residue_idx = None;
        self
    }

    pub fn start_residue(&mut self, id: isize, name: &str) -> &mut Self {
        let chain_idx = self
            .current_chain_idx
            .expect("Must start a chain before a residue");
        let chain = &mut self.system.chains[chain_idx];

        let res_idx = *chain.residue_map.entry(id).or_insert_with(|| {
            let index = chain.residues.len();
            chain.residues.push(Residue::new(id, name));
            index
        });
        self.current_residue_idx = Some(res_idx);
        self
    }

    pub fn add_atom(
        &mut self,
        serial: usize,
        name: &str,
        element: Element,
        position: Point3<f64>,
        charge: f64,
        ff_type: &str,
    ) -> &mut Self {
        let chain_idx = self
            .current_chain_idx
            .expect("Cannot add atom without a current chain");
        let res_idx = self
            .current_residue_idx
            .expect("Cannot add atom without a current residue");

        let atom_idx = self.system.atoms.len();
        let atom = Atom {
            index: atom_idx,
            serial,
            name: name.to_string(),
            element,
            position,
            partial_charge: charge,
            force_field_type: ff_type.to_string(),
        };

        self.system.atoms.push(atom);
        self.system.atom_serial_map.insert(serial, atom_idx);
        self.system.chains[chain_idx].residues[res_idx].add_atom(name, atom_idx);
        self
    }

    pub fn add_bond(&mut self, serial1: usize, serial2: usize, order: BondOrder) -> &mut Self {
        let idx1 = *self
            .system
            .atom_serial_map
            .get(&serial1)
            .expect("Atom 1 not found for bond");
        let idx2 = *self
            .system
            .atom_serial_map
            .get(&serial2)
            .expect("Atom 2 not found for bond");
        self.system.bonds.push(Bond::new(idx1, idx2, order));
        self
    }

    pub fn build(self) -> MolecularSystem {
        self.system
    }
}
