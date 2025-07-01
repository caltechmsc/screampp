use super::atom::{Atom, AtomFlags};
use super::chain::{Chain, ChainType};
use super::residue::Residue;
use super::system::MolecularSystem;
use super::topology::{Bond, BondOrder};
use nalgebra::Point3;
use std::collections::HashMap;

pub struct MolecularSystemBuilder {
    system: MolecularSystem,

    // --- Builder-specific state for efficient construction ---
    atom_serial_map: HashMap<usize, usize>,
    chain_id_map: HashMap<char, usize>,
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
            },
            atom_serial_map: HashMap::new(),
            chain_id_map: HashMap::new(),
            current_chain_idx: None,
            current_residue_idx: None,
        }
    }

    pub fn start_chain(&mut self, id: char, chain_type: ChainType) -> &mut Self {
        let idx = *self.chain_id_map.entry(id).or_insert_with(|| {
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
            .expect("Must start a chain before starting a residue");
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
        res_name: &str,
        position: Point3<f64>,
        charge: Option<f64>,
        ff_type: Option<&str>,
    ) -> &mut Self {
        let chain_idx = self
            .current_chain_idx
            .expect("Cannot add atom without a current chain");
        let res_idx = self
            .current_residue_idx
            .expect("Cannot add atom without a current residue");

        let chain_id = self.system.chains[chain_idx].id;
        let res_id = self.system.chains[chain_idx].residues[res_idx].id;
        let atom_idx = self.system.atoms.len();

        let atom = Atom {
            index: atom_idx,
            serial,
            name: name.to_string(),
            res_name: res_name.to_string(),
            res_id,
            chain_id,
            position,
            partial_charge: charge.unwrap_or(0.0),
            force_field_type: ff_type.unwrap_or("").to_string(),
            flags: AtomFlags::default(),
            delta: 0.0,
            vdw_radius: 0.0,
            vdw_well_depth: 0.0,
            hbond_type_id: -1,
        };

        self.system.atoms.push(atom);
        self.atom_serial_map.insert(serial, atom_idx);
        self.system.chains[chain_idx].residues[res_idx].add_atom(name, atom_idx);
        self
    }

    pub fn add_bond(&mut self, serial1: usize, serial2: usize, order: BondOrder) -> &mut Self {
        let idx1 = *self
            .atom_serial_map
            .get(&serial1)
            .expect("Atom 1 for bond not found");
        let idx2 = *self
            .atom_serial_map
            .get(&serial2)
            .expect("Atom 2 for bond not found");
        self.system.bonds.push(Bond::new(idx1, idx2, order));
        self
    }

    pub fn build(self) -> MolecularSystem {
        self.system
    }
}
