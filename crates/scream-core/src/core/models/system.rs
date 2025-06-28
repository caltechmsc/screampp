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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::atom::{Atom, Element};
    use crate::core::models::chain::ChainType;
    use nalgebra::Point3;

    #[test]
    fn builder_creates_system_with_single_chain_residue_and_atom() {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('A', ChainType::Protein)
            .start_residue(1, "GLY")
            .add_atom(
                100,
                "CA",
                Element::C,
                Point3::new(1.0, 2.0, 3.0),
                0.0,
                "C.3",
            );
        let system = builder.build();

        assert_eq!(system.chains().len(), 1);
        assert_eq!(system.atoms().len(), 1);
        assert_eq!(system.bonds().len(), 0);
        assert_eq!(system.get_chain_by_id('A').unwrap().residues().len(), 1);
        assert_eq!(system.get_atom_by_serial(100).unwrap().name, "CA");
    }

    #[test]
    fn builder_adds_multiple_chains_and_residues() {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('A', ChainType::Protein)
            .start_residue(1, "GLY")
            .add_atom(1, "CA", Element::C, Point3::new(0.0, 0.0, 0.0), 0.0, "C.3")
            .start_residue(2, "ALA")
            .add_atom(2, "CB", Element::C, Point3::new(1.0, 0.0, 0.0), 0.0, "C.3")
            .start_chain('B', ChainType::DNA)
            .start_residue(1, "DA")
            .add_atom(3, "N1", Element::N, Point3::new(0.0, 1.0, 0.0), 0.0, "N.1");
        let system = builder.build();

        assert_eq!(system.chains().len(), 2);
        assert_eq!(system.chains()[0].residues().len(), 2);
        assert_eq!(system.chains()[1].residues().len(), 1);
        assert_eq!(system.atoms().len(), 3);
    }

    #[test]
    fn builder_adds_bond_between_atoms_by_serial() {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('A', ChainType::Protein)
            .start_residue(1, "GLY")
            .add_atom(10, "CA", Element::C, Point3::new(0.0, 0.0, 0.0), 0.0, "C.3")
            .add_atom(11, "CB", Element::C, Point3::new(1.0, 0.0, 0.0), 0.0, "C.3")
            .add_bond(10, 11, BondOrder::Single);
        let system = builder.build();

        assert_eq!(system.bonds().len(), 1);
        let bond = &system.bonds()[0];
        assert_eq!(bond.order, BondOrder::Single);
        assert!(bond.atom1_idx != bond.atom2_idx);
    }

    #[test]
    fn get_atom_and_chain_by_index_and_id() {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('C', ChainType::Other)
            .start_residue(5, "UNK")
            .add_atom(
                42,
                "X",
                Element::Unknown,
                Point3::new(0.0, 0.0, 0.0),
                0.0,
                "X",
            );
        let system = builder.build();

        assert!(system.get_atom(0).is_some());
        assert!(system.get_chain(0).is_some());
        assert!(system.get_chain_by_id('C').is_some());
        assert!(system.get_atom_by_serial(42).is_some());
        assert!(system.get_atom(1).is_none());
        assert!(system.get_chain(1).is_none());
        assert!(system.get_chain_by_id('Z').is_none());
        assert!(system.get_atom_by_serial(999).is_none());
    }

    #[test]
    fn atoms_in_residue_returns_correct_atoms() {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('A', ChainType::Protein)
            .start_residue(1, "GLY")
            .add_atom(1, "CA", Element::C, Point3::new(0.0, 0.0, 0.0), 0.0, "C.3")
            .add_atom(2, "CB", Element::C, Point3::new(1.0, 0.0, 0.0), 0.0, "C.3");
        let system = builder.build();
        let residue = &system.chains()[0].residues()[0];
        let atom_names: Vec<_> = system
            .atoms_in_residue(residue)
            .map(|a| a.name.as_str())
            .collect();

        assert_eq!(atom_names, vec!["CA", "CB"]);
    }

    #[test]
    #[should_panic(expected = "Must start a chain before a residue")]
    fn start_residue_without_chain_panics() {
        let mut builder = MolecularSystemBuilder::new();
        builder.start_residue(1, "GLY");
    }

    #[test]
    #[should_panic(expected = "Cannot add atom without a current chain")]
    fn add_atom_without_chain_panics() {
        let mut builder = MolecularSystemBuilder::new();
        builder.add_atom(1, "CA", Element::C, Point3::new(0.0, 0.0, 0.0), 0.0, "C.3");
    }

    #[test]
    #[should_panic(expected = "Cannot add atom without a current residue")]
    fn add_atom_without_residue_panics() {
        let mut builder = MolecularSystemBuilder::new();
        builder.start_chain('A', ChainType::Protein);
        builder.add_atom(1, "CA", Element::C, Point3::new(0.0, 0.0, 0.0), 0.0, "C.3");
    }

    #[test]
    #[should_panic(expected = "Atom 1 not found for bond")]
    fn add_bond_with_invalid_first_serial_panics() {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('A', ChainType::Protein)
            .start_residue(1, "GLY")
            .add_atom(1, "CA", Element::C, Point3::new(0.0, 0.0, 0.0), 0.0, "C.3");
        builder.add_bond(999, 1, BondOrder::Single);
    }

    #[test]
    #[should_panic(expected = "Atom 2 not found for bond")]
    fn add_bond_with_invalid_second_serial_panics() {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('A', ChainType::Protein)
            .start_residue(1, "GLY")
            .add_atom(1, "CA", Element::C, Point3::new(0.0, 0.0, 0.0), 0.0, "C.3");
        builder.add_bond(1, 999, BondOrder::Single);
    }
}
