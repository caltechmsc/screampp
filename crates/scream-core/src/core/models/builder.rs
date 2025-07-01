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

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    #[test]
    fn builder_creates_system_with_single_chain_residue_atom() {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('A', ChainType::Protein)
            .start_residue(1, "GLY")
            .add_atom(
                100,
                "CA",
                "GLY",
                Point3::new(1.0, 2.0, 3.0),
                Some(-0.1),
                Some("C.3"),
            );
        let system = builder.build();
        assert_eq!(system.chains().len(), 1);
        assert_eq!(system.chains()[0].residues().len(), 1);
        assert_eq!(system.atoms().len(), 1);
        assert_eq!(system.atoms()[0].serial, 100);
        assert_eq!(system.atoms()[0].name, "CA");
        assert_eq!(system.atoms()[0].res_name, "GLY");
        assert_eq!(system.atoms()[0].partial_charge, -0.1);
        assert_eq!(system.atoms()[0].force_field_type, "C.3");
    }

    #[test]
    fn builder_adds_multiple_chains_and_residues() {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('A', ChainType::Protein)
            .start_residue(1, "GLY")
            .add_atom(1, "CA", "GLY", Point3::new(0.0, 0.0, 0.0), None, None)
            .start_residue(2, "ALA")
            .add_atom(2, "CB", "ALA", Point3::new(1.0, 1.0, 1.0), None, None)
            .start_chain('B', ChainType::Protein)
            .start_residue(1, "SER")
            .add_atom(3, "OG", "SER", Point3::new(2.0, 2.0, 2.0), None, None);
        let system = builder.build();
        assert_eq!(system.chains().len(), 2);
        assert_eq!(system.chains()[0].residues().len(), 2);
        assert_eq!(system.chains()[1].residues().len(), 1);
        assert_eq!(system.atoms().len(), 3);
    }

    #[test]
    fn builder_adds_bond_between_atoms() {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('A', ChainType::Protein)
            .start_residue(1, "GLY")
            .add_atom(1, "CA", "GLY", Point3::new(0.0, 0.0, 0.0), None, None)
            .add_atom(2, "CB", "GLY", Point3::new(1.0, 1.0, 1.0), None, None)
            .add_bond(1, 2, BondOrder::Single);
        let system = builder.build();
        assert_eq!(system.bonds().len(), 1);
        let bond = &system.bonds()[0];
        assert_eq!(bond.atom1_idx, 0);
        assert_eq!(bond.atom2_idx, 1);
    }

    #[test]
    fn builder_add_atom_without_charge_and_ff_type_defaults() {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('A', ChainType::Protein)
            .start_residue(1, "GLY")
            .add_atom(1, "CA", "GLY", Point3::new(0.0, 0.0, 0.0), None, None);
        let system = builder.build();
        let atom = &system.atoms()[0];
        assert_eq!(atom.partial_charge, 0.0);
        assert_eq!(atom.force_field_type, "");
    }

    #[test]
    #[should_panic(expected = "Must start a chain before starting a residue")]
    fn builder_panics_if_start_residue_without_chain() {
        let mut builder = MolecularSystemBuilder::new();
        builder.start_residue(1, "GLY");
    }

    #[test]
    #[should_panic(expected = "Cannot add atom without a current chain")]
    fn builder_panics_if_add_atom_without_chain() {
        let mut builder = MolecularSystemBuilder::new();
        builder.add_atom(1, "CA", "GLY", Point3::new(0.0, 0.0, 0.0), None, None);
    }

    #[test]
    #[should_panic(expected = "Cannot add atom without a current residue")]
    fn builder_panics_if_add_atom_without_residue() {
        let mut builder = MolecularSystemBuilder::new();
        builder.start_chain('A', ChainType::Protein);
        builder.add_atom(1, "CA", "GLY", Point3::new(0.0, 0.0, 0.0), None, None);
    }

    #[test]
    #[should_panic(expected = "Atom 1 for bond not found")]
    fn builder_panics_if_add_bond_with_invalid_serial1() {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('A', ChainType::Protein)
            .start_residue(1, "GLY")
            .add_atom(1, "CA", "GLY", Point3::new(0.0, 0.0, 0.0), None, None)
            .add_bond(99, 1, BondOrder::Single);
    }

    #[test]
    #[should_panic(expected = "Atom 2 for bond not found")]
    fn builder_panics_if_add_bond_with_invalid_serial2() {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('A', ChainType::Protein)
            .start_residue(1, "GLY")
            .add_atom(1, "CA", "GLY", Point3::new(0.0, 0.0, 0.0), None, None)
            .add_bond(1, 99, BondOrder::Single);
    }
}
