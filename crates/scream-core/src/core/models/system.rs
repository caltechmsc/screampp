use super::atom::Atom;
use super::chain::Chain;
use super::residue::Residue;
use super::topology::Bond;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct MolecularSystem {
    pub(crate) atoms: Vec<Atom>,
    pub(crate) chains: Vec<Chain>,
    pub(crate) bonds: Vec<Bond>,

    // --- Internal maps for efficient lookups ---
    atom_serial_map: HashMap<usize, usize>, // Maps original serial to internal index
    chain_id_map: HashMap<char, usize>,     // Maps chain ID to internal index
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

    pub fn bonds_mut(&mut self) -> &mut Vec<Bond> {
        &mut self.bonds
    }

    pub fn get_atom(&self, index: usize) -> Option<&Atom> {
        self.atoms.get(index)
    }

    pub fn get_atom_mut(&mut self, index: usize) -> Option<&mut Atom> {
        self.atoms.get_mut(index)
    }

    pub fn get_atom_by_serial(&self, serial: usize) -> Option<&Atom> {
        self.atom_serial_map
            .get(&serial)
            .and_then(|&index| self.atoms.get(index))
    }

    pub fn get_atom_by_serial_mut(&mut self, serial: usize) -> Option<&mut Atom> {
        if let Some(&index) = self.atom_serial_map.get(&serial) {
            self.atoms.get_mut(index)
        } else {
            None
        }
    }

    pub fn get_chain(&self, index: usize) -> Option<&Chain> {
        self.chains.get(index)
    }

    pub fn get_chain_by_id(&self, id: char) -> Option<&Chain> {
        self.chain_id_map
            .get(&id)
            .and_then(|&index| self.chains.get(index))
    }

    pub fn atoms_in_residue<'a>(&'a self, residue: &'a Residue) -> impl Iterator<Item = &'a Atom> {
        residue
            .atom_indices
            .iter()
            .map(move |&idx| &self.atoms[idx])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::atom::{Atom, AtomFlags};
    use crate::core::models::chain::Chain;
    use crate::core::models::residue::Residue;
    use crate::core::models::topology::Bond;
    use nalgebra::Point3;

    fn mock_system() -> MolecularSystem {
        let atom1 = Atom {
            index: 0,
            serial: 10,
            name: "CA".to_string(),
            res_name: "GLY".to_string(),
            res_id: 5,
            chain_id: 'A',
            force_field_type: "C.3".to_string(),
            partial_charge: -0.123,
            position: Point3::new(1.0, 2.0, 3.0),
            flags: AtomFlags::IS_FIXED_ROLE | AtomFlags::IS_VISIBLE_LATTICE,
            delta: 0.5,
            vdw_radius: 1.7,
            vdw_well_depth: 0.2,
            hbond_type_id: 3,
        };
        let atom2 = Atom {
            index: 1,
            serial: 20,
            name: "CB".to_string(),
            res_name: "GLY".to_string(),
            res_id: 5,
            chain_id: 'A',
            force_field_type: "C.3".to_string(),
            partial_charge: -0.456,
            position: Point3::new(2.0, 3.0, 4.0),
            flags: AtomFlags::IS_VISIBLE_INTERACTION,
            delta: 0.6,
            vdw_radius: 1.8,
            vdw_well_depth: 0.3,
            hbond_type_id: 2,
        };
        let residue = Residue::new(5, "GLY");
        let chain = Chain::new('A', crate::core::models::chain::ChainType::Protein);
        let bond = Bond::new(0, 1, crate::core::models::topology::BondOrder::Single);

        let mut atom_serial_map = HashMap::new();
        atom_serial_map.insert(10, 0);
        atom_serial_map.insert(20, 1);

        let mut chain_id_map = HashMap::new();
        chain_id_map.insert('A', 0);

        MolecularSystem {
            atoms: vec![atom1, atom2],
            chains: vec![chain],
            bonds: vec![bond],
            atom_serial_map,
            chain_id_map,
        }
    }

    #[test]
    fn atoms_returns_all_atoms() {
        let system = mock_system();
        assert_eq!(system.atoms().len(), 2);
    }

    #[test]
    fn chains_returns_all_chains() {
        let system = mock_system();
        assert_eq!(system.chains().len(), 1);
    }

    #[test]
    fn bonds_returns_all_bonds() {
        let system = mock_system();
        assert_eq!(system.bonds().len(), 1);
    }

    #[test]
    fn bonds_mut_allows_modification() {
        let mut system = mock_system();
        system.bonds_mut().clear();
        assert_eq!(system.bonds().len(), 0);
    }

    #[test]
    fn get_atom_returns_correct_atom() {
        let system = mock_system();
        assert_eq!(system.get_atom(0).unwrap().serial, 10);
        assert_eq!(system.get_atom(1).unwrap().serial, 20);
        assert!(system.get_atom(2).is_none());
    }

    #[test]
    fn get_atom_mut_allows_modification() {
        let mut system = mock_system();
        if let Some(atom) = system.get_atom_mut(0) {
            atom.serial = 99;
        }
        assert_eq!(system.get_atom(0).unwrap().serial, 99);
    }

    #[test]
    fn get_atom_by_serial_returns_correct_atom() {
        let system = mock_system();
        assert_eq!(system.get_atom_by_serial(10).unwrap().index, 0);
        assert_eq!(system.get_atom_by_serial(20).unwrap().index, 1);
        assert!(system.get_atom_by_serial(999).is_none());
    }

    #[test]
    fn get_chain_returns_correct_chain() {
        let system = mock_system();
        assert!(system.get_chain(0).is_some());
        assert!(system.get_chain(1).is_none());
    }

    #[test]
    fn get_chain_by_id_returns_correct_chain() {
        let system = mock_system();
        assert!(system.get_chain_by_id('A').is_some());
        assert!(system.get_chain_by_id('B').is_none());
    }

    #[test]
    fn atoms_in_residue_iterates_correct_atoms() {
        let system = mock_system();
        let residue = &system.chains()[0].residues()[0];
        let indices: Vec<_> = system.atoms_in_residue(residue).map(|a| a.index).collect();
        assert_eq!(indices, vec![0, 1]);
    }
}
