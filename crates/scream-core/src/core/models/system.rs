use super::atom::Atom;
use super::chain::Chain;
use super::residue::Residue;
use super::topology::Bond;

#[derive(Debug, Clone)]
pub struct MolecularSystem {
    pub(crate) atoms: Vec<Atom>,
    pub(crate) chains: Vec<Chain>,
    pub(crate) bonds: Vec<Bond>,
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
        self.atoms.iter().find(|a| a.serial == serial)
    }

    pub fn get_chain(&self, index: usize) -> Option<&Chain> {
        self.chains.get(index)
    }

    pub fn get_chain_by_id(&self, id: char) -> Option<&Chain> {
        self.chains.iter().find(|c| c.id == id)
    }

    pub fn atoms_in_residue<'a>(&'a self, residue: &'a Residue) -> impl Iterator<Item = &'a Atom> {
        residue
            .atom_indices
            .iter()
            .map(move |&idx| &self.atoms[idx])
    }
}
