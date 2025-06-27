use super::atom::Atom;
use std::collections::HashMap;

pub struct Residue {
    pub id: usize,                         // Unique identifier for the residue
    pub name: String,                      // Name of the residue (e.g., "ALA", "GLY")
    atoms: Vec<Atom>,                      // List of atoms in the residue
    atom_map: HashMap<String, Vec<usize>>, // Map from atom name to list of atom indices
}

impl Residue {
    pub fn new(id: usize, name: &str, atoms: Vec<Atom>) -> Self {
        let atoms_vec = atoms;
        let mut atom_map = HashMap::new();
        for (i, atom) in atoms_vec.iter().enumerate() {
            atom_map
                .entry(atom.name.clone())
                .or_insert_with(Vec::new)
                .push(i);
        }
        Self {
            id,
            name: name.to_string(),
            atoms: atoms_vec,
            atom_map,
        }
    }

    pub fn add_atom(&mut self, atom: Atom) {
        let index = self.atoms.len();
        self.atom_map
            .entry(atom.name.clone())
            .or_insert_with(Vec::new)
            .push(index);
        self.atoms.push(atom);
    }

    pub fn get_atom(&self, index: usize) -> Option<&Atom> {
        self.atoms.get(index)
    }

    pub fn get_atoms_by_name(&self, name: &str) -> impl Iterator<Item = &Atom> {
        self.atom_map
            .get(name)
            .map(|indices| indices.iter().map(move |&i| &self.atoms[i]))
            .into_iter()
            .flatten()
    }

    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }
}
