use crate::core::models::atom::Atom;
use serde::Deserialize;

#[derive(Debug, Clone, Deserialize)]
pub struct RotamerAtomData {
    pub serial: usize,
    pub atom_name: String,
    pub partial_charge: f64,
    pub position: [f64; 3],
    pub force_field_type: String,
}

#[derive(Debug, Clone, Deserialize)]
pub struct RotamerData {
    pub atoms: Vec<RotamerAtomData>,
    pub bonds: Vec<(usize, usize)>,
}

#[derive(Debug, Clone)]
pub struct Rotamer {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<(usize, usize)>,
}

impl Rotamer {
    pub fn new(atoms: Vec<Atom>, bonds: Vec<(usize, usize)>) -> Self {
        Self { atoms, bonds }
    }
}
