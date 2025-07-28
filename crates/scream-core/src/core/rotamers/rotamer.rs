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
    pub bonds: Vec<[usize; 2]>,
}

#[derive(Debug, Clone)]
pub struct Rotamer {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<(usize, usize)>,
}
