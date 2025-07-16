use crate::core::forcefield::term::EnergyTerm;
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
}

#[derive(Debug, Clone)]
pub struct Rotamer {
    pub atoms: Vec<Atom>,
    pub empty_lattice_energy: Option<EnergyTerm>,
}
