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

#[derive(Debug, Deserialize, Clone, Default)]
pub struct EnergyData {
    pub bond: f64,
    pub angle: f64,
    pub torsion: f64,
    pub inversion: f64,
}

#[derive(Debug, Clone, Deserialize)]
pub struct RotamerData {
    pub atoms: Vec<RotamerAtomData>,
    pub bonds: Vec<[usize; 2]>,
    #[serde(default)]
    pub energy: EnergyData,
}

#[derive(Debug, Clone)]
pub struct Rotamer {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<(usize, usize)>,
    pub energy: EnergyTerm,
}
