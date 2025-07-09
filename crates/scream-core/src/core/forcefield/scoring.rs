use super::energy::EnergyCalculator;
use super::params::Forcefield;
use crate::core::models::ids::AtomId;
use crate::core::models::system::MolecularSystem;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ScoringError {
    #[error("Atom with ID {0:?} not found in the system")]
    AtomNotFound(AtomId),
    #[error("Force field type not parameterized for atom {0:?}")]
    ForceFieldTypeMissing(AtomId),
    #[error("Could not find donor for hydrogen atom {0:?}")]
    DonorNotFound(AtomId),
}

pub struct Scorer<'a> {
    system: &'a MolecularSystem,
    forcefield: &'a Forcefield,
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct InteractionEnergy {
    pub vdw: f64,
    pub coulomb: f64,
    pub hbond: f64,
}

impl InteractionEnergy {
    pub fn total(&self) -> f64 {
        self.vdw + self.coulomb + self.hbond
    }
}

impl<'a> Scorer<'a> {
    pub fn new(system: &'a MolecularSystem, forcefield: &'a Forcefield) -> Self {
        Self { system, forcefield }
    }
}
