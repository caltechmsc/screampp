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
