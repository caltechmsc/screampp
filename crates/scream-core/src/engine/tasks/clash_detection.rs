use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::engine::context::{Context, ProvidesResidueSelections};
use crate::engine::error::EngineError;
use crate::engine::progress::Progress;
use itertools::Itertools;
use std::cmp::Ordering;
use tracing::{info, instrument};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Debug, Clone, PartialEq)]
pub struct ClashPair {
    pub residue_a: ResidueId,
    pub residue_b: ResidueId,
    pub energy: EnergyTerm,
}

impl Eq for ClashPair {}

impl PartialOrd for ClashPair {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        other.energy.total().partial_cmp(&self.energy.total())
    }
}

impl Ord for ClashPair {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}
