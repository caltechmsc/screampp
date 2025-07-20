use crate::core::forcefield::scoring::Scorer;
use crate::core::models::ids::ResidueId;
use crate::engine::cache::ELCache;
use crate::engine::context::{Context, ProvidesResidueSelections};
use crate::engine::error::EngineError;
use crate::engine::placement;
use tracing::{debug, instrument, trace};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Debug, Clone, Copy)]
pub struct DoubletResult {
    pub rotamer_idx_a: usize,
    pub rotamer_idx_b: usize,
    pub best_local_energy: f64,
}
