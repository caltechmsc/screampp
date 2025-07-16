use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::ids::ResidueId;
use crate::core::models::residue::ResidueType;
use crate::engine::cache::ELCache;
use crate::engine::config::DesignSpecExt;
use crate::engine::context::{Context, ProvidesResidueSelections};
use crate::engine::error::EngineError;
use crate::engine::placement::place_rotamer_on_system;
use crate::engine::progress::Progress;
use tracing::{info, instrument, warn};

#[cfg(feature = "parallel")]
use rayon::prelude::*;
