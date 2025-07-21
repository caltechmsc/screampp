use crate::core::forcefield::params::Forcefield;
use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::engine::error::EngineError;
use crate::engine::progress::{Progress, ProgressReporter};
use itertools::Itertools;
use std::cmp::Ordering;
use std::collections::HashSet;
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

#[instrument(skip_all, name = "clash_detection_task")]
pub fn run(
    system: &MolecularSystem,
    forcefield: &Forcefield,
    active_residues: &HashSet<ResidueId>,
    clash_threshold_kcal_mol: f64,
    reporter: &ProgressReporter,
) -> Result<Vec<ClashPair>, EngineError> {
    info!(
        threshold = clash_threshold_kcal_mol,
        "Detecting residue clashes."
    );
    reporter.report(Progress::Message("Detecting clashes...".to_string()));

    let residue_pairs: Vec<_> = active_residues.iter().combinations(2).collect();

    if residue_pairs.is_empty() {
        return Ok(Vec::new());
    }

    reporter.report(Progress::TaskStart {
        total_steps: residue_pairs.len() as u64,
    });

    let scorer = Scorer::new(system, forcefield);

    #[cfg(not(feature = "parallel"))]
    let iterator = residue_pairs.iter();

    #[cfg(feature = "parallel")]
    let iterator = residue_pairs.par_iter();

    let clashes: Result<Vec<ClashPair>, EngineError> = iterator
        .filter_map(|pair| {
            reporter.report(Progress::TaskIncrement);
            let res_id_a = *pair[0];
            let res_id_b = *pair[1];

            let atoms_a = system.residue(res_id_a).unwrap().atoms();
            let atoms_b = system.residue(res_id_b).unwrap().atoms();

            match scorer.score_interaction(atoms_a, atoms_b) {
                Ok(energy) if energy.total() > clash_threshold_kcal_mol => Some(Ok(ClashPair {
                    residue_a: res_id_a,
                    residue_b: res_id_b,
                    energy,
                })),
                Ok(_) => None,
                Err(e) => Some(Err(EngineError::from(e))),
            }
        })
        .collect();

    let mut clashes = clashes?;

    reporter.report(Progress::TaskFinish);

    clashes.sort_unstable();

    info!(num_clashes = clashes.len(), "Clash detection complete.");

    Ok(clashes)
}
