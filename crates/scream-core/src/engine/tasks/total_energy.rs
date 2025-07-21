use crate::core::forcefield::params::Forcefield;
use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::engine::cache::ELCache;
use crate::engine::error::EngineError;
use itertools::Itertools;
use std::collections::{HashMap, HashSet};
use std::iter::Sum;
use tracing::instrument;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

impl Sum for EnergyTerm {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::default(), |acc, term| acc + term)
    }
}

#[instrument(skip_all, name = "total_energy_task")]
pub fn run(
    system: &MolecularSystem,
    forcefield: &Forcefield,
    active_residues: &HashSet<ResidueId>,
    current_rotamers: &HashMap<ResidueId, usize>,
    el_cache: &ELCache,
) -> Result<EnergyTerm, EngineError> {
    let total_el_energy: EnergyTerm = active_residues
        .iter()
        .map(|&residue_id| {
            let residue = system.residue(residue_id).unwrap();
            let res_type = residue.res_type.unwrap();
            let rotamer_idx = current_rotamers.get(&residue_id).unwrap();

            el_cache
                .get(residue_id, res_type, *rotamer_idx)
                .copied()
                .unwrap_or_default()
        })
        .sum();

    let scorer = Scorer::new(system, forcefield);
    let residue_pairs: Vec<_> = active_residues.iter().combinations(2).collect();

    #[cfg(not(feature = "parallel"))]
    let iterator = residue_pairs.iter();

    #[cfg(feature = "parallel")]
    let iterator = residue_pairs.par_iter();

    let interaction_energy: EnergyTerm = iterator
        .map(|pair| {
            let res_a_id = *pair[0];
            let res_b_id = *pair[1];
            let atoms_a = system.residue(res_a_id).unwrap().atoms();
            let atoms_b = system.residue(res_b_id).unwrap().atoms();
            scorer
                .score_interaction(atoms_a, atoms_b)
                .map_err(EngineError::from)
        })
        .try_fold(EnergyTerm::default(), |acc, term| term.map(|t| acc + t))?;

    Ok(total_el_energy + interaction_energy)
}
