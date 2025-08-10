use crate::core::forcefield::params::Forcefield;
use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::atom::AtomRole;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::core::models::system::MolecularSystem;
use crate::engine::error::EngineError;
use itertools::Itertools;
use std::collections::{HashMap, HashSet};
use tracing::instrument;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[instrument(skip_all, name = "interaction_energy_task")]
pub fn run(
    system: &MolecularSystem,
    forcefield: &Forcefield,
    active_residues: &HashSet<ResidueId>,
) -> Result<EnergyTerm, EngineError> {
    if active_residues.len() < 2 {
        return Ok(EnergyTerm::default());
    }

    let sidechain_atoms_map = collect_active_sidechain_atoms(system, active_residues);

    let active_residue_vec: Vec<_> = active_residues.iter().collect();
    let residue_pairs = active_residue_vec.into_iter().combinations(2);

    #[cfg(not(feature = "parallel"))]
    let iterator = residue_pairs;

    #[cfg(feature = "parallel")]
    let iterator = residue_pairs.par_bridge();

    let total_interaction_energy = iterator
        .map(|pair| {
            let res_a_id = *pair[0];
            let res_b_id = *pair[1];

            let atoms_a = sidechain_atoms_map
                .get(&res_a_id)
                .map_or([].as_slice(), |v| v.as_slice());
            let atoms_b = sidechain_atoms_map
                .get(&res_b_id)
                .map_or([].as_slice(), |v| v.as_slice());

            if atoms_a.is_empty() || atoms_b.is_empty() {
                return Ok(EnergyTerm::default());
            }

            let scorer = Scorer::new(system, forcefield);
            scorer.score_interaction(atoms_a, atoms_b)
        })
        .try_reduce(EnergyTerm::default, |acc, term| Ok(acc + term))?;

    Ok(total_interaction_energy)
}

fn collect_active_sidechain_atoms(
    system: &MolecularSystem,
    active_residues: &HashSet<ResidueId>,
) -> HashMap<ResidueId, Vec<AtomId>> {
    let mut map = HashMap::with_capacity(active_residues.len());
    for &residue_id in active_residues {
        if let Some(residue) = system.residue(residue_id) {
            let sidechain_ids: Vec<AtomId> = residue
                .atoms()
                .iter()
                .filter_map(|&atom_id| {
                    system.atom(atom_id).and_then(|atom| {
                        if atom.role == AtomRole::Sidechain {
                            Some(atom_id)
                        } else {
                            None
                        }
                    })
                })
                .collect();
            map.insert(residue_id, sidechain_ids);
        }
    }
    map
}
