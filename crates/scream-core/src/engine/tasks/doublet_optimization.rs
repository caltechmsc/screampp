use crate::core::forcefield::scoring::Scorer;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::core::models::system::MolecularSystem;
use crate::engine::cache::ELCache;
use crate::engine::context::{OptimizationContext, ProvidesResidueSelections};
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

#[instrument(skip_all, name = "doublet_optimization_task", fields(res_a = ?res_a_id, res_b = ?res_b_id))]
pub fn run<C>(
    res_a_id: ResidueId,
    res_b_id: ResidueId,
    system: &MolecularSystem,
    el_cache: &ELCache,
    context: &OptimizationContext<C>,
) -> Result<DoubletResult, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    let residue_a = system.residue(res_a_id).unwrap();
    let residue_b = system.residue(res_b_id).unwrap();
    let res_type_a = residue_a.res_type.unwrap();
    let res_type_b = residue_b.res_type.unwrap();

    let rotamers_a = context
        .rotamer_library
        .get_rotamers_for(res_type_a)
        .ok_or_else(|| EngineError::RotamerLibrary {
            residue_type: res_type_a.to_string(),
            message: "No rotamers found for doublet optimization.".to_string(),
        })?;

    let rotamers_b = context
        .rotamer_library
        .get_rotamers_for(res_type_b)
        .ok_or_else(|| EngineError::RotamerLibrary {
            residue_type: res_type_b.to_string(),
            message: "No rotamers found for doublet optimization.".to_string(),
        })?;

    if rotamers_a.is_empty() || rotamers_b.is_empty() {
        return Err(EngineError::RotamerLibrary {
            residue_type: if rotamers_a.is_empty() {
                res_type_a.to_string()
            } else {
                res_type_b.to_string()
            },
            message: "Rotamer list is empty, cannot perform doublet optimization.".to_string(),
        });
    }

    let index_pairs: Vec<(usize, usize)> = (0..rotamers_a.len())
        .flat_map(|i| (0..rotamers_b.len()).map(move |j| (i, j)))
        .collect();

    debug!(
        "Optimizing doublet with {}x{} = {} total pairs.",
        rotamers_a.len(),
        rotamers_b.len(),
        index_pairs.len()
    );

    let placement_info_a = context
        .rotamer_library
        .get_placement_info_for(res_type_a)
        .unwrap();
    let placement_info_b = context
        .rotamer_library
        .get_placement_info_for(res_type_b)
        .unwrap();

    #[cfg(not(feature = "parallel"))]
    let iterator = index_pairs.iter();

    #[cfg(feature = "parallel")]
    let iterator = index_pairs.par_iter();

    let best_pair = iterator
        .map(|&(idx_a, idx_b)| {
            let rot_a = &rotamers_a[idx_a];
            let rot_b = &rotamers_b[idx_b];

            let el_a = el_cache
                .get(res_a_id, res_type_a, idx_a)
                .map_or(0.0, |e| e.total());
            let el_b = el_cache
                .get(res_b_id, res_type_b, idx_b)
                .map_or(0.0, |e| e.total());

            let mut temp_system = system.clone();
            placement::place_rotamer_on_system(
                &mut temp_system,
                res_a_id,
                rot_a,
                placement_info_a,
            )?;
            placement::place_rotamer_on_system(
                &mut temp_system,
                res_b_id,
                rot_b,
                placement_info_b,
            )?;

            let atoms_a_sidechain_ids: Vec<AtomId> = placement_info_a
                .sidechain_atoms
                .iter()
                .filter_map(|name| {
                    temp_system
                        .residue(res_a_id)
                        .unwrap()
                        .get_atom_id_by_name(name)
                })
                .collect();
            let atoms_b_sidechain_ids: Vec<AtomId> = placement_info_b
                .sidechain_atoms
                .iter()
                .filter_map(|name| {
                    temp_system
                        .residue(res_b_id)
                        .unwrap()
                        .get_atom_id_by_name(name)
                })
                .collect();

            let scorer = Scorer::new(&temp_system, context.forcefield);
            let interaction = scorer
                .score_interaction(&atoms_a_sidechain_ids, &atoms_b_sidechain_ids)?
                .total();
            let local_energy = el_a + el_b + interaction;

            trace!(
                "Pair ({}, {}): E_local = {:.2} (EL_A={:.2}, EL_B={:.2}, Int={:.2})",
                idx_a, idx_b, local_energy, el_a, el_b, interaction
            );

            Ok((local_energy, (idx_a, idx_b)))
        })
        .collect::<Result<Vec<_>, EngineError>>()?
        .into_iter()
        .min_by(|(energy_1, _), (energy_2, _)| {
            energy_1
                .partial_cmp(energy_2)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

    match best_pair {
        Some((best_energy, (best_idx_a, best_idx_b))) => Ok(DoubletResult {
            rotamer_idx_a: best_idx_a,
            rotamer_idx_b: best_idx_b,
            best_local_energy: best_energy,
        }),
        None => Err(EngineError::PhaseFailed {
            phase: "Doublet Optimization",
            reason: "No valid rotamer pairs could be evaluated.".to_string(),
        }),
    }
}
