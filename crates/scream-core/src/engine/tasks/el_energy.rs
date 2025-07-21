use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::core::models::residue::ResidueType;
use crate::engine::cache::ELCache;
use crate::engine::config::DesignSpecExt;
use crate::engine::context::{OptimizationContext, ProvidesResidueSelections};
use crate::engine::error::EngineError;
use crate::engine::placement::place_rotamer_on_system;
use crate::engine::progress::Progress;
use std::collections::HashMap;
use tracing::{info, instrument, warn};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Debug)]
struct WorkUnit {
    residue_id: ResidueId,
    residue_type: ResidueType,
}

type WorkResult = Result<((ResidueId, ResidueType), HashMap<usize, EnergyTerm>), EngineError>;

#[instrument(skip_all, name = "el_energy_task")]
pub fn run<C>(context: &OptimizationContext<C>) -> Result<ELCache, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    info!("Starting Empty Lattice energy pre-computation.");
    context.reporter.report(Progress::PhaseStart {
        name: "EL Pre-computation",
    });

    let work_list = build_work_list(context)?;

    if work_list.is_empty() {
        warn!("No work to be done for EL energy calculation. Returning empty cache.");
        return Ok(ELCache::new());
    }

    context.reporter.report(Progress::TaskStart {
        total_steps: work_list.len() as u64,
    });

    #[cfg(not(feature = "parallel"))]
    let iterator = work_list.iter();

    #[cfg(feature = "parallel")]
    let iterator = work_list.par_iter();

    let results: Vec<WorkResult> = iterator
        .map(|unit| compute_energies_for_unit(unit, context))
        .collect();

    context.reporter.report(Progress::TaskFinish);

    let mut cache = ELCache::new();
    for result in results {
        let ((residue_id, residue_type), energy_map) = result?;
        for (rotamer_idx, energy_term) in energy_map {
            cache.insert(residue_id, residue_type, rotamer_idx, energy_term);
        }
    }

    info!(
        cached_combinations = cache.len(),
        "EL pre-computation finished."
    );
    Ok(cache)
}

fn build_work_list<C>(context: &OptimizationContext<C>) -> Result<Vec<WorkUnit>, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    let mut work_list = Vec::new();
    let active_residues = context.resolve_all_active_residues()?;

    for &residue_id in &active_residues {
        let residue = context.system.residue(residue_id).unwrap();

        let mut is_design_site = false;
        if let Some(design_spec) = context.config.design_spec() {
            let chain = context.system.chain(residue.chain_id).unwrap();
            if let Some(allowed_types) = design_spec.get_by_specifier(chain.id, residue.id) {
                for &residue_type in allowed_types {
                    work_list.push(WorkUnit {
                        residue_id,
                        residue_type,
                    });
                }
                is_design_site = true;
            }
        }

        if !is_design_site {
            if let Some(native_type) = residue.res_type {
                work_list.push(WorkUnit {
                    residue_id,
                    residue_type: native_type,
                });
            }
        }
    }
    Ok(work_list)
}

#[instrument(skip_all, fields(residue_id = ?unit.residue_id, residue_type = %unit.residue_type))]
fn compute_energies_for_unit<C>(unit: &WorkUnit, context: &OptimizationContext<C>) -> WorkResult
where
    C: ProvidesResidueSelections + Sync,
{
    let rotamers = context
        .rotamer_library
        .get_rotamers_for(unit.residue_type)
        .ok_or_else(|| EngineError::RotamerLibrary {
            residue_type: unit.residue_type.to_string(),
            message: "No rotamers found for this residue type.".to_string(),
        })?;

    let placement_info = context
        .rotamer_library
        .get_placement_info_for(unit.residue_type)
        .ok_or_else(|| EngineError::RotamerLibrary {
            residue_type: unit.residue_type.to_string(),
            message: "No placement info found for this residue type.".to_string(),
        })?;

    let active_residue_ids = context.resolve_all_active_residues()?;

    let environment_atoms: Vec<AtomId> = context
        .system
        .atoms_iter()
        .filter_map(|(atom_id, atom)| {
            if !active_residue_ids.contains(&atom.residue_id) {
                Some(atom_id)
            } else {
                None
            }
        })
        .collect();

    let mut energy_map = HashMap::with_capacity(rotamers.len());

    for (rotamer_idx, rotamer) in rotamers.iter().enumerate() {
        let mut temp_system = context.system.clone();

        place_rotamer_on_system(&mut temp_system, unit.residue_id, rotamer, placement_info)?;

        let query_atoms: Vec<AtomId> = temp_system
            .residue(unit.residue_id)
            .unwrap()
            .atoms()
            .to_vec();

        let scorer = Scorer::new(&temp_system, context.forcefield);
        let energy = scorer.score_interaction(&query_atoms, &environment_atoms)?;

        energy_map.insert(rotamer_idx, energy);
    }

    context.reporter.report(Progress::TaskIncrement);

    Ok(((unit.residue_id, unit.residue_type), energy_map))
}
