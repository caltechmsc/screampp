use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::atom::AtomRole;
use crate::core::models::ids::AtomId;
use crate::engine::context::{OptimizationContext, ProvidesResidueSelections};
use crate::engine::error::EngineError;
use tracing::{info, instrument};

struct FixedAtomSets {
    sidechain_atoms: Vec<AtomId>,
    non_sidechain_atoms: Vec<AtomId>,
}

#[instrument(skip_all, name = "fixed_energy_task")]
pub fn run<C>(context: &OptimizationContext<C>) -> Result<EnergyTerm, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    info!("Calculating constant energy offset for the fixed parts of the system.");

    let fixed_sets = collect_fixed_atom_sets(context)?;

    if fixed_sets.sidechain_atoms.is_empty() && fixed_sets.non_sidechain_atoms.is_empty() {
        info!("No fixed atoms found in the system; the energy offset is zero.");
        return Ok(EnergyTerm::default());
    }

    let scorer = Scorer::new(context.system, context.forcefield);

    let energy1 = scorer.score_group_internal(&fixed_sets.non_sidechain_atoms)?;

    let energy2 =
        scorer.score_interaction(&fixed_sets.sidechain_atoms, &fixed_sets.non_sidechain_atoms)?;

    // TODO: Implement summation of `E_internal_or_fallback` for `fixed_sets.sidechain_atoms`.
    let energy3 = EnergyTerm::default();

    let total_offset_energy = energy1 + energy2 + energy3;
    info!(
        energy = total_offset_energy.total(),
        vdw = total_offset_energy.vdw,
        coulomb = total_offset_energy.coulomb,
        hbond = total_offset_energy.hbond,
        "Fixed energy offset calculation complete."
    );

    Ok(total_offset_energy)
}

fn collect_fixed_atom_sets<C>(
    context: &OptimizationContext<C>,
) -> Result<FixedAtomSets, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    let active_residues = context.resolve_all_active_residues()?;

    let mut sidechain_atoms = Vec::new();
    let mut non_sidechain_atoms = Vec::new();

    for (atom_id, atom) in context.system.atoms_iter() {
        if active_residues.contains(&atom.residue_id) {
            if atom.role == AtomRole::Backbone {
                non_sidechain_atoms.push(atom_id);
            }
        } else {
            match atom.role {
                AtomRole::Sidechain => sidechain_atoms.push(atom_id),
                _ => non_sidechain_atoms.push(atom_id),
            }
        }
    }

    Ok(FixedAtomSets {
        sidechain_atoms,
        non_sidechain_atoms,
    })
}
