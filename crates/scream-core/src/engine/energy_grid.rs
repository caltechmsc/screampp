use crate::core::forcefield::params::Forcefield;
use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::atom::AtomRole;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::engine::cache::ELCache;
use crate::engine::error::EngineError;
use crate::engine::transaction::SystemView;
use itertools::Itertools;
use std::collections::{HashMap, HashSet};
use tracing::{info, trace};

#[derive(Debug, Clone)]
pub struct EnergyGrid {
    pair_interactions: HashMap<(ResidueId, ResidueId), EnergyTerm>,
    total_residue_interactions: HashMap<ResidueId, EnergyTerm>,
    total_interaction_energy: EnergyTerm,
    current_el_energies: HashMap<ResidueId, EnergyTerm>,
    current_optimization_score: f64,
}

#[derive(Debug, Clone)]
pub struct MoveDelta {
    pub res_id: ResidueId,
    pub new_rotamer_idx: usize,
    pub new_el: EnergyTerm,
    pub new_total_interaction: EnergyTerm,
    pub new_pair_interactions: HashMap<ResidueId, EnergyTerm>,
    pub delta_score: f64,
}

impl EnergyGrid {
    pub fn new(
        system: &MolecularSystem,
        forcefield: &Forcefield,
        active_residues: &HashSet<ResidueId>,
        el_cache: &ELCache,
        initial_rotamers: &HashMap<ResidueId, usize>,
    ) -> Result<Self, EngineError> {
        info!("Initializing EnergyGrid with full energy calculation...");

        let mut pair_interactions = HashMap::new();
        let mut total_residue_interactions: HashMap<ResidueId, EnergyTerm> = active_residues
            .iter()
            .map(|&id| (id, EnergyTerm::default()))
            .collect();

        let scorer = Scorer::new(system, forcefield);

        for pair in active_residues.iter().combinations(2) {
            let res_a_id = *pair[0];
            let res_b_id = *pair[1];

            let atoms_a = crate::engine::tasks::interaction_energy::collect_active_sidechain_atoms(
                system,
                &HashSet::from([res_a_id]),
            );
            let atoms_b = crate::engine::tasks::interaction_energy::collect_active_sidechain_atoms(
                system,
                &HashSet::from([res_b_id]),
            );

            let atoms_a_slice = atoms_a
                .get(&res_a_id)
                .map_or([].as_slice(), |v| v.as_slice());
            let atoms_b_slice = atoms_b
                .get(&res_b_id)
                .map_or([].as_slice(), |v| v.as_slice());

            if atoms_a_slice.is_empty() || atoms_b_slice.is_empty() {
                continue;
            }

            let interaction = scorer.score_interaction(atoms_a_slice, atoms_b_slice)?;

            let key = if res_a_id < res_b_id {
                (res_a_id, res_b_id)
            } else {
                (res_b_id, res_a_id)
            };
            pair_interactions.insert(key, interaction);

            *total_residue_interactions.get_mut(&res_a_id).unwrap() += interaction;
            *total_residue_interactions.get_mut(&res_b_id).unwrap() += interaction;
        }

        let total_interaction_energy = total_residue_interactions
            .values()
            .fold(EnergyTerm::default(), |acc, term| acc + *term)
            * 0.5;

        let mut current_el_energies = HashMap::with_capacity(active_residues.len());
        let mut total_el_energy = EnergyTerm::default();

        for &residue_id in active_residues {
            let residue = system.residue(residue_id).unwrap();
            if let (Some(residue_type), Some(rotamer_idx)) =
                (residue.residue_type, initial_rotamers.get(&residue_id))
            {
                if let Some(el_energy) = el_cache.get(residue_id, residue_type, *rotamer_idx) {
                    current_el_energies.insert(residue_id, *el_energy);
                    total_el_energy += *el_energy;
                } else {
                    current_el_energies.insert(residue_id, EnergyTerm::default());
                }
            }
        }

        let current_optimization_score = total_interaction_energy.total() + total_el_energy.total();

        info!(
            "EnergyGrid initialized. Total optimization score: {:.4}",
            current_optimization_score
        );

        Ok(Self {
            pair_interactions,
            total_residue_interactions,
            total_interaction_energy,
            current_el_energies,
            current_optimization_score,
        })
    }

    pub fn total_score(&self) -> f64 {
        self.current_optimization_score
    }

    pub fn calculate_delta_for_move<'a, 'ctx, C>(
        &self,
        res_id: ResidueId,
        new_rotamer_idx: usize,
        system_view: &mut SystemView<'a, 'ctx, C>,
        el_cache: &ELCache,
        active_residues: &HashSet<ResidueId>,
    ) -> Result<MoveDelta, EngineError>
    where
        C: crate::engine::context::ProvidesResidueSelections + Sync,
    {
        let (new_total_interaction, new_pair_interactions) =
            system_view.transaction(res_id, |view| {
                view.apply_move(res_id, new_rotamer_idx)?;

                let scorer = Scorer::new(view.system, view.context.forcefield);

                let new_sc_atoms = view
                    .system
                    .residue(res_id)
                    .unwrap()
                    .atoms()
                    .iter()
                    .filter_map(|&id| {
                        view.system.atom(id).and_then(|a| {
                            if a.role == AtomRole::Sidechain {
                                Some(id)
                            } else {
                                None
                            }
                        })
                    })
                    .collect::<Vec<_>>();

                let mut interaction_sum = EnergyTerm::default();
                let mut pair_interactions_map = HashMap::new();

                for &other_res_id in active_residues {
                    if res_id == other_res_id {
                        continue;
                    }

                    let other_sc_atoms = view
                        .system
                        .residue(other_res_id)
                        .unwrap()
                        .atoms()
                        .iter()
                        .filter_map(|&id| {
                            view.system.atom(id).and_then(|a| {
                                if a.role == AtomRole::Sidechain {
                                    Some(id)
                                } else {
                                    None
                                }
                            })
                        })
                        .collect::<Vec<_>>();

                    let pair_interaction =
                        scorer.score_interaction(&new_sc_atoms, &other_sc_atoms)?;
                    interaction_sum += pair_interaction;
                    pair_interactions_map.insert(other_res_id, pair_interaction);
                }
                Ok((interaction_sum, pair_interactions_map))
            })?;

        let residue_type = system_view
            .system
            .residue(res_id)
            .unwrap()
            .residue_type
            .unwrap();

        let old_el = *self
            .current_el_energies
            .get(&res_id)
            .unwrap_or(&EnergyTerm::default());
        let new_el = *el_cache
            .get(res_id, residue_type, new_rotamer_idx)
            .unwrap_or(&EnergyTerm::default());
        let delta_el = new_el - old_el;

        let old_total_interaction = *self
            .total_residue_interactions
            .get(&res_id)
            .unwrap_or(&EnergyTerm::default());
        let delta_interaction = new_total_interaction - old_total_interaction;
        let delta_score = delta_el.total() + delta_interaction.total();

        trace!(
            "ΔE calculation for res {:?}, rot {}: ΔE_total={:.2} (ΔE_EL={:.2}, ΔE_int={:.2})",
            res_id,
            new_rotamer_idx,
            delta_score,
            delta_el.total(),
            delta_interaction.total()
        );

        Ok(MoveDelta {
            res_id,
            new_rotamer_idx,
            new_el,
            new_total_interaction,
            new_pair_interactions,
            delta_score,
        })
    }

    pub fn apply_move(&mut self, move_delta: MoveDelta) {
        let res_id = move_delta.res_id;

        // 1. Update total score with the pre-calculated delta
        self.current_optimization_score += move_delta.delta_score;

        // 2. Get old interaction terms for diffing
        let _old_el = self
            .current_el_energies
            .get(&res_id)
            .copied()
            .unwrap_or_default();
        let old_total_interaction = self
            .total_residue_interactions
            .get(&res_id)
            .copied()
            .unwrap_or_default();

        // 3. Update EL energy map
        self.current_el_energies.insert(res_id, move_delta.new_el);

        // 4. Update total interaction energy (global sum)
        self.total_interaction_energy = self.total_interaction_energy - old_total_interaction
            + move_delta.new_total_interaction;

        // 5. Update the interaction sum for the moved residue itself
        self.total_residue_interactions
            .insert(res_id, move_delta.new_total_interaction);

        // 6. Update interaction sums for all *other* residues and the pair map
        for (other_res_id, pair_interaction_with_new) in move_delta.new_pair_interactions {
            let key = if res_id < other_res_id {
                (res_id, other_res_id)
            } else {
                (other_res_id, res_id)
            };
            let old_pair_interaction = self
                .pair_interactions
                .get(&key)
                .copied()
                .unwrap_or_default();

            // Update the other residue's total interaction
            if let Some(other_total_int) = self.total_residue_interactions.get_mut(&other_res_id) {
                *other_total_int =
                    *other_total_int - old_pair_interaction + pair_interaction_with_new;
            }

            self.pair_interactions
                .insert(key, pair_interaction_with_new);
        }

        trace!(
            "Applied move for res {:?}. New total score: {:.4}",
            res_id, self.current_optimization_score
        );
    }
}

impl std::ops::Sub for EnergyTerm {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            vdw: self.vdw - rhs.vdw,
            coulomb: self.coulomb - rhs.coulomb,
            hbond: self.hbond - rhs.hbond,
        }
    }
}

impl std::ops::Mul<f64> for EnergyTerm {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            vdw: self.vdw * rhs,
            coulomb: self.coulomb * rhs,
            hbond: self.hbond * rhs,
        }
    }
}
