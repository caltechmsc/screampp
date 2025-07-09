use super::energy::EnergyCalculator;
use super::params::Forcefield;
use crate::core::models::ids::AtomId;
use crate::core::models::system::MolecularSystem;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ScoringError {
    #[error("Atom with ID {0:?} not found in the system")]
    AtomNotFound(AtomId),
    #[error("Force field type not parameterized for atom {0:?}")]
    ForceFieldTypeMissing(AtomId),
    #[error("Could not find donor for hydrogen atom {0:?}")]
    DonorNotFound(AtomId),
}

pub struct Scorer<'a> {
    system: &'a MolecularSystem,
    forcefield: &'a Forcefield,
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct InteractionEnergy {
    pub vdw: f64,
    pub coulomb: f64,
    pub hbond: f64,
}

impl InteractionEnergy {
    pub fn total(&self) -> f64 {
        self.vdw + self.coulomb + self.hbond
    }
}

impl<'a> Scorer<'a> {
    pub fn new(system: &'a MolecularSystem, forcefield: &'a Forcefield) -> Self {
        Self { system, forcefield }
    }

    pub fn score_interaction(
        &self,
        query_atom_ids: &[AtomId],
        environment_atom_ids: &[AtomId],
    ) -> Result<InteractionEnergy, ScoringError> {
        let mut energy = InteractionEnergy::default();

        for &query_id in query_atom_ids {
            let query_atom = self
                .system
                .atom(query_id)
                .ok_or(ScoringError::AtomNotFound(query_id))?;

            for &env_id in environment_atom_ids {
                if query_id == env_id {
                    continue;
                }

                let env_atom = self
                    .system
                    .atom(env_id)
                    .ok_or(ScoringError::AtomNotFound(env_id))?;

                if query_atom.residue_id == env_atom.residue_id {
                    continue;
                }

                let vdw_param1 = self
                    .forcefield
                    .non_bonded
                    .vdw
                    .get(&query_atom.force_field_type)
                    .ok_or_else(|| ScoringError::ForceFieldTypeMissing(query_id))?;
                let vdw_param2 = self
                    .forcefield
                    .non_bonded
                    .vdw
                    .get(&env_atom.force_field_type)
                    .ok_or_else(|| ScoringError::ForceFieldTypeMissing(env_id))?;
                energy.vdw +=
                    EnergyCalculator::calculate_vdw(query_atom, env_atom, vdw_param1, vdw_param2);

                energy.coulomb += EnergyCalculator::calculate_coulomb(
                    query_atom,
                    env_atom,
                    self.forcefield.non_bonded.globals.dielectric_constant,
                );
            }
        }

        let all_ids: Vec<_> = query_atom_ids
            .iter()
            .chain(environment_atom_ids.iter())
            .copied()
            .collect();

        for &h_id in &all_ids {
            let hydrogen = self
                .system
                .atom(h_id)
                .ok_or(ScoringError::AtomNotFound(h_id))?;

            if hydrogen.hbond_type_id == 0 {
                let donor_id = *self
                    .system
                    .get_bonded_neighbors(h_id)
                    .and_then(|neighbors| neighbors.first())
                    .ok_or(ScoringError::DonorNotFound(h_id))?;
                let donor = self
                    .system
                    .atom(donor_id)
                    .ok_or(ScoringError::AtomNotFound(donor_id))?;

                for &a_id in &all_ids {
                    if a_id == h_id || a_id == donor_id {
                        continue;
                    }

                    let acceptor = self
                        .system
                        .atom(a_id)
                        .ok_or(ScoringError::AtomNotFound(a_id))?;

                    if acceptor.hbond_type_id > 0 {
                        let is_query_h = query_atom_ids.contains(&h_id);
                        let is_env_a = environment_atom_ids.contains(&a_id);
                        let is_env_h = environment_atom_ids.contains(&h_id);
                        let is_query_a = query_atom_ids.contains(&a_id);

                        if !((is_query_h && is_env_a) || (is_env_h && is_query_a)) {
                            continue;
                        }

                        let hbond_key =
                            format!("{}-{}", acceptor.force_field_type, donor.force_field_type);
                        if let Some(hbond_param) = self.forcefield.non_bonded.hbond.get(&hbond_key)
                        {
                            energy.hbond += EnergyCalculator::calculate_hbond(
                                acceptor,
                                hydrogen,
                                donor,
                                hbond_param.equilibrium_dist,
                                hbond_param.well_depth,
                            );
                        }
                    }
                }
            }
        }

        Ok(energy)
    }
}
