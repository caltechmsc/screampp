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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::params::{GlobalParams, HBondParam, NonBondedParams, VdwParam};
    use crate::core::models::atom::Atom;
    use crate::core::models::ids::ResidueId;
    use crate::core::models::topology::BondOrder;
    use nalgebra::Point3;
    use std::collections::HashMap;

    fn create_test_forcefield() -> Forcefield {
        let globals = GlobalParams {
            dielectric_constant: 1.0,
            potential_function: "mixed".to_string(),
        };

        let mut vdw = HashMap::new();
        vdw.insert(
            "C".to_string(),
            VdwParam::LennardJones {
                radius: 4.0,
                well_depth: 0.1,
            },
        );
        vdw.insert(
            "N".to_string(),
            VdwParam::LennardJones {
                radius: 3.5,
                well_depth: 0.2,
            },
        );
        vdw.insert(
            "O".to_string(),
            VdwParam::LennardJones {
                radius: 3.2,
                well_depth: 0.3,
            },
        );
        vdw.insert(
            "H".to_string(),
            VdwParam::LennardJones {
                radius: 1.0,
                well_depth: 0.01,
            },
        );

        let mut hbond = HashMap::new();
        hbond.insert(
            "O-N".to_string(),
            HBondParam {
                equilibrium_dist: 2.8,
                well_depth: 5.0,
            },
        );

        let non_bonded = NonBondedParams {
            globals,
            vdw,
            hbond,
        };

        Forcefield {
            non_bonded,
            deltas: HashMap::new(),
        }
    }

    fn create_atom(
        serial: usize,
        name: &str,
        residue_id: ResidueId,
        pos: Point3<f64>,
        ff_type: &str,
        charge: f64,
        hbond_type_id: i8,
    ) -> Atom {
        let mut atom = Atom::new(serial, name, residue_id, pos);
        atom.force_field_type = ff_type.to_string();
        atom.partial_charge = charge;
        atom.hbond_type_id = hbond_type_id;
        atom
    }

    #[test]
    fn scores_vdw_and_coulomb_for_simple_query_and_environment() {
        let mut system = MolecularSystem::new();
        let ff = create_test_forcefield();

        let chain_id = system.add_chain('A', crate::core::models::chain::ChainType::Protein);
        let res1_id = system.add_residue(chain_id, 1, "RES", None).unwrap();
        let res2_id = system.add_residue(chain_id, 2, "RES", None).unwrap();

        let query_atom = create_atom(1, "C1", res1_id, Point3::origin(), "C", 0.5, -1);
        let env_atom = create_atom(2, "C2", res2_id, Point3::new(4.5, 0.0, 0.0), "C", -0.5, -1);

        let query_id = system.add_atom_to_residue(res1_id, query_atom).unwrap();
        let env_id = system.add_atom_to_residue(res2_id, env_atom).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let energy = scorer.score_interaction(&[query_id], &[env_id]).unwrap();

        assert!(energy.vdw < 0.0);
        assert!(energy.coulomb < 0.0);
        assert_eq!(energy.hbond, 0.0);
    }

    #[test]
    fn ignores_interactions_within_the_same_residue() {
        let mut system = MolecularSystem::new();
        let ff = create_test_forcefield();

        let chain_id = system.add_chain('A', crate::core::models::chain::ChainType::Protein);
        let res1_id = system.add_residue(chain_id, 1, "RES", None).unwrap();

        let atom1 = create_atom(1, "C1", res1_id, Point3::origin(), "C", 0.5, -1);
        let atom2 = create_atom(2, "C2", res1_id, Point3::new(3.0, 0.0, 0.0), "C", -0.5, -1);

        let id1 = system.add_atom_to_residue(res1_id, atom1).unwrap();
        let id2 = system.add_atom_to_residue(res1_id, atom2).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let energy = scorer.score_interaction(&[id1], &[id2]).unwrap();

        assert_eq!(energy.total(), 0.0);
    }

    #[test]
    fn scores_hydrogen_bond_between_query_and_environment() {
        let mut system = MolecularSystem::new();
        let ff = create_test_forcefield();

        let chain_id = system.add_chain('A', crate::core::models::chain::ChainType::Protein);
        let res1_id = system.add_residue(chain_id, 1, "RES", None).unwrap();
        let res2_id = system.add_residue(chain_id, 2, "RES", None).unwrap();

        let donor = create_atom(1, "N", res1_id, Point3::new(0.0, 0.0, 0.0), "N", -0.3, 1);
        let hydrogen = create_atom(2, "H", res1_id, Point3::new(1.0, 0.0, 0.0), "H", 0.3, 0);
        let donor_id = system.add_atom_to_residue(res1_id, donor).unwrap();
        let h_id = system.add_atom_to_residue(res1_id, hydrogen).unwrap();
        system.add_bond(donor_id, h_id, BondOrder::Single).unwrap();

        let acceptor = create_atom(3, "O", res2_id, Point3::new(2.8, 0.0, 0.0), "O", -0.5, 1);
        let acceptor_id = system.add_atom_to_residue(res2_id, acceptor).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let energy = scorer
            .score_interaction(&[acceptor_id], &[donor_id, h_id])
            .unwrap();

        assert!(energy.hbond < 0.0);
        assert_ne!(energy.vdw, 0.0);
        assert_ne!(energy.coulomb, 0.0);
    }

    #[test]
    fn ignores_hydrogen_bond_within_query_set() {
        let mut system = MolecularSystem::new();
        let ff = create_test_forcefield();

        let chain_id = system.add_chain('A', crate::core::models::chain::ChainType::Protein);
        let res1_id = system.add_residue(chain_id, 1, "RES", None).unwrap();
        let res2_id = system.add_residue(chain_id, 2, "RES", None).unwrap();

        let donor = create_atom(1, "N", res1_id, Point3::new(0.0, 0.0, 0.0), "N", -0.3, 1);
        let hydrogen = create_atom(2, "H", res1_id, Point3::new(1.0, 0.0, 0.0), "H", 0.3, 0);
        let acceptor = create_atom(3, "O", res2_id, Point3::new(2.8, 0.0, 0.0), "O", -0.5, 1);

        let donor_id = system.add_atom_to_residue(res1_id, donor).unwrap();
        let h_id = system.add_atom_to_residue(res1_id, hydrogen).unwrap();
        let acceptor_id = system.add_atom_to_residue(res2_id, acceptor).unwrap();
        system.add_bond(donor_id, h_id, BondOrder::Single).unwrap();

        let res3_id = system.add_residue(chain_id, 3, "RES", None).unwrap();
        let env_atom = create_atom(4, "C", res3_id, Point3::new(5.0, 5.0, 5.0), "C", 0.1, -1);
        let env_id = system.add_atom_to_residue(res3_id, env_atom).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let energy = scorer
            .score_interaction(&[donor_id, h_id, acceptor_id], &[env_id])
            .unwrap();

        assert_eq!(energy.hbond, 0.0);
        assert_ne!(energy.vdw, 0.0);
        assert_ne!(energy.coulomb, 0.0);
    }

    #[test]
    fn returns_error_for_missing_atom() {
        let system = MolecularSystem::new();
        let ff = create_test_forcefield();
        let scorer = Scorer::new(&system, &ff);

        let fake_id = AtomId::default();
        let result = scorer.score_interaction(&[fake_id], &[]);

        assert!(matches!(result, Err(ScoringError::AtomNotFound(_))));
    }

    #[test]
    fn returns_error_for_unparameterized_force_field_type() {
        let mut system = MolecularSystem::new();
        let ff = create_test_forcefield();

        let chain_id = system.add_chain('A', crate::core::models::chain::ChainType::Protein);
        let res1_id = system.add_residue(chain_id, 1, "RES", None).unwrap();
        let res2_id = system.add_residue(chain_id, 2, "RES", None).unwrap();

        let query_atom = create_atom(1, "C1", res1_id, Point3::origin(), "Unknown", 0.0, -1);
        let env_atom = create_atom(2, "C2", res2_id, Point3::new(3.0, 0.0, 0.0), "C", 0.0, -1);

        let query_id = system.add_atom_to_residue(res1_id, query_atom).unwrap();
        let env_id = system.add_atom_to_residue(res2_id, env_atom).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let result = scorer.score_interaction(&[query_id], &[env_id]);

        assert!(matches!(
            result,
            Err(ScoringError::ForceFieldTypeMissing(_))
        ));
    }

    #[test]
    fn returns_error_for_hydrogen_without_donor() {
        let mut system = MolecularSystem::new();
        let ff = create_test_forcefield();

        let chain_id = system.add_chain('A', crate::core::models::chain::ChainType::Protein);
        let res1_id = system.add_residue(chain_id, 1, "RES", None).unwrap();

        let hydrogen = create_atom(1, "H", res1_id, Point3::new(1.0, 0.0, 0.0), "H", 0.3, 0);
        let h_id = system.add_atom_to_residue(res1_id, hydrogen).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let result = scorer.score_interaction(&[h_id], &[]);

        assert!(matches!(result, Err(ScoringError::DonorNotFound(_))));
    }
}
