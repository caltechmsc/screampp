use super::energy::{EnergyCalculationError, EnergyCalculator};
use super::params::Forcefield;
use super::term::EnergyTerm;
use crate::core::models::ids::AtomId;
use crate::core::models::system::MolecularSystem;
use std::collections::HashSet;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ScoringError {
    #[error("Atom with ID {0:?} not found in the system")]
    AtomNotFound(AtomId),
    #[error("Could not find donor for hydrogen atom {0:?}")]
    DonorNotFound(AtomId),
    #[error("Energy calculation failed: {source}")]
    EnergyCalculation {
        #[from]
        source: EnergyCalculationError,
    },
}

pub struct Scorer<'a> {
    system: &'a MolecularSystem,
    forcefield: &'a Forcefield,
}

impl<'a> Scorer<'a> {
    pub fn new(system: &'a MolecularSystem, forcefield: &'a Forcefield) -> Self {
        Self { system, forcefield }
    }

    pub fn score_interaction(
        &self,
        query_atom_ids: &[AtomId],
        environment_atom_ids: &[AtomId],
    ) -> Result<EnergyTerm, ScoringError> {
        let mut energy = self.score_vdw_coulomb(query_atom_ids, environment_atom_ids)?;
        energy += self.score_hbond(query_atom_ids, environment_atom_ids)?;
        Ok(energy)
    }

    pub fn score_group_internal(&self, group_ids: &[AtomId]) -> Result<EnergyTerm, ScoringError> {
        let mut energy = self.score_vdw_coulomb(group_ids, group_ids)?;
        energy += self.score_hbond(group_ids, group_ids)?;
        Ok(energy)
    }

    fn score_vdw_coulomb(
        &self,
        group1_ids: &[AtomId],
        group2_ids: &[AtomId],
    ) -> Result<EnergyTerm, ScoringError> {
        let mut energy = EnergyTerm::default();
        let is_internal = group1_ids.as_ptr() == group2_ids.as_ptr();

        for (i, &id1) in group1_ids.iter().enumerate() {
            let atom1 = self
                .system
                .atom(id1)
                .ok_or(ScoringError::AtomNotFound(id1))?;
            let neighbors1 = self.system.get_bonded_neighbors(id1).unwrap_or(&[]);
            let neighbors1_set: HashSet<_> = neighbors1.iter().copied().collect();

            let start_index = if is_internal { i + 1 } else { 0 };
            for &id2 in &group2_ids[start_index..] {
                if neighbors1_set.contains(&id2) {
                    continue;
                }

                let mut is_1_3 = false;
                for &neighbor_id in neighbors1 {
                    if let Some(neighbors_of_neighbor) =
                        self.system.get_bonded_neighbors(neighbor_id)
                    {
                        if neighbors_of_neighbor.contains(&id2) {
                            is_1_3 = true;
                            break;
                        }
                    }
                }
                if is_1_3 {
                    continue;
                }

                let atom2 = self
                    .system
                    .atom(id2)
                    .ok_or(ScoringError::AtomNotFound(id2))?;

                energy.vdw += EnergyCalculator::calculate_vdw(atom1, atom2)?;
                energy.coulomb += EnergyCalculator::calculate_coulomb(
                    atom1,
                    atom2,
                    self.forcefield.non_bonded.globals.dielectric_constant,
                );
            }
        }
        Ok(energy)
    }

    fn score_hbond(
        &self,
        group1_ids: &[AtomId],
        group2_ids: &[AtomId],
    ) -> Result<EnergyTerm, ScoringError> {
        let mut energy = EnergyTerm::default();
        let is_internal = group1_ids.as_ptr() == group2_ids.as_ptr();

        energy += self.calculate_hbond_one_way(group1_ids, group2_ids)?;

        if !is_internal {
            energy += self.calculate_hbond_one_way(group2_ids, group1_ids)?;
        }

        Ok(energy)
    }

    fn calculate_hbond_one_way(
        &self,
        donor_group_ids: &[AtomId],
        acceptor_group_ids: &[AtomId],
    ) -> Result<EnergyTerm, ScoringError> {
        let mut hbond_energy = 0.0;

        for &h_id in donor_group_ids {
            let hydrogen = self
                .system
                .atom(h_id)
                .ok_or(ScoringError::AtomNotFound(h_id))?;

            if hydrogen.hbond_type_id == 0 {
                let donor_id = *self
                    .system
                    .get_bonded_neighbors(h_id)
                    .and_then(|n| n.first())
                    .ok_or(ScoringError::DonorNotFound(h_id))?;

                if !donor_group_ids.contains(&donor_id) {
                    continue;
                }
                let donor = self
                    .system
                    .atom(donor_id)
                    .ok_or(ScoringError::AtomNotFound(donor_id))?;

                for &a_id in acceptor_group_ids {
                    if a_id == donor_id || a_id == h_id {
                        continue;
                    }

                    let acceptor = self
                        .system
                        .atom(a_id)
                        .ok_or(ScoringError::AtomNotFound(a_id))?;

                    if acceptor.hbond_type_id > 0 {
                        let hbond_key =
                            format!("{}-{}", acceptor.force_field_type, donor.force_field_type);
                        if let Some(hbond_param) = self.forcefield.non_bonded.hbond.get(&hbond_key)
                        {
                            hbond_energy += EnergyCalculator::calculate_hbond(
                                acceptor,
                                hydrogen,
                                donor,
                                hbond_param.equilibrium_distance,
                                hbond_param.well_depth,
                            );
                        }
                    }
                }
            }
        }

        Ok(EnergyTerm {
            hbond: hbond_energy,
            ..Default::default()
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::params::{GlobalParams, HBondParam, NonBondedParams, VdwParam};
    use crate::core::models::atom::{Atom, CachedVdwParam};
    use crate::core::models::chain::ChainType;
    use crate::core::models::ids::ResidueId;
    use crate::core::models::system::MolecularSystem;
    use crate::core::models::topology::BondOrder;
    use nalgebra::Point3;
    use std::collections::HashMap;

    const TOLERANCE: f64 = 1e-9;

    fn create_generic_forcefield() -> Forcefield {
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
                equilibrium_distance: 2.8,
                well_depth: 5.0,
            },
        );
        Forcefield {
            non_bonded: NonBondedParams {
                globals,
                vdw,
                hbond,
            },
            deltas: HashMap::new(),
        }
    }

    fn create_hbond_forcefield(include_hbond_params: bool) -> Forcefield {
        let mut vdw = HashMap::new();
        vdw.insert(
            "O_R".to_string(),
            VdwParam::LennardJones {
                radius: 3.0,
                well_depth: 0.2,
            },
        );
        vdw.insert(
            "H_A".to_string(),
            VdwParam::LennardJones {
                radius: 1.0,
                well_depth: 0.01,
            },
        );
        vdw.insert(
            "O_2".to_string(),
            VdwParam::LennardJones {
                radius: 3.2,
                well_depth: 0.3,
            },
        );

        let mut hbond = HashMap::new();
        if include_hbond_params {
            hbond.insert(
                "O_2-O_R".to_string(),
                HBondParam {
                    equilibrium_distance: 2.7,
                    well_depth: 6.0,
                },
            );
        }
        Forcefield {
            non_bonded: NonBondedParams {
                globals: GlobalParams {
                    dielectric_constant: 4.0,
                    potential_function: "lennard-jones-12-6".to_string(),
                },
                vdw,
                hbond,
            },
            deltas: HashMap::new(),
        }
    }

    fn create_generic_atom(
        name: &str,
        residue_id: ResidueId,
        pos: Point3<f64>,
        ff_type: &str,
        charge: f64,
        hbond_type_id: i8,
    ) -> Atom {
        let mut atom = Atom::new(name, residue_id, pos);
        atom.force_field_type = ff_type.to_string();
        atom.partial_charge = charge;
        atom.hbond_type_id = hbond_type_id;
        atom.vdw_param = match ff_type {
            "C" => CachedVdwParam::LennardJones {
                radius: 4.0,
                well_depth: 0.1,
            },
            "N" => CachedVdwParam::LennardJones {
                radius: 3.5,
                well_depth: 0.2,
            },
            "O" => CachedVdwParam::LennardJones {
                radius: 3.2,
                well_depth: 0.3,
            },
            "H" => CachedVdwParam::LennardJones {
                radius: 1.0,
                well_depth: 0.01,
            },
            _ => CachedVdwParam::None,
        };
        atom
    }

    fn create_hbond_atom(
        name: &str,
        residue_id: ResidueId,
        pos: Point3<f64>,
        ff_type: &str,
        hbond_type_id: i8,
    ) -> Atom {
        let mut atom = Atom::new(name, residue_id, pos);
        atom.force_field_type = ff_type.to_string();
        atom.hbond_type_id = hbond_type_id;

        atom.vdw_param = match ff_type {
            "O_R" => CachedVdwParam::LennardJones {
                radius: 3.0,
                well_depth: 0.2,
            },
            "H_A" => CachedVdwParam::LennardJones {
                radius: 1.0,
                well_depth: 0.01,
            },
            "O_2" => CachedVdwParam::LennardJones {
                radius: 3.2,
                well_depth: 0.3,
            },
            _ => CachedVdwParam::None,
        };
        atom
    }

    #[test]
    fn scores_vdw_and_coulomb_for_simple_interaction() {
        let mut system = MolecularSystem::new();
        let ff = create_generic_forcefield();

        let chain_id = system.add_chain('A', ChainType::Protein);
        let res1_id = system.add_residue(chain_id, 1, "RES", None).unwrap();
        let res2_id = system.add_residue(chain_id, 2, "RES", None).unwrap();

        let query_atom = create_generic_atom("C1", res1_id, Point3::origin(), "C", 0.5, -1);
        let env_atom =
            create_generic_atom("C2", res2_id, Point3::new(4.5, 0.0, 0.0), "C", -0.5, -1);

        let query_id = system.add_atom_to_residue(res1_id, query_atom).unwrap();
        let env_id = system.add_atom_to_residue(res2_id, env_atom).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let energy = scorer.score_interaction(&[query_id], &[env_id]).unwrap();

        assert!(energy.vdw < 0.0);
        assert!(energy.coulomb < 0.0);
        assert!((energy.hbond).abs() < TOLERANCE);
    }

    #[test]
    fn ignores_1_2_and_1_3_bonded_interactions() {
        let mut system = MolecularSystem::new();
        let ff = create_generic_forcefield();

        let chain_id = system.add_chain('A', ChainType::Protein);
        let res1_id = system.add_residue(chain_id, 1, "RES", None).unwrap();
        let res2_id = system.add_residue(chain_id, 2, "RES", None).unwrap();

        let atom1 = create_generic_atom("C", res1_id, Point3::new(0.0, 0.0, 0.0), "C", 0.5, -1);
        let atom2 = create_generic_atom("N", res2_id, Point3::new(1.3, 0.0, 0.0), "N", -0.5, -1);
        let atom3 = create_generic_atom("CA", res2_id, Point3::new(1.8, 1.2, 0.0), "C", 0.1, -1);

        let id1 = system.add_atom_to_residue(res1_id, atom1).unwrap();
        let id2 = system.add_atom_to_residue(res2_id, atom2).unwrap();
        let id3 = system.add_atom_to_residue(res2_id, atom3).unwrap();

        system.add_bond(id1, id2, BondOrder::Single).unwrap();
        system.add_bond(id2, id3, BondOrder::Single).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let energy_1_2 = scorer.score_interaction(&[id1], &[id2]).unwrap();
        assert!(
            (energy_1_2.total()).abs() < TOLERANCE,
            "1-2 interactions should be ignored"
        );

        let energy_1_3 = scorer.score_interaction(&[id1], &[id3]).unwrap();
        assert!(
            (energy_1_3.total()).abs() < TOLERANCE,
            "1-3 interactions should be ignored"
        );
    }

    #[test]
    fn scores_ideal_hbond_with_correct_energy() {
        let ff = create_hbond_forcefield(true);
        let hbond_params = ff.non_bonded.hbond.get("O_2-O_R").unwrap();
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let ser_id = system.add_residue(chain_id, 1, "SER", None).unwrap();
        let asp_id = system.add_residue(chain_id, 2, "ASP", None).unwrap();

        let og = create_hbond_atom("OG", ser_id, Point3::new(0.0, 0.0, 0.0), "O_R", 1);
        let hg = create_hbond_atom("HG", ser_id, Point3::new(1.0, 0.0, 0.0), "H_A", 0);
        let od1 = create_hbond_atom(
            "OD1",
            asp_id,
            Point3::new(hbond_params.equilibrium_distance, 0.0, 0.0),
            "O_2",
            1,
        );

        let og_id = system.add_atom_to_residue(ser_id, og).unwrap();
        let hg_id = system.add_atom_to_residue(ser_id, hg).unwrap();
        system.add_atom_to_residue(asp_id, od1).unwrap();
        system.add_bond(og_id, hg_id, BondOrder::Single).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let energy = scorer
            .score_interaction(
                system.residue(ser_id).unwrap().atoms(),
                system.residue(asp_id).unwrap().atoms(),
            )
            .unwrap();

        assert!(
            (energy.hbond + hbond_params.well_depth).abs() < TOLERANCE,
            "Ideal H-bond energy should be -well_depth"
        );
    }

    #[test]
    fn hbond_energy_is_zero_for_unfavorable_angle() {
        let ff = create_hbond_forcefield(true);
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let ser_id = system.add_residue(chain_id, 1, "SER", None).unwrap();
        let asp_id = system.add_residue(chain_id, 2, "ASP", None).unwrap();

        let og = create_hbond_atom("OG", ser_id, Point3::new(0.0, 0.0, 0.0), "O_R", 1);
        let hg = create_hbond_atom("HG", ser_id, Point3::new(1.0, 0.0, 0.0), "H_A", 0);
        let od1 = create_hbond_atom("OD1", asp_id, Point3::new(1.0, 2.0, 0.0), "O_2", 1);

        let og_id = system.add_atom_to_residue(ser_id, og).unwrap();
        let hg_id = system.add_atom_to_residue(ser_id, hg).unwrap();
        system.add_atom_to_residue(asp_id, od1).unwrap();
        system.add_bond(og_id, hg_id, BondOrder::Single).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let energy = scorer
            .score_interaction(
                system.residue(ser_id).unwrap().atoms(),
                system.residue(asp_id).unwrap().atoms(),
            )
            .unwrap();

        assert!(
            (energy.hbond).abs() < TOLERANCE,
            "H-bond energy should be zero for angles <= 90 degrees"
        );
    }

    #[test]
    fn hbond_energy_is_zero_if_parameters_are_missing() {
        let ff = create_hbond_forcefield(false);
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let ser_id = system.add_residue(chain_id, 1, "SER", None).unwrap();
        let asp_id = system.add_residue(chain_id, 2, "ASP", None).unwrap();

        let og = create_hbond_atom("OG", ser_id, Point3::new(0.0, 0.0, 0.0), "O_R", 1);
        let hg = create_hbond_atom("HG", ser_id, Point3::new(1.0, 0.0, 0.0), "H_A", 0);
        let od1 = create_hbond_atom("OD1", asp_id, Point3::new(2.7, 0.0, 0.0), "O_2", 1);

        let og_id = system.add_atom_to_residue(ser_id, og).unwrap();
        let hg_id = system.add_atom_to_residue(ser_id, hg).unwrap();
        system.add_atom_to_residue(asp_id, od1).unwrap();
        system.add_bond(og_id, hg_id, BondOrder::Single).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let energy = scorer
            .score_interaction(
                system.residue(ser_id).unwrap().atoms(),
                system.residue(asp_id).unwrap().atoms(),
            )
            .unwrap();

        assert!(
            (energy.hbond).abs() < TOLERANCE,
            "H-bond energy should be zero if parameters are missing"
        );
    }

    #[test]
    fn hbond_respects_flat_bottom_potential() {
        let ff = create_hbond_forcefield(true);
        let hbond_params = ff.non_bonded.hbond.get("O_2-O_R").unwrap();
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let ser_id = system.add_residue(chain_id, 1, "SER", None).unwrap();
        let asp_id = system.add_residue(chain_id, 2, "ASP", None).unwrap();

        let delta = 0.2;
        let dist_inside_well = hbond_params.equilibrium_distance - 0.1;

        let mut og = create_hbond_atom("OG", ser_id, Point3::new(0.0, 0.0, 0.0), "O_R", 1);
        og.delta = delta;
        let hg = create_hbond_atom("HG", ser_id, Point3::new(1.0, 0.0, 0.0), "H_A", 0);
        let mut od1 = create_hbond_atom(
            "OD1",
            asp_id,
            Point3::new(dist_inside_well, 0.0, 0.0),
            "O_2",
            1,
        );
        od1.delta = delta;

        let og_id = system.add_atom_to_residue(ser_id, og).unwrap();
        let hg_id = system.add_atom_to_residue(ser_id, hg).unwrap();
        system.add_atom_to_residue(asp_id, od1).unwrap();
        system.add_bond(og_id, hg_id, BondOrder::Single).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let energy = scorer
            .score_interaction(
                system.residue(ser_id).unwrap().atoms(),
                system.residue(asp_id).unwrap().atoms(),
            )
            .unwrap();

        assert!(
            (energy.hbond + hbond_params.well_depth).abs() < TOLERANCE,
            "H-bond energy should remain at -well_depth inside the flat-bottom region"
        );
    }

    #[test]
    fn ignores_hbond_within_same_group() {
        let mut system = MolecularSystem::new();
        let ff = create_generic_forcefield();

        let chain_id = system.add_chain('A', ChainType::Protein);
        let res1_id = system.add_residue(chain_id, 1, "RES", None).unwrap();
        let res2_id = system.add_residue(chain_id, 2, "RES", None).unwrap();

        let donor = create_generic_atom("N", res1_id, Point3::new(0.0, 0.0, 0.0), "N", -0.3, 1);
        let hydrogen = create_generic_atom("H", res1_id, Point3::new(1.0, 0.0, 0.0), "H", 0.3, 0);
        let acceptor = create_generic_atom("O", res2_id, Point3::new(2.8, 0.0, 0.0), "O", -0.5, 1);

        let donor_id = system.add_atom_to_residue(res1_id, donor).unwrap();
        let h_id = system.add_atom_to_residue(res1_id, hydrogen).unwrap();
        let acceptor_id = system.add_atom_to_residue(res2_id, acceptor).unwrap();
        system.add_bond(donor_id, h_id, BondOrder::Single).unwrap();

        let res3_id = system.add_residue(chain_id, 3, "RES", None).unwrap();
        let env_atom = create_generic_atom("C", res3_id, Point3::new(5.0, 5.0, 5.0), "C", 0.1, -1);
        let env_id = system.add_atom_to_residue(res3_id, env_atom).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let energy = scorer
            .score_interaction(&[donor_id, h_id, acceptor_id], &[env_id])
            .unwrap();

        assert!((energy.hbond).abs() < TOLERANCE);
        assert_ne!(energy.vdw, 0.0);
        assert_ne!(energy.coulomb, 0.0);
    }

    #[test]
    fn returns_error_for_missing_atom() {
        let system = MolecularSystem::new();
        let ff = create_generic_forcefield();
        let scorer = Scorer::new(&system, &ff);

        let fake_id = AtomId::default();
        let result = scorer.score_interaction(&[fake_id], &[]);

        assert!(matches!(result, Err(ScoringError::AtomNotFound(_))));
    }

    #[test]
    fn returns_error_for_hydrogen_without_donor() {
        let mut system = MolecularSystem::new();
        let ff = create_generic_forcefield();

        let chain_id = system.add_chain('A', ChainType::Protein);
        let res1_id = system.add_residue(chain_id, 1, "RES", None).unwrap();

        let hydrogen = create_generic_atom("H", res1_id, Point3::new(1.0, 0.0, 0.0), "H", 0.3, 0);
        let h_id = system.add_atom_to_residue(res1_id, hydrogen).unwrap();

        let scorer = Scorer::new(&system, &ff);
        let result = scorer.score_interaction(&[h_id], &[]);

        assert!(matches!(result, Err(ScoringError::DonorNotFound(_))));
    }
}
