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
                            format!("{}-{}", donor.force_field_type, acceptor.force_field_type);

                        if let Some(hbond_param) = self.forcefield.non_bonded.hbond.get(&hbond_key)
                        {
                            let energy_contribution = EnergyCalculator::calculate_hbond(
                                donor,
                                hydrogen,
                                acceptor,
                                hbond_param.equilibrium_distance,
                                hbond_param.well_depth,
                            );

                            hbond_energy += energy_contribution;
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
    use crate::core::forcefield::params::{GlobalParams, HBondParam, NonBondedParams};
    use crate::core::models::atom::{Atom, AtomRole, CachedVdwParam};
    use crate::core::models::chain::ChainType;
    use crate::core::models::ids::ResidueId;
    use crate::core::models::system::MolecularSystem;
    use crate::core::models::topology::BondOrder;
    use nalgebra::Point3;
    use std::collections::HashMap;

    const TOLERANCE: f64 = 1e-9;

    struct TestSetup {
        system: MolecularSystem,
        forcefield: Forcefield,
        res_ala_id: ResidueId,
        res_gly_id: ResidueId,
        res_ser_id: ResidueId,
        res_leu_id: ResidueId,
    }

    impl TestSetup {
        fn new() -> Self {
            let mut system = MolecularSystem::new();
            let forcefield = Self::create_test_forcefield();

            let chain_a = system.add_chain('A', ChainType::Protein);

            let res_ala_id = system.add_residue(chain_a, 1, "ALA", None).unwrap();
            let ala_atoms = vec![
                (
                    "N",
                    Point3::new(-1.2, 0.5, 0.0),
                    "N_R",
                    -0.3,
                    AtomRole::Backbone,
                    -1,
                ),
                (
                    "CA",
                    Point3::new(0.0, 0.0, 0.0),
                    "C_3",
                    0.1,
                    AtomRole::Backbone,
                    -1,
                ),
                (
                    "C",
                    Point3::new(0.8, -0.8, 0.0),
                    "C_R",
                    0.5,
                    AtomRole::Backbone,
                    -1,
                ),
                (
                    "O",
                    Point3::new(1.5, -1.2, 0.0),
                    "O_2",
                    -0.5,
                    AtomRole::Backbone,
                    1,
                ),
                (
                    "CB",
                    Point3::new(-0.5, -1.0, 1.2),
                    "C_M",
                    -0.1,
                    AtomRole::Sidechain,
                    -1,
                ),
            ];
            Self::add_atoms_to_residue(&mut system, res_ala_id, &ala_atoms);
            system
                .add_bond(
                    system
                        .residue(res_ala_id)
                        .unwrap()
                        .get_first_atom_id_by_name("N")
                        .unwrap(),
                    system
                        .residue(res_ala_id)
                        .unwrap()
                        .get_first_atom_id_by_name("CA")
                        .unwrap(),
                    BondOrder::Single,
                )
                .unwrap();
            system
                .add_bond(
                    system
                        .residue(res_ala_id)
                        .unwrap()
                        .get_first_atom_id_by_name("CA")
                        .unwrap(),
                    system
                        .residue(res_ala_id)
                        .unwrap()
                        .get_first_atom_id_by_name("C")
                        .unwrap(),
                    BondOrder::Single,
                )
                .unwrap();
            system
                .add_bond(
                    system
                        .residue(res_ala_id)
                        .unwrap()
                        .get_first_atom_id_by_name("C")
                        .unwrap(),
                    system
                        .residue(res_ala_id)
                        .unwrap()
                        .get_first_atom_id_by_name("O")
                        .unwrap(),
                    BondOrder::Single,
                )
                .unwrap();
            system
                .add_bond(
                    system
                        .residue(res_ala_id)
                        .unwrap()
                        .get_first_atom_id_by_name("CA")
                        .unwrap(),
                    system
                        .residue(res_ala_id)
                        .unwrap()
                        .get_first_atom_id_by_name("CB")
                        .unwrap(),
                    BondOrder::Single,
                )
                .unwrap();

            let res_gly_id = system.add_residue(chain_a, 2, "GLY", None).unwrap();

            let res_ser_id = system.add_residue(chain_a, 3, "SER", None).unwrap();
            let ser_atoms = vec![
                (
                    "OG",
                    Point3::new(5.0, 5.0, 0.0),
                    "O_H",
                    -0.6,
                    AtomRole::Sidechain,
                    1,
                ),
                (
                    "HG",
                    Point3::new(5.0, 5.8, 0.0),
                    "H_O",
                    0.4,
                    AtomRole::Sidechain,
                    0,
                ),
            ];
            Self::add_atoms_to_residue(&mut system, res_ser_id, &ser_atoms);
            system
                .add_bond(
                    system
                        .residue(res_ser_id)
                        .unwrap()
                        .get_first_atom_id_by_name("OG")
                        .unwrap(),
                    system
                        .residue(res_ser_id)
                        .unwrap()
                        .get_first_atom_id_by_name("HG")
                        .unwrap(),
                    BondOrder::Single,
                )
                .unwrap();

            let res_leu_id = system.add_residue(chain_a, 4, "LEU", None).unwrap();
            let leu_atoms = vec![(
                "CD1",
                Point3::new(0.0, -1.0, 5.0),
                "C_M",
                -0.1,
                AtomRole::Sidechain,
                -1,
            )];
            Self::add_atoms_to_residue(&mut system, res_leu_id, &leu_atoms);

            Self {
                system,
                forcefield,
                res_ala_id,
                res_gly_id,
                res_ser_id,
                res_leu_id,
            }
        }

        fn add_atoms_to_residue(
            system: &mut MolecularSystem,
            res_id: ResidueId,
            atom_data: &[(&str, Point3<f64>, &str, f64, AtomRole, i8)],
        ) {
            for (name, pos, ff_type, charge, role, hbond_id) in atom_data {
                let mut atom = Atom::new(name, res_id, *pos);
                atom.force_field_type = ff_type.to_string();
                atom.partial_charge = *charge;
                atom.role = *role;
                atom.hbond_type_id = *hbond_id;
                atom.vdw_param = CachedVdwParam::LennardJones {
                    radius: 3.5,
                    well_depth: 0.1,
                };
                system.add_atom_to_residue(res_id, atom).unwrap();
            }
        }

        fn create_test_forcefield() -> Forcefield {
            let mut hbond = HashMap::new();
            hbond.insert(
                "O_2-O_H".to_string(),
                HBondParam {
                    equilibrium_distance: 2.8,
                    well_depth: 5.0,
                },
            );
            Forcefield {
                non_bonded: NonBondedParams {
                    globals: GlobalParams {
                        dielectric_constant: 4.0,
                        potential_function: "lj-12-6".to_string(),
                    },
                    vdw: HashMap::new(),
                    hbond,
                    hbond_donors: HashSet::new(),
                    hbond_acceptors: HashSet::new(),
                },
                deltas: HashMap::new(),
            }
        }
    }

    mod score_group_internal_tests {
        use super::*;

        #[test]
        fn calculates_non_zero_energy_for_a_residue() {
            let setup = TestSetup::new();
            let scorer = Scorer::new(&setup.system, &setup.forcefield);
            let ala_atoms = setup.system.residue(setup.res_ala_id).unwrap().atoms();

            let energy = scorer.score_group_internal(ala_atoms).unwrap();

            assert!(
                energy.total().abs() > 1e-6,
                "Internal energy should not be zero"
            );
        }

        #[test]
        fn excludes_1_2_and_1_3_interactions_correctly() {
            let setup = TestSetup::new();
            let scorer = Scorer::new(&setup.system, &setup.forcefield);

            let n_id = setup
                .system
                .residue(setup.res_ala_id)
                .unwrap()
                .get_first_atom_id_by_name("N")
                .unwrap();
            let ca_id = setup
                .system
                .residue(setup.res_ala_id)
                .unwrap()
                .get_first_atom_id_by_name("CA")
                .unwrap();
            let energy_1_2 = scorer.score_group_internal(&[n_id, ca_id]).unwrap();
            assert_eq!(energy_1_2.total(), 0.0, "1-2 interaction should be zero");

            let c_id = setup
                .system
                .residue(setup.res_ala_id)
                .unwrap()
                .get_first_atom_id_by_name("C")
                .unwrap();
            let energy_1_3 = scorer.score_group_internal(&[n_id, c_id]).unwrap();
            assert_eq!(energy_1_3.total(), 0.0, "1-3 interaction should be zero");
        }

        #[test]
        fn includes_1_4_interactions() {
            let setup = TestSetup::new();
            let scorer = Scorer::new(&setup.system, &setup.forcefield);

            let n_id = setup
                .system
                .residue(setup.res_ala_id)
                .unwrap()
                .get_first_atom_id_by_name("N")
                .unwrap();
            let o_id = setup
                .system
                .residue(setup.res_ala_id)
                .unwrap()
                .get_first_atom_id_by_name("O")
                .unwrap();

            let energy_1_4 = scorer.score_group_internal(&[n_id, o_id]).unwrap();
            assert!(
                energy_1_4.total().abs() > 1e-6,
                "1-4 interaction (N-O) should be non-zero"
            );
        }

        #[test]
        fn internal_hbond_is_not_counted() {
            let setup = TestSetup::new();
            let scorer = Scorer::new(&setup.system, &setup.forcefield);
            let ser_atoms = setup.system.residue(setup.res_ser_id).unwrap().atoms();

            let energy = scorer.score_group_internal(ser_atoms).unwrap();
            assert!(
                (energy.hbond).abs() < TOLERANCE,
                "Internal H-bond should not be counted in score_group_internal"
            );
        }
    }

    mod score_interaction_tests {
        use super::*;
        use nalgebra::Vector3;

        #[test]
        fn calculates_attractive_coulomb_for_distant_residues() {
            let mut setup = TestSetup::new();

            let ala_c_id = setup
                .system
                .residue(setup.res_ala_id)
                .unwrap()
                .get_first_atom_id_by_name("C")
                .unwrap();
            let leu_cd1_id = setup
                .system
                .residue(setup.res_leu_id)
                .unwrap()
                .get_first_atom_id_by_name("CD1")
                .unwrap();

            setup.system.atom_mut(leu_cd1_id).unwrap().partial_charge = -0.5;

            let scorer = Scorer::new(&setup.system, &setup.forcefield);
            let energy = scorer
                .score_interaction(&[ala_c_id], &[leu_cd1_id])
                .unwrap();
            assert!(
                energy.coulomb < 0.0,
                "Coulomb energy for a single (+,-) pair should be attractive"
            );
        }

        #[test]
        fn calculates_repulsive_vdw_for_clashing_residues() {
            let mut setup = TestSetup::new();
            let leu_cd1_id = setup
                .system
                .residue(setup.res_leu_id)
                .unwrap()
                .get_first_atom_id_by_name("CD1")
                .unwrap();
            let ala_cb_id = setup
                .system
                .residue(setup.res_ala_id)
                .unwrap()
                .get_first_atom_id_by_name("CB")
                .unwrap();
            let ala_cb_pos = setup.system.atom(ala_cb_id).unwrap().position;
            setup.system.atom_mut(leu_cd1_id).unwrap().position =
                ala_cb_pos + Vector3::new(0.1, 0.0, 0.0);

            let scorer = Scorer::new(&setup.system, &setup.forcefield);
            let ala_atoms = setup.system.residue(setup.res_ala_id).unwrap().atoms();
            let leu_atoms = setup.system.residue(setup.res_leu_id).unwrap().atoms();

            let energy = scorer.score_interaction(ala_atoms, leu_atoms).unwrap();
            assert!(
                energy.vdw > 10.0,
                "VDW energy for clashing atoms should be highly repulsive"
            );
        }

        #[test]
        fn calculates_inter_residue_hbond_correctly() {
            let mut setup = TestSetup::new();
            let ala_o_id = setup
                .system
                .residue(setup.res_ala_id)
                .unwrap()
                .get_first_atom_id_by_name("O")
                .unwrap();

            let ser_og_id = setup
                .system
                .residue(setup.res_ser_id)
                .unwrap()
                .get_first_atom_id_by_name("OG")
                .unwrap();
            let ser_hg_id = setup
                .system
                .residue(setup.res_ser_id)
                .unwrap()
                .get_first_atom_id_by_name("HG")
                .unwrap();
            let donor_pos = setup.system.atom(ser_og_id).unwrap().position;
            let hydrogen_pos = setup.system.atom(ser_hg_id).unwrap().position;

            let dh_vector = hydrogen_pos - donor_pos;
            let ideal_dist = setup
                .forcefield
                .non_bonded
                .hbond
                .get("O_2-O_H")
                .unwrap()
                .equilibrium_distance;

            let acceptor_pos = donor_pos + dh_vector.normalize() * ideal_dist;

            setup.system.atom_mut(ala_o_id).unwrap().position = acceptor_pos;

            let scorer = Scorer::new(&setup.system, &setup.forcefield);
            let ala_atoms = setup.system.residue(setup.res_ala_id).unwrap().atoms();
            let ser_atoms = setup.system.residue(setup.res_ser_id).unwrap().atoms();

            let energy = scorer.score_interaction(ala_atoms, ser_atoms).unwrap();
            let expected_hbond_energy = -setup
                .forcefield
                .non_bonded
                .hbond
                .get("O_2-O_H")
                .unwrap()
                .well_depth;

            assert!(
                (energy.hbond - expected_hbond_energy).abs() < 0.1,
                "H-bond energy is incorrect. Got {}, expected ~{}",
                energy.hbond,
                expected_hbond_energy
            );
        }

        #[test]
        fn hbond_is_symmetric() {
            let mut setup = TestSetup::new();
            let ala_o_id = setup
                .system
                .residue(setup.res_ala_id)
                .unwrap()
                .get_first_atom_id_by_name("O")
                .unwrap();

            let ser_og_id = setup
                .system
                .residue(setup.res_ser_id)
                .unwrap()
                .get_first_atom_id_by_name("OG")
                .unwrap();
            let ser_hg_id = setup
                .system
                .residue(setup.res_ser_id)
                .unwrap()
                .get_first_atom_id_by_name("HG")
                .unwrap();

            let donor_pos = setup.system.atom(ser_og_id).unwrap().position;
            let hydrogen_pos = setup.system.atom(ser_hg_id).unwrap().position;
            let dh_vector = hydrogen_pos - donor_pos;
            let ideal_dist = 2.8;
            let acceptor_pos = donor_pos + dh_vector.normalize() * ideal_dist;
            setup.system.atom_mut(ala_o_id).unwrap().position = acceptor_pos;

            let scorer = Scorer::new(&setup.system, &setup.forcefield);
            let ala_atoms = setup.system.residue(setup.res_ala_id).unwrap().atoms();
            let ser_atoms = setup.system.residue(setup.res_ser_id).unwrap().atoms();

            let energy_ala_ser = scorer
                .score_interaction(ala_atoms, ser_atoms)
                .unwrap()
                .hbond;
            let energy_ser_ala = scorer
                .score_interaction(ser_atoms, ala_atoms)
                .unwrap()
                .hbond;

            assert!(
                (energy_ala_ser - energy_ser_ala).abs() < TOLERANCE,
                "H-bond calculation should be symmetric"
            );
            assert!(
                energy_ala_ser < -4.0,
                "H-bond energy should be significantly attractive"
            );
        }

        #[test]
        fn interaction_is_symmetric_for_vdw_coulomb() {
            let setup = TestSetup::new();
            let scorer = Scorer::new(&setup.system, &setup.forcefield);
            let ala_atoms = setup.system.residue(setup.res_ala_id).unwrap().atoms();
            let leu_atoms = setup.system.residue(setup.res_leu_id).unwrap().atoms();

            let energy1 = scorer.score_interaction(ala_atoms, leu_atoms).unwrap();
            let energy2 = scorer.score_interaction(leu_atoms, ala_atoms).unwrap();

            assert!((energy1.vdw - energy2.vdw).abs() < TOLERANCE);
            assert!((energy1.coulomb - energy2.coulomb).abs() < TOLERANCE);
        }

        #[test]
        fn interaction_with_empty_group_is_zero() {
            let setup = TestSetup::new();
            let scorer = Scorer::new(&setup.system, &setup.forcefield);
            let ala_atoms = setup.system.residue(setup.res_ala_id).unwrap().atoms();

            let energy = scorer.score_interaction(ala_atoms, &[]).unwrap();
            assert_eq!(energy.total(), 0.0);
        }
    }

    mod error_handling {
        use super::*;

        #[test]
        fn returns_atom_not_found_error() {
            let setup = TestSetup::new();
            let scorer = Scorer::new(&setup.system, &setup.forcefield);
            let fake_id = AtomId::default();

            let result = scorer.score_interaction(&[fake_id], &[]);
            assert!(matches!(result, Err(ScoringError::AtomNotFound(_))));

            let result_internal = scorer.score_group_internal(&[fake_id]);
            assert!(matches!(
                result_internal,
                Err(ScoringError::AtomNotFound(_))
            ));
        }

        #[test]
        fn returns_donor_not_found_error_for_hbond() {
            let mut setup = TestSetup::new();
            let ser_hg_id = setup
                .system
                .residue(setup.res_ser_id)
                .unwrap()
                .get_first_atom_id_by_name("HG")
                .unwrap();
            let ser_og_id = setup
                .system
                .residue(setup.res_ser_id)
                .unwrap()
                .get_first_atom_id_by_name("OG")
                .unwrap();
            setup.system.remove_atom(ser_og_id);

            let scorer = Scorer::new(&setup.system, &setup.forcefield);
            let ala_atoms = setup.system.residue(setup.res_ala_id).unwrap().atoms();

            let result = scorer.score_hbond(&[ser_hg_id], ala_atoms);
            assert!(matches!(result, Err(ScoringError::DonorNotFound(_))));
        }
    }
}
