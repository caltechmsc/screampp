use crate::core::forcefield::params::Forcefield;
use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::engine::error::EngineError;
use crate::engine::progress::{Progress, ProgressReporter};
use itertools::Itertools;
use std::cmp::Ordering;
use std::collections::HashSet;
use tracing::{info, instrument};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Debug, Clone, PartialEq)]
pub struct ClashPair {
    pub residue_a: ResidueId,
    pub residue_b: ResidueId,
    pub energy: EnergyTerm,
}

impl Eq for ClashPair {}

impl PartialOrd for ClashPair {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        other.energy.total().partial_cmp(&self.energy.total())
    }
}

impl Ord for ClashPair {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}

#[instrument(skip_all, name = "clash_detection_task")]
pub fn run(
    system: &MolecularSystem,
    forcefield: &Forcefield,
    active_residues: &HashSet<ResidueId>,
    clash_threshold_kcal_mol: f64,
    reporter: &ProgressReporter,
) -> Result<Vec<ClashPair>, EngineError> {
    info!(
        threshold = clash_threshold_kcal_mol,
        "Detecting residue clashes."
    );

    let residue_pairs: Vec<_> = active_residues.iter().combinations(2).collect();

    if residue_pairs.is_empty() {
        return Ok(Vec::new());
    }

    reporter.report(Progress::TaskStart {
        total: residue_pairs.len() as u64,
    });

    let scorer = Scorer::new(system, forcefield);

    #[cfg(not(feature = "parallel"))]
    let iterator = residue_pairs.iter();

    #[cfg(feature = "parallel")]
    let iterator = residue_pairs.par_iter();

    let clashes: Result<Vec<ClashPair>, EngineError> = iterator
        .filter_map(|pair| {
            reporter.report(Progress::TaskIncrement { amount: 1 });
            let res_id_a = *pair[0];
            let res_id_b = *pair[1];

            let atoms_a = system.residue(res_id_a).unwrap().atoms();
            let atoms_b = system.residue(res_id_b).unwrap().atoms();

            match scorer.score_interaction(atoms_a, atoms_b) {
                Ok(energy) if energy.total() > clash_threshold_kcal_mol => Some(Ok(ClashPair {
                    residue_a: res_id_a,
                    residue_b: res_id_b,
                    energy,
                })),
                Ok(_) => None,
                Err(e) => Some(Err(EngineError::from(e))),
            }
        })
        .collect();

    let mut clashes = clashes?;

    reporter.report(Progress::TaskFinish);

    clashes.sort_unstable();

    info!(num_clashes = clashes.len(), "Clash detection complete.");

    Ok(clashes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::params::{Forcefield, GlobalParams, NonBondedParams, VdwParam};
    use crate::core::models::atom::{Atom, CachedVdwParam};
    use crate::core::models::chain::ChainType;
    use crate::core::models::residue::ResidueType;
    use nalgebra::Point3;
    use std::collections::HashMap;

    fn create_test_forcefield() -> Forcefield {
        let mut vdw = HashMap::new();
        vdw.insert(
            "C".to_string(),
            VdwParam::LennardJones {
                radius: 2.0,
                well_depth: 0.1,
            },
        );
        let non_bonded = NonBondedParams {
            globals: GlobalParams {
                dielectric_constant: 1.0,
                potential_function: "lennard-jones-12-6".to_string(),
            },
            vdw,
            hbond: HashMap::new(),
        };
        Forcefield {
            non_bonded,
            deltas: HashMap::new(),
        }
    }

    fn create_test_system() -> (MolecularSystem, HashSet<ResidueId>, ResidueId, ResidueId) {
        let mut system = MolecularSystem::new();
        let chain_a = system.add_chain('A', ChainType::Protein);

        let res1_id = system
            .add_residue(chain_a, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let mut atom1 = Atom::new("CA", res1_id, Point3::new(0.0, 0.0, 0.0));
        atom1.force_field_type = "C".to_string();
        atom1.vdw_param = CachedVdwParam::LennardJones {
            radius: 1.7,
            well_depth: 0.2,
        };
        system.add_atom_to_residue(res1_id, atom1).unwrap();

        let res2_id = system
            .add_residue(chain_a, 2, "LEU", Some(ResidueType::Leucine))
            .unwrap();
        let mut atom2 = Atom::new("CA", res2_id, Point3::new(10.0, 0.0, 0.0));
        atom2.force_field_type = "C".to_string();
        atom2.vdw_param = CachedVdwParam::LennardJones {
            radius: 1.7,
            well_depth: 0.2,
        };
        system.add_atom_to_residue(res2_id, atom2).unwrap();

        let active_residues = vec![res1_id, res2_id].into_iter().collect();

        (system, active_residues, res1_id, res2_id)
    }

    #[test]
    fn run_detects_clash_when_residues_are_close() {
        let (mut system, active_residues, _, res2_id) = create_test_system();
        let res2_atom_id = system.residue(res2_id).unwrap().atoms()[0];
        system.atom_mut(res2_atom_id).unwrap().position = Point3::new(0.1, 0.0, 0.0);

        let ff = create_test_forcefield();
        let reporter = ProgressReporter::default();

        let clashes = run(&system, &ff, &active_residues, 1.0, &reporter).unwrap();

        assert_eq!(clashes.len(), 1);
        assert!(clashes[0].energy.total() > 1.0);
    }

    #[test]
    fn run_detects_no_clash_when_residues_are_far() {
        let (system, active_residues, _, _) = create_test_system();
        let ff = create_test_forcefield();
        let reporter = ProgressReporter::default();

        let clashes = run(&system, &ff, &active_residues, 1.0, &reporter).unwrap();

        assert!(clashes.is_empty());
    }

    #[test]
    fn run_respects_clash_threshold() {
        let (mut system, active_residues, _, res2_id) = create_test_system();
        let res2_atom_id = system.residue(res2_id).unwrap().atoms()[0];
        system.atom_mut(res2_atom_id).unwrap().position = Point3::new(0.1, 0.0, 0.0);
        let ff = create_test_forcefield();
        let reporter = ProgressReporter::default();

        let clashes_low_threshold = run(&system, &ff, &active_residues, 1.0, &reporter).unwrap();
        assert_eq!(clashes_low_threshold.len(), 1);

        let clashes_high_threshold = run(&system, &ff, &active_residues, 1e15, &reporter).unwrap();
        assert!(clashes_high_threshold.is_empty());
    }

    #[test]
    fn run_returns_empty_vec_for_no_pairs() {
        let (mut system, mut active_residues, _res1_id, res2_id) = create_test_system();
        system.remove_residue(res2_id);
        active_residues.remove(&res2_id);
        assert_eq!(active_residues.len(), 1);

        let ff = create_test_forcefield();
        let reporter = ProgressReporter::default();

        let clashes = run(&system, &ff, &active_residues, 1.0, &reporter).unwrap();

        assert!(clashes.is_empty());
    }

    #[test]
    fn run_sorts_clashes_by_energy_descending() {
        let (mut system, _, res1_id, res2_id) = create_test_system();
        let chain_a = system.find_chain_by_id('A').unwrap();

        let res3_id = system
            .add_residue(chain_a, 3, "LEU", Some(ResidueType::Leucine))
            .unwrap();
        let mut atom3 = Atom::new("CA", res3_id, Point3::new(1.5, 0.0, 0.0));
        atom3.force_field_type = "C".to_string();
        atom3.vdw_param = CachedVdwParam::LennardJones {
            radius: 1.7,
            well_depth: 0.2,
        };
        system.add_atom_to_residue(res3_id, atom3).unwrap();

        let res2_atom_id = system.residue(res2_id).unwrap().atoms()[0];
        system.atom_mut(res2_atom_id).unwrap().position = Point3::new(0.1, 0.0, 0.0);

        let active_residues = vec![res1_id, res2_id, res3_id].into_iter().collect();

        let ff = create_test_forcefield();
        let reporter = ProgressReporter::default();

        let clashes = run(&system, &ff, &active_residues, 0.01, &reporter).unwrap();

        assert_eq!(
            clashes.len(),
            3,
            "Expected 3 clashes between the three residues"
        );

        assert!(
            clashes[0].energy.total() >= clashes[1].energy.total(),
            "Clash[0] should have energy >= Clash[1]"
        );
        assert!(
            clashes[1].energy.total() >= clashes[2].energy.total(),
            "Clash[1] should have energy >= Clash[2]"
        );

        let severe_clash_pair = &clashes[0];
        let severe_clash_res_ids = [severe_clash_pair.residue_a, severe_clash_pair.residue_b];
        assert!(
            severe_clash_res_ids.contains(&res1_id),
            "The worst clash should involve res1"
        );
        assert!(
            severe_clash_res_ids.contains(&res2_id),
            "The worst clash should involve res2"
        );
    }
}
