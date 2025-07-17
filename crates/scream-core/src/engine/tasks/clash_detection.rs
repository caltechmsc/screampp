use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::engine::context::{Context, ProvidesResidueSelections};
use crate::engine::error::EngineError;
use crate::engine::progress::Progress;
use itertools::Itertools;
use std::cmp::Ordering;
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
pub fn run<C: ProvidesResidueSelections + Sync>(
    system: &MolecularSystem,
    context: &Context<C>,
    clash_threshold_kcal_mol: f64,
) -> Result<Vec<ClashPair>, EngineError> {
    info!(
        threshold = clash_threshold_kcal_mol,
        "Detecting residue clashes."
    );
    context
        .reporter
        .report(Progress::Message("Detecting clashes...".to_string()));

    let active_residues = context.resolve_all_active_residues()?;
    let residue_pairs: Vec<_> = active_residues.into_iter().combinations(2).collect();

    if residue_pairs.is_empty() {
        return Ok(Vec::new());
    }

    context.reporter.report(Progress::TaskStart {
        total_steps: residue_pairs.len() as u64,
    });

    let scorer = Scorer::new(system, context.forcefield);

    #[cfg(not(feature = "parallel"))]
    let iterator = residue_pairs.iter();

    #[cfg(feature = "parallel")]
    let iterator = residue_pairs.par_iter();

    let mut clashes: Vec<ClashPair> = iterator
        .filter_map(|pair| {
            let res_id_a = pair[0];
            let res_id_b = pair[1];

            let atoms_a = system.residue(res_id_a).unwrap().atoms().to_vec();
            let atoms_b = system.residue(res_id_b).unwrap().atoms().to_vec();

            match scorer.score_interaction(&atoms_a, &atoms_b) {
                Ok(energy) if energy.total() > clash_threshold_kcal_mol => {
                    context.reporter.report(Progress::TaskIncrement);
                    Some(ClashPair {
                        residue_a: res_id_a,
                        residue_b: res_id_b,
                        energy,
                    })
                }
                _ => {
                    context.reporter.report(Progress::TaskIncrement);
                    None
                }
            }
        })
        .collect();

    context.reporter.report(Progress::TaskFinish);

    clashes.sort_unstable_by(|a, b| a.cmp(b));

    info!(num_clashes = clashes.len(), "Clash detection complete.");

    Ok(clashes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::params::{Forcefield, GlobalParams, NonBondedParams, VdwParam};
    use crate::core::models::atom::Atom;
    use crate::core::models::chain::ChainType;
    use crate::core::models::residue::ResidueType;
    use crate::core::rotamers::library::RotamerLibrary;
    use crate::engine::config::{PlacementConfig, ResidueSelection};
    use crate::engine::progress::ProgressReporter;
    use nalgebra::Point3;
    use std::collections::HashMap;

    fn create_test_forcefield() -> Forcefield {
        let globals = GlobalParams {
            dielectric_constant: 1.0,
            potential_function: "lj".to_string(),
        };
        let mut vdw = HashMap::new();
        vdw.insert(
            "C".to_string(),
            VdwParam::LennardJones {
                radius: 2.0,
                well_depth: 0.1,
            },
        );
        let non_bonded = NonBondedParams {
            globals,
            vdw,
            hbond: HashMap::new(),
        };
        Forcefield {
            non_bonded,
            deltas: HashMap::new(),
        }
    }

    fn create_test_rotamer_library() -> RotamerLibrary {
        let mut lib = RotamerLibrary::default();
        lib.rotamers
            .insert(ResidueType::Alanine, Default::default());
        lib.rotamers
            .insert(ResidueType::Leucine, Default::default());
        lib
    }

    fn create_test_placement_config(selection: ResidueSelection) -> PlacementConfig {
        PlacementConfig {
            residues_to_optimize: selection,
            scoring: crate::engine::config::ScoringConfig {
                forcefield_path: "".into(),
                rotamer_library_path: "".into(),
                delta_params_path: "".into(),
                s_factor: 0.0,
            },
            optimization: crate::engine::config::OptimizationConfig {
                max_iterations: 1,
                convergence_threshold: 0.1,
                num_solutions: 1,
                include_input_conformation: false,
            },
        }
    }

    fn create_clashing_system() -> MolecularSystem {
        let mut system = MolecularSystem::new();
        let chain_a = system.add_chain('A', ChainType::Protein);

        let res1_id = system
            .add_residue(chain_a, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let mut atom1 = Atom::new(1, "CA", res1_id, Point3::new(0.0, 0.0, 0.0));
        atom1.force_field_type = "C".to_string();
        system.add_atom_to_residue(res1_id, atom1).unwrap();

        let res2_id = system
            .add_residue(chain_a, 2, "LEU", Some(ResidueType::Leucine))
            .unwrap();
        let mut atom2 = Atom::new(2, "CA", res2_id, Point3::new(0.1, 0.0, 0.0));
        atom2.force_field_type = "C".to_string();
        system.add_atom_to_residue(res2_id, atom2).unwrap();

        system
    }

    fn create_non_clashing_system() -> MolecularSystem {
        let mut system = MolecularSystem::new();
        let chain_a = system.add_chain('A', ChainType::Protein);

        let res1_id = system
            .add_residue(chain_a, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let mut atom1 = Atom::new(1, "CA", res1_id, Point3::new(0.0, 0.0, 0.0));
        atom1.force_field_type = "C".to_string();
        system.add_atom_to_residue(res1_id, atom1).unwrap();

        let res2_id = system
            .add_residue(chain_a, 2, "LEU", Some(ResidueType::Leucine))
            .unwrap();
        let mut atom2 = Atom::new(2, "CA", res2_id, Point3::new(10.0, 0.0, 0.0));
        atom2.force_field_type = "C".to_string();
        system.add_atom_to_residue(res2_id, atom2).unwrap();

        system
    }

    #[test]
    fn run_detects_clash_when_residues_are_close() {
        let system = create_clashing_system();
        let ff = create_test_forcefield();
        let rot_lib = create_test_rotamer_library();
        let reporter = ProgressReporter::default();
        let config = create_test_placement_config(ResidueSelection::All);
        let context = Context::new(&system, &config, &ff, &rot_lib, &reporter);

        let clashes = run(&system, &context, 1.0).unwrap();

        assert_eq!(clashes.len(), 1);
        assert!(clashes[0].energy.total() > 1.0);
    }

    #[test]
    fn run_detects_no_clash_when_residues_are_far() {
        let system = create_non_clashing_system();
        let ff = create_test_forcefield();
        let rot_lib = create_test_rotamer_library();
        let reporter = ProgressReporter::default();
        let config = create_test_placement_config(ResidueSelection::All);
        let context = Context::new(&system, &config, &ff, &rot_lib, &reporter);

        let clashes = run(&system, &context, 1.0).unwrap();

        assert!(clashes.is_empty());
    }

    #[test]
    fn run_respects_clash_threshold() {
        let system = create_clashing_system();
        let ff = create_test_forcefield();
        let rot_lib = create_test_rotamer_library();
        let reporter = ProgressReporter::default();
        let config = create_test_placement_config(ResidueSelection::All);
        let context = Context::new(&system, &config, &ff, &rot_lib, &reporter);

        let clashes_low_threshold = run(&system, &context, 1.0).unwrap();
        assert_eq!(clashes_low_threshold.len(), 1);

        let clashes_high_threshold = run(&system, &context, 1e15).unwrap();
        assert!(
            clashes_high_threshold.is_empty(),
            "Expected no clashes with a very high threshold"
        );
    }

    #[test]
    fn run_returns_empty_vec_for_no_pairs() {
        let mut system = create_non_clashing_system();
        let chain_a = system.find_chain_by_id('A').unwrap();
        let res2_id = system.find_residue_by_id(chain_a, 2).unwrap();
        system.remove_residue(res2_id);

        let ff = create_test_forcefield();
        let rot_lib = create_test_rotamer_library();
        let reporter = ProgressReporter::default();
        let config = create_test_placement_config(ResidueSelection::All);
        let context = Context::new(&system, &config, &ff, &rot_lib, &reporter);

        let clashes = run(&system, &context, 1.0).unwrap();

        assert!(clashes.is_empty());
    }

    #[test]
    fn run_sorts_clashes_by_energy_descending() {
        let mut system = MolecularSystem::new();
        let chain_a = system.add_chain('A', ChainType::Protein);

        let res1_id = system
            .add_residue(chain_a, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let mut atom1 = Atom::new(1, "CA", res1_id, Point3::new(0.0, 0.0, 0.0));
        atom1.force_field_type = "C".to_string();
        system.add_atom_to_residue(res1_id, atom1).unwrap();

        let res2_id = system
            .add_residue(chain_a, 2, "LEU", Some(ResidueType::Leucine))
            .unwrap();
        let mut atom2 = Atom::new(2, "CA", res2_id, Point3::new(1.5, 0.0, 0.0));
        atom2.force_field_type = "C".to_string();
        system.add_atom_to_residue(res2_id, atom2).unwrap();

        let res3_id = system
            .add_residue(chain_a, 3, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let mut atom3 = Atom::new(3, "CA", res3_id, Point3::new(0.1, 0.0, 0.0));
        atom3.force_field_type = "C".to_string();
        system.add_atom_to_residue(res3_id, atom3).unwrap();

        let ff = create_test_forcefield();
        let rot_lib = create_test_rotamer_library();
        let reporter = ProgressReporter::default();
        let config = create_test_placement_config(ResidueSelection::All);
        let context = Context::new(&system, &config, &ff, &rot_lib, &reporter);

        let clashes = run(&system, &context, 0.01).unwrap();

        assert!(clashes.len() >= 2);
        assert!(clashes[0].energy.total() > clashes[1].energy.total());

        let severe_clash_pair = &clashes[0];
        let severe_clash_res_ids = [severe_clash_pair.residue_a, severe_clash_pair.residue_b];
        assert!(severe_clash_res_ids.contains(&res1_id));
        assert!(severe_clash_res_ids.contains(&res3_id));
    }
}
