use crate::core::forcefield::params::Forcefield;
use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::atom::AtomRole;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::core::models::system::MolecularSystem;
use crate::engine::error::EngineError;
use itertools::Itertools;
use std::collections::{HashMap, HashSet};

#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::sync::Mutex;
use tracing::instrument;

#[instrument(skip_all, name = "interaction_energy_task")]
pub fn run(
    system: &MolecularSystem,
    forcefield: &Forcefield,
    active_residues: &HashSet<ResidueId>,
) -> Result<EnergyTerm, EngineError> {
    if active_residues.len() < 2 {
        return Ok(EnergyTerm::default());
    }

    let sidechain_atoms_map = collect_active_sidechain_atoms(system, active_residues);

    let active_residue_vec: Vec<_> = active_residues.iter().collect();
    let residue_pairs = active_residue_vec.into_iter().combinations(2);

    #[cfg(not(feature = "parallel"))]
    let iterator = residue_pairs;

    #[cfg(feature = "parallel")]
    let iterator = residue_pairs.par_bridge();

    #[cfg(not(feature = "parallel"))]
    let total_interaction_energy = iterator
        .map(|pair| {
            let res_a_id = *pair[0];
            let res_b_id = *pair[1];

            let atoms_a = sidechain_atoms_map
                .get(&res_a_id)
                .map_or([].as_slice(), |v| v.as_slice());
            let atoms_b = sidechain_atoms_map
                .get(&res_b_id)
                .map_or([].as_slice(), |v| v.as_slice());

            if atoms_a.is_empty() || atoms_b.is_empty() {
                return Ok(EnergyTerm::default());
            }

            let scorer = Scorer::new(system, forcefield);
            scorer.score_interaction(atoms_a, atoms_b)
        })
        .try_fold(
            EnergyTerm::default(),
            |acc, term_res| -> Result<EnergyTerm, EngineError> {
                let term = term_res?;
                Ok(acc + term)
            },
        )?;

    #[cfg(feature = "parallel")]
    let total_interaction_energy = {
        let acc = Mutex::new(EnergyTerm::default());

        iterator
            .map(|pair| {
                let res_a_id = *pair[0];
                let res_b_id = *pair[1];

                let atoms_a = sidechain_atoms_map
                    .get(&res_a_id)
                    .map_or([].as_slice(), |v| v.as_slice());
                let atoms_b = sidechain_atoms_map
                    .get(&res_b_id)
                    .map_or([].as_slice(), |v| v.as_slice());

                if atoms_a.is_empty() || atoms_b.is_empty() {
                    return Ok(EnergyTerm::default());
                }

                let scorer = Scorer::new(system, forcefield);
                scorer.score_interaction(atoms_a, atoms_b)
            })
            .try_for_each(|term_res| -> Result<(), EngineError> {
                let term = term_res?;
                let mut guard = acc.lock().expect("mutex poisoned");
                *guard = *guard + term;
                Ok(())
            })?;

        acc.into_inner().expect("mutex poisoned")
    };

    Ok(total_interaction_energy)
}

fn collect_active_sidechain_atoms(
    system: &MolecularSystem,
    active_residues: &HashSet<ResidueId>,
) -> HashMap<ResidueId, Vec<AtomId>> {
    let mut map = HashMap::with_capacity(active_residues.len());
    for &residue_id in active_residues {
        if let Some(residue) = system.residue(residue_id) {
            let sidechain_ids: Vec<AtomId> = residue
                .atoms()
                .iter()
                .filter_map(|&atom_id| {
                    system.atom(atom_id).and_then(|atom| {
                        if atom.role == AtomRole::Sidechain {
                            Some(atom_id)
                        } else {
                            None
                        }
                    })
                })
                .collect();
            map.insert(residue_id, sidechain_ids);
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{
        forcefield::parameterization::Parameterizer,
        models::{
            atom::{Atom, AtomRole},
            chain::ChainType,
            residue::ResidueType,
        },
        topology::registry::TopologyRegistry,
    };
    use nalgebra::{Point3, Vector3};
    use std::fs::File;
    use std::io::Write;
    use tempfile::TempDir;

    const TOLERANCE: f64 = 1e-9;

    struct TestSetup {
        system: MolecularSystem,
        forcefield: Forcefield,
        ala_id: ResidueId,
        gly_id: ResidueId,
        leu_id: ResidueId,
        _temp_dir: TempDir,
    }

    fn setup() -> TestSetup {
        let temp_dir = TempDir::new().unwrap();
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);

        let ala_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let gly_id = system
            .add_residue(chain_id, 2, "GLY", Some(ResidueType::Glycine))
            .unwrap();
        let leu_id = system
            .add_residue(chain_id, 3, "LEU", Some(ResidueType::Leucine))
            .unwrap();

        add_atoms_to_residue(
            &mut system,
            ala_id,
            &[
                ("N", [-1.2, 0.5, 0.0]),
                ("CA", [0.0, 0.0, 0.0]),
                ("C", [0.8, -0.8, 0.0]),
                ("CB", [0.5, 1.2, 0.0]),
            ],
        );
        add_atoms_to_residue(
            &mut system,
            gly_id,
            &[
                ("N", [10.0, 1.0, 0.0]),
                ("CA", [11.2, 0.5, 0.0]),
                ("C", [12.0, 1.2, 0.0]),
            ],
        );
        add_atoms_to_residue(
            &mut system,
            leu_id,
            &[
                ("N", [20.0, 0.0, 0.0]),
                ("CA", [21.2, 0.0, 0.0]),
                ("C", [22.0, 0.0, 0.0]),
                ("CB", [21.5, 1.5, 0.0]),
                ("CG", [21.8, 2.5, 1.0]),
            ],
        );

        let forcefield = create_test_forcefield(&temp_dir);
        let topology_registry = create_test_topology_registry(&temp_dir);

        let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 0.0);
        parameterizer.parameterize_system(&mut system).unwrap();

        TestSetup {
            system,
            forcefield,
            ala_id,
            gly_id,
            leu_id,
            _temp_dir: temp_dir,
        }
    }

    fn add_atoms_to_residue(
        system: &mut MolecularSystem,
        res_id: ResidueId,
        atom_data: &[(&str, [f64; 3])],
    ) {
        for (name, pos) in atom_data {
            let mut atom = Atom::new(name, res_id, Point3::from_slice(pos));
            atom.force_field_type = if name.len() == 1 {
                "C_BB".to_string()
            } else {
                "C_SC".to_string()
            };
            system.add_atom_to_residue(res_id, atom).unwrap();
        }
    }

    fn create_test_forcefield(temp_dir: &TempDir) -> Forcefield {
        let ff_path = temp_dir.path().join("test.ff.toml");
        let mut ff_file = File::create(&ff_path).unwrap();
        write!(
            ff_file,
            r#"
            [globals]
            dielectric_constant = 4.0
            potential_function = "lennard-jones-12-6"
            [vdw]
            C_BB = {{ radius = 3.8, well_depth = 0.1 }}
            C_SC = {{ radius = 4.0, well_depth = 0.12 }}
            [hbond]
            "#
        )
        .unwrap();

        let delta_path = temp_dir.path().join("test.delta.csv");
        File::create(&delta_path).unwrap();

        Forcefield::load(&ff_path, &delta_path).unwrap()
    }

    fn create_test_topology_registry(temp_dir: &TempDir) -> TopologyRegistry {
        let topo_path = temp_dir.path().join("test.topo.toml");
        let mut topo_file = File::create(&topo_path).unwrap();
        write!(
            topo_file,
            r#"
            [ALA]
            anchor_atoms = ["N", "CA", "C"]
            sidechain_atoms = ["CB"]

            [GLY]
            anchor_atoms = ["N", "CA", "C"]
            sidechain_atoms = []

            [LEU]
            anchor_atoms = ["N", "CA", "C"]
            sidechain_atoms = ["CB", "CG"]
            "#
        )
        .unwrap();

        TopologyRegistry::load(&topo_path).unwrap()
    }

    #[test]
    fn run_returns_zero_for_no_active_residues() {
        let setup = setup();

        let active_residues = HashSet::new();
        let energy = run(&setup.system, &setup.forcefield, &active_residues).unwrap();
        assert_eq!(energy.total(), 0.0);
    }

    #[test]
    fn run_returns_zero_for_one_active_residue() {
        let setup = setup();

        let mut active_residues = HashSet::new();
        active_residues.insert(setup.ala_id);
        let energy = run(&setup.system, &setup.forcefield, &active_residues).unwrap();
        assert_eq!(energy.total(), 0.0);
    }

    #[test]
    fn run_excludes_interactions_with_residue_without_sidechain() {
        let setup = setup();

        let mut active_residues = HashSet::new();
        active_residues.insert(setup.ala_id);
        active_residues.insert(setup.gly_id);
        let energy = run(&setup.system, &setup.forcefield, &active_residues).unwrap();
        assert_eq!(energy.total(), 0.0, "Interaction with GLY should be zero");
    }

    #[test]
    fn run_calculates_interaction_between_sidechains_only() {
        let mut setup = setup();

        let ala_c_id = setup
            .system
            .residue(setup.ala_id)
            .unwrap()
            .get_first_atom_id_by_name("C")
            .unwrap();
        let leu_n_id = setup
            .system
            .residue(setup.leu_id)
            .unwrap()
            .get_first_atom_id_by_name("N")
            .unwrap();
        let leu_n_pos = setup.system.atom(leu_n_id).unwrap().position;
        setup.system.atom_mut(ala_c_id).unwrap().position = leu_n_pos + Vector3::new(0.1, 0.0, 0.0);

        let ala_cb_id = setup
            .system
            .residue(setup.ala_id)
            .unwrap()
            .get_first_atom_id_by_name("CB")
            .unwrap();
        setup.system.atom_mut(ala_cb_id).unwrap().position = Point3::new(0.0, 100.0, 0.0);

        let active_residues: HashSet<_> = [setup.ala_id, setup.leu_id].iter().cloned().collect();
        let energy = run(&setup.system, &setup.forcefield, &active_residues).unwrap();

        assert!(
            energy.total() < 0.1,
            "Energy should be near zero as only sidechains interact, but got {}",
            energy.total()
        );
    }

    #[test]
    fn run_calculates_correct_energy_for_a_single_pair() {
        let mut setup = setup();

        let ala_cb_id = setup
            .system
            .residue(setup.ala_id)
            .unwrap()
            .get_first_atom_id_by_name("CB")
            .unwrap();
        let leu_cb_id = setup
            .system
            .residue(setup.leu_id)
            .unwrap()
            .get_first_atom_id_by_name("CB")
            .unwrap();
        let leu_cg_id = setup
            .system
            .residue(setup.leu_id)
            .unwrap()
            .get_first_atom_id_by_name("CG")
            .unwrap();

        setup.system.atom_mut(ala_cb_id).unwrap().position = Point3::new(0.0, 0.0, 0.0);
        setup.system.atom_mut(leu_cb_id).unwrap().position = Point3::new(4.0, 0.0, 0.0);
        setup.system.atom_mut(leu_cg_id).unwrap().position = Point3::new(100.0, 100.0, 100.0);

        let expected_vdw = crate::core::forcefield::potentials::lennard_jones_12_6(4.0, 4.0, 0.12);
        assert!(
            (expected_vdw - (-0.12)).abs() < TOLERANCE,
            "Manual calculation sanity check failed"
        );

        let active_residues: HashSet<_> = [setup.ala_id, setup.leu_id].iter().cloned().collect();
        let energy = run(&setup.system, &setup.forcefield, &active_residues).unwrap();

        assert!(
            (energy.vdw - expected_vdw).abs() < TOLERANCE,
            "VDW energy mismatch: expected {}, got {}",
            expected_vdw,
            energy.vdw
        );
    }

    #[test]
    fn run_correctly_sums_multiple_pair_interactions() {
        let setup = setup();
        let active_residues: HashSet<_> = [setup.ala_id, setup.gly_id, setup.leu_id]
            .iter()
            .cloned()
            .collect();

        let scorer = Scorer::new(&setup.system, &setup.forcefield);

        let ala_sc_atoms = setup
            .system
            .residue(setup.ala_id)
            .unwrap()
            .atoms()
            .iter()
            .filter(|id| setup.system.atom(**id).unwrap().role == AtomRole::Sidechain)
            .copied()
            .collect::<Vec<_>>();
        let leu_sc_atoms = setup
            .system
            .residue(setup.leu_id)
            .unwrap()
            .atoms()
            .iter()
            .filter(|id| setup.system.atom(**id).unwrap().role == AtomRole::Sidechain)
            .copied()
            .collect::<Vec<_>>();
        let expected_energy = scorer
            .score_interaction(&ala_sc_atoms, &leu_sc_atoms)
            .unwrap();

        let total_energy = run(&setup.system, &setup.forcefield, &active_residues).unwrap();

        assert!((total_energy.total() - expected_energy.total()).abs() < TOLERANCE);
    }

    #[test]
    fn collect_active_sidechain_atoms_works_correctly() {
        let setup = setup();
        let active_residues: HashSet<_> = [setup.ala_id, setup.gly_id, setup.leu_id]
            .iter()
            .cloned()
            .collect();

        let map = collect_active_sidechain_atoms(&setup.system, &active_residues);

        assert_eq!(map.len(), 3);
        assert_eq!(
            map.get(&setup.ala_id).unwrap().len(),
            1,
            "ALA should have 1 sidechain atom"
        );
        assert!(
            map.get(&setup.gly_id).unwrap().is_empty(),
            "GLY should have 0 sidechain atoms"
        );
        assert_eq!(
            map.get(&setup.leu_id).unwrap().len(),
            2,
            "LEU should have 2 sidechain atoms"
        );

        let ala_cb_id = setup
            .system
            .residue(setup.ala_id)
            .unwrap()
            .get_first_atom_id_by_name("CB")
            .unwrap();
        assert_eq!(map.get(&setup.ala_id).unwrap()[0], ala_cb_id);
    }
}
