use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::atom::AtomRole;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::core::models::residue::ResidueType;
use crate::core::models::system::MolecularSystem;
use crate::engine::cache::ELCache;
use crate::engine::config::DesignSpecExt;
use crate::engine::context::{OptimizationContext, ProvidesResidueSelections};
use crate::engine::error::EngineError;
use crate::engine::placement::place_rotamer_on_system;
use crate::engine::progress::Progress;
use crate::engine::utils::query::precompute_environment_atoms;
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

#[instrument(skip_all, name = "el_energy_cache_generation_task")]
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
        context.reporter.report(Progress::PhaseFinish);
        return Ok(ELCache::new());
    }

    let active_residue_ids = context.resolve_all_active_residues()?;
    let environment_atom_ids = precompute_environment_atoms(context.system, &active_residue_ids);

    context.reporter.report(Progress::TaskStart {
        total: work_list.len() as u64,
    });

    #[cfg(not(feature = "parallel"))]
    let iterator = work_list.iter();

    #[cfg(feature = "parallel")]
    let iterator = work_list.par_iter();

    let results: Vec<WorkResult> = iterator
        .map(|unit| {
            #[cfg(feature = "parallel")]
            {
                let mut local_system = context.system.clone();
                compute_energies_for_unit(unit, &environment_atom_ids, context, &mut local_system)
            }

            #[cfg(not(feature = "parallel"))]
            {
                let mut temp_system = context.system.clone();
                compute_energies_for_unit(unit, &environment_atom_ids, context, &mut temp_system)
            }
        })
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
        "EL pre-computation finished. Cached {} residue-type combinations.",
        cache.len()
    );
    context.reporter.report(Progress::PhaseFinish);
    Ok(cache)
}

#[instrument(skip_all, name = "current_el_energy_task")]
pub fn calculate_current<C>(context: &OptimizationContext<C>) -> Result<EnergyTerm, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    info!("Calculating current total EL energy for all active sidechains.");

    let active_residues = context.resolve_all_active_residues()?;
    if active_residues.is_empty() {
        info!("No active sidechains found. Current EL energy is zero.");
        return Ok(EnergyTerm::default());
    }

    let scorer = Scorer::new(context.system, context.forcefield);
    let environment_atom_ids = precompute_environment_atoms(context.system, &active_residues);

    let mut total_el_energy = EnergyTerm::default();

    for residue_id in active_residues {
        let sidechain_atoms: Vec<AtomId> = context
            .system
            .residue(residue_id)
            .unwrap()
            .atoms()
            .iter()
            .filter_map(|&atom_id| {
                context.system.atom(atom_id).and_then(|atom| {
                    if atom.role == AtomRole::Sidechain {
                        Some(atom_id)
                    } else {
                        None
                    }
                })
            })
            .collect();

        if sidechain_atoms.is_empty() {
            continue;
        }

        // 1. Interaction with fixed environment for THIS sidechain
        let interaction_energy =
            scorer.score_interaction(&sidechain_atoms, &environment_atom_ids)?;

        // 2. Internal non-bonded energy for THIS sidechain
        let internal_energy = scorer.score_group_internal(&sidechain_atoms)?;

        // 3. Add this residue's complete EL energy to the total
        total_el_energy += interaction_energy + internal_energy;
    }

    info!(
        energy = total_el_energy.total(),
        "Sum of individual EL energies for the current conformation has been calculated."
    );
    Ok(total_el_energy)
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
            if let Some(allowed_types) =
                design_spec.get_by_specifier(chain.id, residue.residue_number)
            {
                is_design_site = true;
                for &residue_type in allowed_types {
                    if context
                        .rotamer_library
                        .get_rotamers_for(residue_type)
                        .is_some()
                    {
                        work_list.push(WorkUnit {
                            residue_id,
                            residue_type,
                        });
                    }
                }
            }
        }

        if !is_design_site {
            if let Some(native_type) = residue.residue_type {
                if context
                    .rotamer_library
                    .get_rotamers_for(native_type)
                    .is_some()
                {
                    work_list.push(WorkUnit {
                        residue_id,
                        residue_type: native_type,
                    });
                }
            }
        }
    }
    Ok(work_list)
}

#[instrument(skip_all, fields(residue_id = ?unit.residue_id, residue_type = %unit.residue_type))]
fn compute_energies_for_unit<C>(
    unit: &WorkUnit,
    environment_atom_ids: &[AtomId],
    context: &OptimizationContext<C>,
    system: &mut MolecularSystem,
) -> WorkResult
where
    C: ProvidesResidueSelections + Sync,
{
    let rotamers = context
        .rotamer_library
        .get_rotamers_for(unit.residue_type)
        .ok_or_else(|| EngineError::RotamerLibrary {
            residue_type: unit.residue_type.to_string(),
            message: "No rotamers found for this residue type during EL calculation.".to_string(),
        })?;

    let residue_name = unit.residue_type.to_three_letter();
    let topology = context.topology_registry.get(residue_name).ok_or_else(|| {
        EngineError::TopologyNotFound {
            residue_name: residue_name.to_string(),
        }
    })?;

    let mut energy_map = HashMap::with_capacity(rotamers.len());

    for (rotamer_idx, rotamer) in rotamers.iter().enumerate() {
        place_rotamer_on_system(system, unit.residue_id, rotamer, topology)?;

        let query_atoms: Vec<AtomId> = system
            .residue(unit.residue_id)
            .unwrap()
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

        if query_atoms.is_empty() {
            energy_map.insert(rotamer_idx, EnergyTerm::default());
            continue;
        }

        let scorer = Scorer::new(system, context.forcefield);

        let interaction_energy = scorer.score_interaction(&query_atoms, environment_atom_ids)?;
        let internal_energy = scorer.score_group_internal(&query_atoms)?;

        let total_el_energy = interaction_energy + internal_energy;

        energy_map.insert(rotamer_idx, total_el_energy);
    }

    context
        .reporter
        .report(Progress::TaskIncrement { amount: 1 });

    Ok(((unit.residue_id, unit.residue_type), energy_map))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{
        forcefield::params::{EnergyWeights, Forcefield},
        models::{atom::Atom, chain::ChainType, residue::ResidueType, system::MolecularSystem},
        rotamers::{library::RotamerLibrary, rotamer::Rotamer},
        topology::registry::TopologyRegistry,
    };
    use crate::engine::{
        config::{
            ConvergenceConfig, PlacementConfig, PlacementConfigBuilder, ResidueSelection,
            ResidueSpecifier,
        },
        context::OptimizationContext,
        progress::{Progress, ProgressReporter},
    };
    use std::collections::{HashMap, HashSet};
    use std::sync::{Arc, Mutex};
    use tempfile::tempdir;

    struct TestSetup {
        system: MolecularSystem,
        forcefield: Forcefield,
        rotamer_library: RotamerLibrary,
        topology_registry: TopologyRegistry,
        progress_events: Arc<Mutex<Vec<Progress>>>,
        _temp_dir: tempfile::TempDir,
    }

    impl TestSetup {
        fn new() -> Self {
            let temp_dir = tempdir().unwrap();
            let mut system = MolecularSystem::new();
            let forcefield = Self::create_forcefield(&temp_dir);
            let topology_registry = Self::create_topology_registry(&temp_dir);

            let chain_a = system.add_chain('A', ChainType::Protein);

            let ala_id = system
                .add_residue(chain_a, 1, "ALA", Some(ResidueType::Alanine))
                .unwrap();
            Self::add_atom(&mut system, ala_id, "N", [0.0, 1.4, 0.0]);
            Self::add_atom(&mut system, ala_id, "CA", [0.0, 0.0, 0.0]);
            Self::add_atom(&mut system, ala_id, "C", [1.4, 0.0, 0.0]);
            Self::add_atom(&mut system, ala_id, "CB", [0.0, -1.5, 0.0]);

            let ser_id = system
                .add_residue(chain_a, 2, "SER", Some(ResidueType::Serine))
                .unwrap();
            Self::add_atom(&mut system, ser_id, "N", [10.0, 0.0, 0.0]);
            Self::add_atom(&mut system, ser_id, "CA", [11.0, 0.0, 0.0]);
            Self::add_atom(&mut system, ser_id, "C", [12.0, 0.0, 0.0]);
            Self::add_atom(&mut system, ser_id, "CB", [11.0, 1.5, 0.0]);

            let leu_id = system
                .add_residue(chain_a, 3, "LEU", Some(ResidueType::Leucine))
                .unwrap();
            Self::add_atom(&mut system, leu_id, "N", [19.0, 0.0, 0.0]);
            Self::add_atom(&mut system, leu_id, "CA", [20.0, 0.0, 0.0]);
            Self::add_atom(&mut system, leu_id, "C", [21.0, 0.0, 0.0]);

            let chain_b = system.add_chain('B', ChainType::Water);
            let hoh_id = system.add_residue(chain_b, 1, "HOH", None).unwrap();
            Self::add_atom(&mut system, hoh_id, "O", [5.0, 5.0, 5.0]);

            let rotamer_library = Self::create_rotamer_library(&topology_registry, &forcefield);

            let parameterizer = crate::core::forcefield::parameterization::Parameterizer::new(
                &forcefield,
                &topology_registry,
                0.0,
            );
            parameterizer.parameterize_system(&mut system).unwrap();

            let progress_events = Arc::new(Mutex::new(Vec::new()));

            Self {
                system,
                forcefield,
                rotamer_library,
                topology_registry,
                progress_events,
                _temp_dir: temp_dir,
            }
        }

        fn create_forcefield(dir: &tempfile::TempDir) -> Forcefield {
            let ff_path = dir.path().join("ff.toml");
            std::fs::write(
                &ff_path,
                r#"
                [globals]
                dielectric_constant = 4.0
                potential_function = "lennard-jones-12-6"
                [vdw]
                N_BB = { radius = 2.8, well_depth = 0.15 }
                C_BB = { radius = 3.0, well_depth = 0.1 }
                C_SC = { radius = 4.0, well_depth = 0.2 }
                O_W  = { radius = 3.2, well_depth = 0.3 }
                [hbond]
            "#,
            )
            .unwrap();
            let delta_path = dir.path().join("delta.csv");
            std::fs::write(&delta_path, "residue_type,atom_name,mu,sigma\n").unwrap();
            Forcefield::load(&ff_path, &delta_path, &EnergyWeights::default()).unwrap()
        }

        fn create_topology_registry(dir: &tempfile::TempDir) -> TopologyRegistry {
            let file_path = dir.path().join("registry.toml");
            std::fs::write(
                &file_path,
                r#"
                [ALA]
                anchor_atoms = ["N", "CA", "C"]
                sidechain_atoms = ["CB"]

                [SER]
                anchor_atoms = ["N", "CA", "C"]
                sidechain_atoms = ["CB"]

                [LEU]
                anchor_atoms = ["N", "CA", "C"]
                sidechain_atoms = ["CB", "CG"]
            "#,
            )
            .unwrap();
            TopologyRegistry::load(&file_path).unwrap()
        }

        fn create_rotamer_library(
            topology_registry: &TopologyRegistry,
            forcefield: &Forcefield,
        ) -> RotamerLibrary {
            let parameterizer = crate::core::forcefield::parameterization::Parameterizer::new(
                forcefield,
                topology_registry,
                0.0,
            );
            let mut rotamers = HashMap::new();
            let ala_rotamers = vec![
                Self::build_parameterized_rotamer(
                    "ALA",
                    &parameterizer,
                    topology_registry,
                    &[("CB", [0.0, -0.7, 1.2])],
                ),
                Self::build_parameterized_rotamer(
                    "ALA",
                    &parameterizer,
                    topology_registry,
                    &[("CB", [0.0, -1.5, 0.0])],
                ),
            ];
            rotamers.insert(ResidueType::Alanine, ala_rotamers);

            let ser_rotamers = vec![Self::build_parameterized_rotamer(
                "SER",
                &parameterizer,
                topology_registry,
                &[("CB", [11.0, 1.5, 0.0])],
            )];
            rotamers.insert(ResidueType::Serine, ser_rotamers);

            RotamerLibrary { rotamers }
        }

        fn build_parameterized_rotamer(
            res_name: &str,
            parameterizer: &crate::core::forcefield::parameterization::Parameterizer,
            topology_registry: &TopologyRegistry,
            sidechain_atoms_data: &[(&str, [f64; 3])],
        ) -> Rotamer {
            let mut atoms = vec![
                Atom::new("N", ResidueId::default(), [0.0, 1.4, 0.0].into()),
                Atom::new("CA", ResidueId::default(), [0.0, 0.0, 0.0].into()),
                Atom::new("C", ResidueId::default(), [1.4, 0.0, 0.0].into()),
            ];
            for (name, pos) in sidechain_atoms_data {
                atoms.push(Atom::new(name, ResidueId::default(), (*pos).into()));
            }

            atoms.iter_mut().for_each(|a| {
                a.force_field_type = match a.name.as_str() {
                    "N" => "N_BB".to_string(),
                    "CA" | "C" => "C_BB".to_string(),
                    _ => "C_SC".to_string(),
                }
            });

            let bonds = match res_name {
                "ALA" | "SER" => vec![(0, 1), (1, 2), (1, 3)],
                _ => vec![],
            };

            let topology = topology_registry.get(res_name).unwrap();
            let mut rotamer = Rotamer { atoms, bonds };
            parameterizer
                .parameterize_rotamer(&mut rotamer, res_name, topology)
                .unwrap();
            rotamer
        }

        fn add_atom(system: &mut MolecularSystem, res_id: ResidueId, name: &str, pos: [f64; 3]) {
            let mut atom = Atom::new(name, res_id, pos.into());
            atom.force_field_type = match name {
                "N" => "N_BB".to_string(),
                "CA" | "C" => "C_BB".to_string(),
                "CB" => "C_SC".to_string(),
                "O" => "O_W".to_string(),
                _ => "UNKNOWN".to_string(),
            };
            system.add_atom_to_residue(res_id, atom).unwrap();
        }

        fn create_placement_config(&self, selection: ResidueSelection) -> PlacementConfig {
            PlacementConfigBuilder::new()
                .forcefield_path("dummy.ff")
                .delta_params_path("dummy.delta")
                .s_factor(0.0)
                .rotamer_library_path("dummy.rotlib")
                .topology_registry_path("dummy.topo")
                .max_iterations(10)
                .num_solutions(1)
                .convergence_config(ConvergenceConfig {
                    energy_threshold: 0.01,
                    patience_iterations: 3,
                })
                .final_refinement_iterations(2)
                .include_input_conformation(false)
                .residues_to_optimize(selection)
                .build()
                .unwrap()
        }

        fn get_residue_id(&self, chain: char, res_num: isize) -> ResidueId {
            let chain_id = self.system.find_chain_by_id(chain).unwrap();
            self.system.find_residue_by_id(chain_id, res_num).unwrap()
        }

        fn reporter(&self) -> ProgressReporter {
            let events = self.progress_events.clone();
            ProgressReporter::with_callback(Box::new(move |p: Progress| {
                events.lock().unwrap().push(p);
            }))
        }
    }

    #[test]
    fn run_calculates_correct_el_energy_value() {
        let setup = TestSetup::new();
        let ala_id = setup.get_residue_id('A', 1);
        let config = setup.create_placement_config(ResidueSelection::List {
            include: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 1,
            }],
            exclude: vec![],
        });
        let reporter = setup.reporter();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let cache = run(&context).unwrap();

        let energies = cache
            .get_energies_for(ala_id, ResidueType::Alanine)
            .unwrap();
        let calculated_energy = energies.get(&1).unwrap();

        let rotamer = &setup
            .rotamer_library
            .get_rotamers_for(ResidueType::Alanine)
            .unwrap()[1];

        let mut temp_system = setup.system.clone();
        let topology = setup.topology_registry.get("ALA").unwrap();
        place_rotamer_on_system(&mut temp_system, ala_id, rotamer, topology).unwrap();

        let query_atom_ids: Vec<AtomId> = temp_system
            .residue(ala_id)
            .unwrap()
            .atoms()
            .iter()
            .filter_map(|&atom_id| {
                temp_system.atom(atom_id).and_then(|atom| {
                    if atom.role == AtomRole::Sidechain {
                        Some(atom_id)
                    } else {
                        None
                    }
                })
            })
            .collect();

        let env_atom_ids = precompute_environment_atoms(&temp_system, &HashSet::from([ala_id]));

        let scorer = Scorer::new(&temp_system, &setup.forcefield);

        let interaction_energy = scorer
            .score_interaction(&query_atom_ids, &env_atom_ids)
            .unwrap();
        let internal_energy = scorer.score_group_internal(&query_atom_ids).unwrap();
        let expected_energy = interaction_energy + internal_energy;

        assert!((calculated_energy.vdw - expected_energy.vdw).abs() < 1e-9);
        assert!((calculated_energy.coulomb - expected_energy.coulomb).abs() < 1e-9);
        assert!((calculated_energy.hbond - expected_energy.hbond).abs() < 1e-9);
        assert!((calculated_energy.total() - expected_energy.total()).abs() < 1e-9);
    }

    #[test]
    fn calculate_current_returns_zero_for_no_active_residues() {
        let setup = TestSetup::new();
        let config = setup.create_placement_config(ResidueSelection::List {
            include: vec![],
            exclude: vec![],
        });
        let reporter = setup.reporter();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let energy = calculate_current(&context).unwrap();

        assert_eq!(energy.total(), 0.0);
    }

    #[test]
    fn calculate_current_returns_zero_for_active_residue_with_no_sidechain() {
        let setup = TestSetup::new();
        let config = setup.create_placement_config(ResidueSelection::List {
            include: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 2,
            }],
            exclude: vec![],
        });
        let reporter = setup.reporter();

        let mut modified_system = setup.system.clone();
        let ser_id = setup.get_residue_id('A', 2);
        let ser_cb_id = modified_system
            .residue(ser_id)
            .unwrap()
            .get_first_atom_id_by_name("CB")
            .unwrap();
        modified_system.remove_atom(ser_cb_id);

        let context = OptimizationContext::new(
            &modified_system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let energy = calculate_current(&context).unwrap();
        assert_eq!(
            energy.total(),
            0.0,
            "Energy should be zero as the only active residue has no sidechain"
        );
    }

    #[test]
    fn calculate_current_calculates_energy_for_a_single_active_residue() {
        let setup = TestSetup::new();
        let config = setup.create_placement_config(ResidueSelection::List {
            include: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 1,
            }],
            exclude: vec![],
        });
        let reporter = setup.reporter();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let scorer = Scorer::new(&setup.system, &setup.forcefield);
        let ala_id = setup.get_residue_id('A', 1);
        let ala_sc_atoms = setup
            .system
            .residue(ala_id)
            .unwrap()
            .atoms()
            .iter()
            .filter(|id| setup.system.atom(**id).unwrap().role == AtomRole::Sidechain)
            .copied()
            .collect::<Vec<_>>();
        let all_other_atoms = setup
            .system
            .atoms_iter()
            .map(|(id, _)| id)
            .filter(|id| !ala_sc_atoms.contains(id))
            .collect::<Vec<_>>();
        let interaction = scorer
            .score_interaction(&ala_sc_atoms, &all_other_atoms)
            .unwrap();
        let internal = scorer.score_group_internal(&ala_sc_atoms).unwrap();
        let expected_energy = interaction + internal;

        let calculated_energy = calculate_current(&context).unwrap();

        assert!((calculated_energy.total() - expected_energy.total()).abs() < 1e-9);
        assert!(
            expected_energy.total().abs() > 1e-6,
            "Expected a non-zero energy"
        );
    }

    #[test]
    fn calculate_current_calculates_energy_for_multiple_active_residues() {
        let setup = TestSetup::new();
        let config = setup.create_placement_config(ResidueSelection::List {
            include: vec![
                ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: 1,
                },
                ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: 2,
                },
            ],
            exclude: vec![],
        });
        let reporter = setup.reporter();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let scorer = Scorer::new(&setup.system, &setup.forcefield);
        let active_residue_ids = context.resolve_all_active_residues().unwrap();
        let env_atoms = precompute_environment_atoms(&setup.system, &active_residue_ids);
        let mut expected_energy = EnergyTerm::default();
        let active_residue_ids = context.resolve_all_active_residues().unwrap();
        for res_id in active_residue_ids {
            let sc_atoms: Vec<AtomId> = setup
                .system
                .residue(res_id)
                .unwrap()
                .atoms()
                .iter()
                .filter(|id| setup.system.atom(**id).unwrap().role == AtomRole::Sidechain)
                .copied()
                .collect();
            if sc_atoms.is_empty() {
                continue;
            }
            let interaction = scorer.score_interaction(&sc_atoms, &env_atoms).unwrap();
            let internal = scorer.score_group_internal(&sc_atoms).unwrap();
            expected_energy = expected_energy + interaction + internal;
        }

        let calculated_energy = calculate_current(&context).unwrap();

        assert!((calculated_energy.total() - expected_energy.total()).abs() < 1e-9);
        assert!(
            expected_energy.total().abs() > 1e-6,
            "Expected a non-zero energy for multiple sidechains"
        );
    }
}
