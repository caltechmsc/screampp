use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::atom::AtomRole;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::core::models::residue::ResidueType;
use crate::engine::cache::ELCache;
use crate::engine::config::DesignSpecExt;
use crate::engine::context::{OptimizationContext, ProvidesResidueSelections};
use crate::engine::error::EngineError;
use crate::engine::placement::place_rotamer_on_system;
use crate::engine::progress::Progress;
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

    let environment_atom_ids = precompute_environment_atoms(context)?;

    context.reporter.report(Progress::TaskStart {
        total: work_list.len() as u64,
    });

    #[cfg(not(feature = "parallel"))]
    let iterator = work_list.iter();

    #[cfg(feature = "parallel")]
    let iterator = work_list.par_iter();

    let results: Vec<WorkResult> = iterator
        .map(|unit| compute_energies_for_unit(unit, &environment_atom_ids, context))
        .collect();

    context.reporter.report(Progress::TaskFinish);

    let mut cache = ELCache::new();
    for result in results {
        let ((residue_id, residue_type), energy_map) = result?;
        for (rotamer_idx, energy_term) in energy_map {
            // TODO: Add pre-calculated internal energy from rotamer library to `energy_term`.
            cache.insert(residue_id, residue_type, rotamer_idx, energy_term);
        }
    }

    info!(
        cached_combinations = cache.len(),
        "EL pre-computation finished."
    );
    context.reporter.report(Progress::PhaseFinish);
    Ok(cache)
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

fn precompute_environment_atoms<C>(
    context: &OptimizationContext<C>,
) -> Result<Vec<AtomId>, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    let active_residue_ids = context.resolve_all_active_residues()?;

    let environment_atom_ids = context
        .system
        .atoms_iter()
        .filter_map(|(atom_id, atom)| match atom.role {
            AtomRole::Ligand | AtomRole::Water | AtomRole::Other => Some(atom_id),
            AtomRole::Backbone => Some(atom_id),
            AtomRole::Sidechain => {
                if !active_residue_ids.contains(&atom.residue_id) {
                    Some(atom_id)
                } else {
                    None
                }
            }
        })
        .collect();

    Ok(environment_atom_ids)
}

fn collect_active_sidechain_atoms<C>(
    context: &OptimizationContext<C>,
) -> Result<Vec<AtomId>, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    let active_residue_ids = context.resolve_all_active_residues()?;
    Ok(context
        .system
        .atoms_iter()
        .filter_map(|(atom_id, atom)| {
            if atom.role == AtomRole::Sidechain && active_residue_ids.contains(&atom.residue_id) {
                Some(atom_id)
            } else {
                None
            }
        })
        .collect())
}

#[instrument(skip_all, fields(residue_id = ?unit.residue_id, residue_type = %unit.residue_type))]
fn compute_energies_for_unit<C>(
    unit: &WorkUnit,
    environment_atom_ids: &[AtomId],
    context: &OptimizationContext<C>,
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
        let mut temp_system = context.system.clone();

        place_rotamer_on_system(&mut temp_system, unit.residue_id, rotamer, topology)?;

        let query_atoms: Vec<AtomId> = temp_system
            .residue(unit.residue_id)
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

        if query_atoms.is_empty() {
            energy_map.insert(rotamer_idx, EnergyTerm::default());
            continue;
        }

        let scorer = Scorer::new(&temp_system, context.forcefield);
        let interaction_energy = scorer.score_interaction(&query_atoms, environment_atom_ids)?;

        energy_map.insert(rotamer_idx, interaction_energy);
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
        forcefield::params::{Forcefield, GlobalParams, NonBondedParams, VdwParam},
        models::{atom::Atom, chain::ChainType, residue::ResidueType, system::MolecularSystem},
        rotamers::{library::RotamerLibrary, rotamer::Rotamer},
        topology::registry::TopologyRegistry,
    };
    use crate::engine::{
        config::{
            ConvergenceConfig, DesignConfig, DesignConfigBuilder, DesignSpec, PlacementConfig,
            PlacementConfigBuilder, ResidueSelection, ResidueSpecifier,
        },
        context::OptimizationContext,
        progress::{Progress, ProgressReporter},
    };
    use nalgebra::Point3;
    use std::collections::{HashMap, HashSet};
    use std::sync::{Arc, Mutex};
    use tempfile::tempdir;

    struct TestSetup {
        system: MolecularSystem,
        forcefield: Forcefield,
        rotamer_library: RotamerLibrary,
        topology_registry: TopologyRegistry,
        progress_events: Arc<Mutex<Vec<Progress>>>,
    }

    impl TestSetup {
        fn new() -> Self {
            let mut system = MolecularSystem::new();
            let forcefield = Self::create_forcefield();
            let topology_registry = Self::create_topology_registry();

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
            }
        }

        fn create_forcefield() -> Forcefield {
            let mut vdw = HashMap::new();
            vdw.insert(
                "N_BB".to_string(),
                VdwParam::LennardJones {
                    radius: 2.8,
                    well_depth: 0.15,
                },
            );
            vdw.insert(
                "C_BB".to_string(),
                VdwParam::LennardJones {
                    radius: 3.0,
                    well_depth: 0.1,
                },
            );
            vdw.insert(
                "C_SC".to_string(),
                VdwParam::LennardJones {
                    radius: 4.0,
                    well_depth: 0.2,
                },
            );
            vdw.insert(
                "O_W".to_string(),
                VdwParam::LennardJones {
                    radius: 3.2,
                    well_depth: 0.3,
                },
            );
            Forcefield {
                non_bonded: NonBondedParams {
                    globals: GlobalParams {
                        dielectric_constant: 1.0,
                        ..Default::default()
                    },
                    vdw,
                    hbond: HashMap::new(),
                },
                deltas: HashMap::new(),
            }
        }

        fn create_topology_registry() -> TopologyRegistry {
            let dir = tempdir().unwrap();
            let file_path = dir.path().join("registry.toml");
            let registry_content = r#"
                [ALA]
                anchor_atoms = ["N", "CA", "C"]
                sidechain_atoms = ["CB"]

                [SER]
                anchor_atoms = ["N", "CA", "C"]
                sidechain_atoms = ["CB"]

                [LEU]
                anchor_atoms = ["N", "CA", "C"]
                sidechain_atoms = ["CB", "CG"]
            "#;
            std::fs::write(&file_path, registry_content).unwrap();
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
                "ALA" => vec![(0, 1), (1, 2), (1, 3)],
                "SER" => vec![(0, 1), (1, 2), (1, 3)],
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

        fn create_placement_config(selection: ResidueSelection) -> PlacementConfig {
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

        fn create_design_config(
            design_spec: DesignSpec,
            repack_selection: ResidueSelection,
        ) -> DesignConfig {
            DesignConfigBuilder::new()
                .forcefield_path("dummy.ff")
                .delta_params_path("dummy.delta")
                .s_factor(0.0)
                .rotamer_library_path("dummy.rotlib")
                .topology_registry_path("dummy.topo")
                .max_iterations(10)
                .num_solutions(1)
                .include_input_conformation(false)
                .convergence_config(ConvergenceConfig {
                    energy_threshold: 0.01,
                    patience_iterations: 3,
                })
                .final_refinement_iterations(2)
                .design_spec(design_spec)
                .neighbors_to_repack(repack_selection)
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
    fn precompute_environment_atoms_is_correct() {
        let setup = TestSetup::new();
        let config = TestSetup::create_placement_config(ResidueSelection::All);
        let reporter = setup.reporter();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let env_atoms = precompute_environment_atoms(&context).unwrap();

        let env_atom_info: HashSet<(isize, String)> = env_atoms
            .iter()
            .map(|id| {
                let atom = context.system.atom(*id).unwrap();
                let res = context.system.residue(atom.residue_id).unwrap();
                (res.residue_number, atom.name.clone())
            })
            .collect();

        assert!(!env_atom_info.contains(&(1, "CB".to_string())));
        assert!(!env_atom_info.contains(&(2, "CB".to_string())));

        assert!(env_atom_info.contains(&(1, "N".to_string())));
        assert!(env_atom_info.contains(&(1, "CA".to_string())));
        assert!(env_atom_info.contains(&(1, "C".to_string())));
        assert!(env_atom_info.contains(&(2, "N".to_string())));
        assert!(env_atom_info.contains(&(2, "CA".to_string())));
        assert!(env_atom_info.contains(&(1, "O".to_string())));
    }

    #[test]
    fn run_calculates_correct_el_energy_value() {
        let setup = TestSetup::new();
        let ala_id = setup.get_residue_id('A', 1);
        let config = TestSetup::create_placement_config(ResidueSelection::List {
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

        let rotamer = &setup
            .rotamer_library
            .get_rotamers_for(ResidueType::Alanine)
            .unwrap()[1];
        let rot_cb_pos = rotamer
            .atoms
            .iter()
            .find(|a| a.name == "CB")
            .unwrap()
            .position;

        let ser_res_id = setup.get_residue_id('A', 2);
        let ser_n = setup
            .system
            .atom(
                setup
                    .system
                    .residue(ser_res_id)
                    .unwrap()
                    .get_first_atom_id_by_name("N")
                    .unwrap(),
            )
            .unwrap();
        let ser_ca = setup
            .system
            .atom(
                setup
                    .system
                    .residue(ser_res_id)
                    .unwrap()
                    .get_first_atom_id_by_name("CA")
                    .unwrap(),
            )
            .unwrap();
        let ser_c = setup
            .system
            .atom(
                setup
                    .system
                    .residue(ser_res_id)
                    .unwrap()
                    .get_first_atom_id_by_name("C")
                    .unwrap(),
            )
            .unwrap();
        let ser_cb = setup
            .system
            .atom(
                setup
                    .system
                    .residue(ser_res_id)
                    .unwrap()
                    .get_first_atom_id_by_name("CB")
                    .unwrap(),
            )
            .unwrap();

        let leu_res_id = setup.get_residue_id('A', 3);
        let leu_n = setup
            .system
            .atom(
                setup
                    .system
                    .residue(leu_res_id)
                    .unwrap()
                    .get_first_atom_id_by_name("N")
                    .unwrap(),
            )
            .unwrap();
        let leu_ca = setup
            .system
            .atom(
                setup
                    .system
                    .residue(leu_res_id)
                    .unwrap()
                    .get_first_atom_id_by_name("CA")
                    .unwrap(),
            )
            .unwrap();
        let leu_c = setup
            .system
            .atom(
                setup
                    .system
                    .residue(leu_res_id)
                    .unwrap()
                    .get_first_atom_id_by_name("C")
                    .unwrap(),
            )
            .unwrap();

        let hoh_o = setup
            .system
            .atom(
                setup
                    .system
                    .residue(setup.get_residue_id('B', 1))
                    .unwrap()
                    .get_first_atom_id_by_name("O")
                    .unwrap(),
            )
            .unwrap();

        let params = HashMap::from([
            ("N_BB", (2.8, 0.15)),
            ("C_BB", (3.0, 0.1)),
            ("C_SC", (4.0, 0.2)),
            ("O_W", (3.2, 0.3)),
        ]);
        let (cb_r, cb_d) = params["C_SC"];

        let env_atoms = [ser_n, ser_ca, ser_c, ser_cb, leu_n, leu_ca, leu_c, hoh_o];

        let mut expected_vdw_by_scorer = 0.0;

        fn calc_vdw(
            pos1: Point3<f64>,
            pos2: Point3<f64>,
            r1: f64,
            d1: f64,
            r2: f64,
            d2: f64,
        ) -> f64 {
            let dist = (pos1 - pos2).norm();
            let r_min = (r1 + r2) / 2.0;
            let well_depth = (d1 * d2).sqrt();
            let rho = r_min / dist;
            let rho6 = rho.powi(6);
            well_depth * (rho6 * rho6 - 2.0 * rho6)
        }

        for env_atom in &env_atoms {
            expected_vdw_by_scorer += calc_vdw(
                rot_cb_pos,
                env_atom.position,
                cb_r,
                cb_d,
                params[env_atom.force_field_type.as_str()].0,
                params[env_atom.force_field_type.as_str()].1,
            );
        }

        let calculated_energy = energies.get(&1).unwrap();

        assert!(
            (calculated_energy.vdw - expected_vdw_by_scorer).abs() < 1e-9,
            "Calculated VDW energy differs from manual calculation. Got {}, expected {}",
            calculated_energy.vdw,
            expected_vdw_by_scorer
        );
        assert!(
            (calculated_energy.coulomb).abs() < 1e-9,
            "Coulomb energy should be near zero for this test setup"
        );
    }

    #[test]
    fn build_work_list_handles_design_site_correctly() {
        let setup = TestSetup::new();
        let reporter = setup.reporter();

        let mut design_spec = DesignSpec::new();
        design_spec.insert(
            ResidueSpecifier {
                chain_id: 'A',
                residue_number: 1,
            },
            vec![ResidueType::Serine],
        );
        let config = TestSetup::create_design_config(
            design_spec,
            ResidueSelection::List {
                include: vec![],
                exclude: vec![],
            },
        );

        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let work_list = build_work_list(&context).unwrap();

        assert_eq!(
            work_list.len(),
            1,
            "Should only have one work unit for the design site"
        );
        let work_unit = &work_list[0];
        assert_eq!(work_unit.residue_id, setup.get_residue_id('A', 1));
        assert_eq!(
            work_unit.residue_type,
            ResidueType::Serine,
            "Work unit should be for SER, not the native ALA"
        );
    }

    #[test]
    fn run_with_no_active_residues_returns_empty_cache() {
        let setup = TestSetup::new();
        let config = TestSetup::create_placement_config(ResidueSelection::List {
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

        let cache = run(&context).unwrap();
        assert!(cache.is_empty());
        assert_eq!(setup.progress_events.lock().unwrap().len(), 2);
    }

    #[test]
    fn run_handles_residue_with_no_rotamers_in_library() {
        let setup = TestSetup::new();
        let config = TestSetup::create_placement_config(ResidueSelection::List {
            include: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 3,
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

        assert!(
            cache.is_empty(),
            "Cache should be empty as LEU has no rotamers to process"
        );
    }
}
