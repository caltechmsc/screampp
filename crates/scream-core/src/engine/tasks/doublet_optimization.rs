use crate::core::forcefield::scoring::Scorer;
use crate::core::models::atom::AtomRole;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::core::models::system::MolecularSystem;
use crate::engine::cache::ELCache;
use crate::engine::context::{OptimizationContext, ProvidesResidueSelections};
use crate::engine::error::EngineError;
use crate::engine::transaction::SystemView;
use std::collections::{HashMap, HashSet};
use tracing::{debug, instrument};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Debug, Clone, Copy)]
pub struct DoubletResult {
    pub rotamer_idx_a: usize,
    pub rotamer_idx_b: usize,
    pub best_local_energy: f64,
}

#[instrument(skip_all, name = "doublet_optimization_task", fields(res_a = ?res_a_id, res_b = ?res_b_id))]
pub fn run<C>(
    res_a_id: ResidueId,
    res_b_id: ResidueId,
    system: &MolecularSystem,
    el_cache: &ELCache,
    context: &OptimizationContext<C>,
    active_residues: &HashSet<ResidueId>,
) -> Result<DoubletResult, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    // --- 1. Prepare data ---
    let residue_a = system.residue(res_a_id).ok_or_else(|| {
        EngineError::Internal(format!("Residue {:?} not found in system", res_a_id))
    })?;
    let residue_b = system.residue(res_b_id).ok_or_else(|| {
        EngineError::Internal(format!("Residue {:?} not found in system", res_b_id))
    })?;

    let res_type_a = residue_a.residue_type.ok_or_else(|| {
        EngineError::Internal(format!(
            "Residue {:?} is not a standard residue type",
            res_a_id
        ))
    })?;
    let res_type_b = residue_b.residue_type.ok_or_else(|| {
        EngineError::Internal(format!(
            "Residue {:?} is not a standard residue type",
            res_b_id
        ))
    })?;

    let rotamers_a = context
        .rotamer_library
        .get_rotamers_for(res_type_a)
        .ok_or_else(|| EngineError::RotamerLibrary {
            residue_type: res_type_a.to_string(),
            message: "No rotamers found for doublet optimization.".to_string(),
        })?;
    let rotamers_b = context
        .rotamer_library
        .get_rotamers_for(res_type_b)
        .ok_or_else(|| EngineError::RotamerLibrary {
            residue_type: res_type_b.to_string(),
            message: "No rotamers found for doublet optimization.".to_string(),
        })?;

    if rotamers_a.is_empty() || rotamers_b.is_empty() {
        return Err(EngineError::RotamerLibrary {
            residue_type: if rotamers_a.is_empty() {
                res_type_a.to_string()
            } else {
                res_type_b.to_string()
            },
            message: "Rotamer list is empty, cannot perform doublet optimization.".to_string(),
        });
    }

    let index_pairs: Vec<(usize, usize)> = (0..rotamers_a.len())
        .flat_map(|i| (0..rotamers_b.len()).map(move |j| (i, j)))
        .collect();

    debug!(
        "Optimizing doublet with {}x{} = {} total pairs.",
        rotamers_a.len(),
        rotamers_b.len(),
        index_pairs.len()
    );

    let other_active_residue_ids: Vec<ResidueId> = active_residues
        .iter()
        .filter(|&&id| id != res_a_id && id != res_b_id)
        .cloned()
        .collect();

    // --- Phase 2: High-Performance Parallel Evaluation with Thread-Local State ---
    let best_pair_result = {
        let iterator = {
            #[cfg(feature = "parallel")]
            {
                index_pairs.par_iter()
            }
            #[cfg(not(feature = "parallel"))]
            {
                index_pairs.iter()
            }
        };

        iterator.fold(
            || -> Result<(Option<DoubletResult>, MolecularSystem, HashMap<ResidueId, usize>, &OptimizationContext<C>), EngineError> {
                let thread_local_system = system.clone();
                let mut thread_local_rotamers = HashMap::new();
                thread_local_rotamers.insert(res_a_id, 0);
                thread_local_rotamers.insert(res_b_id, 0);
                let thread_local_context = context;
                Ok((None, thread_local_system, thread_local_rotamers, thread_local_context))
            },

            |mut acc, &(idx_a, idx_b)| {
                if acc.is_err() { return acc; }
                let (thread_best_result, thread_system, thread_rotamers, thread_context) = acc.as_mut().unwrap();

                let mut system_view = SystemView::new(thread_system, thread_context, thread_rotamers);

                match system_view.transaction_doublet(res_a_id, res_b_id, |view| {
                    view.apply_move(res_a_id, idx_a)?;
                    view.apply_move(res_b_id, idx_b)?;

                    let el_a = el_cache.get(res_a_id, res_type_a, idx_a).copied().unwrap_or_default();
                    let el_b = el_cache.get(res_b_id, res_type_b, idx_b).copied().unwrap_or_default();

                    let scorer = Scorer::new(view.system, thread_context.forcefield);

                    let get_sc_atoms = |sys: &MolecularSystem, res_id: ResidueId| -> Vec<AtomId> {
                        sys.residue(res_id)
                            .unwrap()
                            .atoms()
                            .iter()
                            .filter_map(|&id| {
                                sys.atom(id).and_then(|a| {
                                    if a.role == AtomRole::Sidechain {
                                        Some(id)
                                    } else {
                                        None
                                    }
                                })
                            })
                            .collect()
                    };

                    let atoms_a_sc = get_sc_atoms(view.system, res_a_id);
                    let atoms_b_sc = get_sc_atoms(view.system, res_b_id);
                    let other_active_sc_atoms: Vec<AtomId> = other_active_residue_ids
                        .iter()
                        .flat_map(|&id| get_sc_atoms(view.system, id))
                        .collect();

                    let inter_ab = scorer.score_interaction(&atoms_a_sc, &atoms_b_sc)?;
                    let inter_a_others = scorer.score_interaction(&atoms_a_sc, &other_active_sc_atoms)?;
                    let inter_b_others = scorer.score_interaction(&atoms_b_sc, &other_active_sc_atoms)?;

                    let total_interaction_energy = inter_ab + inter_a_others + inter_b_others;

                    Ok(el_a.total() + el_b.total() + total_interaction_energy.total())
                }) {
                    Ok(current_energy) => {
                        if thread_best_result.is_none() || current_energy < thread_best_result.as_ref().unwrap().best_local_energy {
                            *thread_best_result = Some(DoubletResult {
                                rotamer_idx_a: idx_a,
                                rotamer_idx_b: idx_b,
                                best_local_energy: current_energy,
                            });
                        }
                    }
                    Err(e) => return Err(e),
                }
                acc
            },
        )
        .filter_map(|res| res.ok())
        .map(|(res, _, _, _)| res)
        .reduce_with(|best, current| {
            match (best, current) {
                (Some(b), Some(c)) => Some(if b.best_local_energy < c.best_local_energy { b } else { c }),
                (Some(b), None) => Some(b),
                (None, Some(c)) => Some(c),
                (None, None) => None,
            }
        })
        .flatten()
    };

    match best_pair_result {
        Some(result) => {
            debug!(
                "Found best pair (A:{}, B:{}) with local energy {:.4} kcal/mol",
                result.rotamer_idx_a, result.rotamer_idx_b, result.best_local_energy
            );
            Ok(result)
        }
        None => Err(EngineError::PhaseFailed {
            phase: "Doublet Optimization",
            reason: "No valid rotamer pairs could be evaluated.".to_string(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{
        forcefield::{
            parameterization::Parameterizer,
            params::{Forcefield, GlobalParams, NonBondedParams, VdwParam},
        },
        models::{atom::Atom, chain::ChainType, residue::ResidueType, system::MolecularSystem},
        rotamers::{library::RotamerLibrary, rotamer::Rotamer},
        topology::registry::TopologyRegistry,
    };
    use crate::engine::{
        config::{ConvergenceConfig, PlacementConfigBuilder, ResidueSelection},
        context::OptimizationContext,
        progress::ProgressReporter,
    };
    use nalgebra::Point3;
    use std::{
        collections::{HashMap, HashSet},
        fs,
        path::Path,
    };
    use tempfile::TempDir;

    struct TestSetupBasic {
        system: MolecularSystem,
        res_a_id: ResidueId,
        res_b_id: ResidueId,
        el_cache: ELCache,
        forcefield: Forcefield,
        rotamer_library: RotamerLibrary,
        topology_registry: TopologyRegistry,
        _temp_dir: TempDir,
    }

    struct TestSetupWithEnv {
        system: MolecularSystem,
        res_a_id: ResidueId,
        res_b_id: ResidueId,
        res_c_id: ResidueId,
        el_cache: ELCache,
        forcefield: Forcefield,
        rotamer_library: RotamerLibrary,
        topology_registry: TopologyRegistry,
        _temp_dir: TempDir,
    }

    fn write_file(path: &Path, content: &str) {
        fs::write(path, content).expect("Failed to write temporary file for test");
    }

    fn setup_basic_test_data() -> TestSetupBasic {
        let temp_dir = tempfile::tempdir().expect("Failed to create temp dir");
        let dir_path = temp_dir.path();

        let topology_path = dir_path.join("topology.toml");
        write_file(
            &topology_path,
            r#"
[ALA]
anchor_atoms = ["N", "CA", "C"]
sidechain_atoms = ["CB"]

[LEU]
anchor_atoms = ["N", "CA", "C"]
sidechain_atoms = ["CB"]
"#,
        );

        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_a_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let res_b_id = system
            .add_residue(chain_id, 2, "LEU", Some(ResidueType::Leucine))
            .unwrap();

        let add_residue_atoms = |system: &mut MolecularSystem, res_id: ResidueId, offset: f64| {
            let backbone_atoms_data = vec![
                ("N", Point3::new(offset, 1.0, 0.0)),
                ("CA", Point3::new(offset, 0.0, 0.0)),
                ("C", Point3::new(offset + 1.0, 0.0, 0.0)),
            ];
            for (name, pos) in backbone_atoms_data {
                let mut atom = Atom::new(name, res_id, pos);
                atom.force_field_type = "BB".to_string();
                system.add_atom_to_residue(res_id, atom).unwrap();
            }
            let mut cb_atom = Atom::new("CB", res_id, Point3::new(offset, -0.5, 1.2));
            cb_atom.force_field_type = "C_SC".to_string();
            system.add_atom_to_residue(res_id, cb_atom).unwrap();
        };

        add_residue_atoms(&mut system, res_a_id, 0.0);
        add_residue_atoms(&mut system, res_b_id, 2.0);

        let mut vdw = HashMap::new();
        vdw.insert(
            "BB".to_string(),
            VdwParam::LennardJones {
                radius: 1.0,
                well_depth: 0.0,
            },
        );
        vdw.insert(
            "C_SC".to_string(),
            VdwParam::LennardJones {
                radius: 3.8,
                well_depth: 0.1,
            },
        );
        let forcefield = Forcefield {
            non_bonded: NonBondedParams {
                globals: GlobalParams {
                    dielectric_constant: 1.0,
                    potential_function: "lennard-jones-12-6".to_string(),
                },
                vdw,
                hbond: HashMap::new(),
                hbond_donors: HashSet::new(),
                hbond_acceptors: HashSet::new(),
            },
            deltas: HashMap::new(),
            weight_map: HashMap::new(),
        };

        let topology_registry = TopologyRegistry::load(&topology_path).unwrap();
        let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 0.0);

        let create_rotamer = |residue_id, cb_pos: Point3<f64>| -> Rotamer {
            let placeholder_residue_id = ResidueId::default();
            let atoms = vec![
                Atom::new("N", placeholder_residue_id, Point3::new(0.0, 1.0, 0.0)),
                Atom::new("CA", placeholder_residue_id, Point3::new(0.0, 0.0, 0.0)),
                Atom::new("C", placeholder_residue_id, Point3::new(1.0, 0.0, 0.0)),
                Atom::new("CB", placeholder_residue_id, cb_pos),
            ];

            let mut rotamer = Rotamer {
                atoms,
                bonds: vec![(0, 1), (1, 2), (1, 3)],
            };

            rotamer.atoms.iter_mut().for_each(|a| {
                a.force_field_type = if a.name == "CB" {
                    "C_SC".to_string()
                } else {
                    "BB".to_string()
                };
            });

            let res_name = system.residue(residue_id).unwrap().name.as_str();
            let topo = topology_registry.get(res_name).unwrap();

            parameterizer
                .parameterize_rotamer(&mut rotamer, res_name, topo)
                .unwrap();

            rotamer
        };

        let rotamer_a0 = create_rotamer(res_a_id, Point3::new(-0.5, -0.8, 0.0));
        let rotamer_b0 = create_rotamer(res_b_id, Point3::new(0.5, 0.8, 0.0));
        let rotamer_b1 = create_rotamer(res_b_id, Point3::new(-0.5, -0.8, 0.0));

        let mut rot_lib_map = HashMap::new();
        rot_lib_map.insert(ResidueType::Alanine, vec![rotamer_a0]);
        rot_lib_map.insert(ResidueType::Leucine, vec![rotamer_b0, rotamer_b1]);

        let rotamer_library = RotamerLibrary {
            rotamers: rot_lib_map,
        };

        parameterizer.parameterize_system(&mut system).unwrap();

        let mut el_cache = ELCache::new();
        el_cache.insert(res_a_id, ResidueType::Alanine, 0, Default::default());
        el_cache.insert(res_b_id, ResidueType::Leucine, 0, Default::default());
        el_cache.insert(res_b_id, ResidueType::Leucine, 1, Default::default());

        TestSetupBasic {
            system,
            res_a_id,
            res_b_id,
            el_cache,
            forcefield,
            rotamer_library,
            topology_registry,
            _temp_dir: temp_dir,
        }
    }

    fn setup_env_test_data() -> TestSetupWithEnv {
        let temp_dir = tempfile::tempdir().expect("Failed to create temp dir");
        let dir_path = temp_dir.path();

        let topology_path = dir_path.join("topology.toml");
        write_file(
            &topology_path,
            r#"
            [ALA]
            anchor_atoms = ["N", "CA", "C"]
            sidechain_atoms = ["CB"]
            "#,
        );

        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_a_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let res_b_id = system
            .add_residue(chain_id, 2, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let res_c_id = system
            .add_residue(chain_id, 3, "ALA", Some(ResidueType::Alanine))
            .unwrap();

        let add_residue_atoms = |system: &mut MolecularSystem, res_id: ResidueId, offset: f64| {
            let backbone_pos = vec![
                ("N", Point3::new(offset, 1.0, 0.0)),
                ("CA", Point3::new(offset, 0.0, 0.0)),
                ("C", Point3::new(offset + 1.0, 0.0, 0.0)),
            ];
            for (name, pos) in backbone_pos {
                let mut atom = Atom::new(name, res_id, pos);
                atom.force_field_type = "BB".to_string();
                system.add_atom_to_residue(res_id, atom).unwrap();
            }
            let mut cb_atom = Atom::new("CB", res_id, Point3::new(offset, -0.5, 1.2));
            cb_atom.force_field_type = "C_SC".to_string();
            system.add_atom_to_residue(res_id, cb_atom).unwrap();
        };

        add_residue_atoms(&mut system, res_a_id, 0.0);
        add_residue_atoms(&mut system, res_b_id, 5.0);
        add_residue_atoms(&mut system, res_c_id, 10.0);

        let mut vdw = HashMap::new();
        vdw.insert(
            "BB".to_string(),
            VdwParam::LennardJones {
                radius: 1.0,
                well_depth: 0.0,
            },
        );
        vdw.insert(
            "C_SC".to_string(),
            VdwParam::LennardJones {
                radius: 3.8,
                well_depth: 0.1,
            },
        );
        let forcefield = Forcefield {
            non_bonded: NonBondedParams {
                globals: GlobalParams {
                    dielectric_constant: 1.0,
                    potential_function: "lennard-jones-12-6".to_string(),
                },
                vdw,
                hbond: HashMap::new(),
                hbond_donors: HashSet::new(),
                hbond_acceptors: HashSet::new(),
            },
            deltas: HashMap::new(),
            weight_map: HashMap::new(),
        };

        let topology_registry = TopologyRegistry::load(&topology_path).unwrap();
        let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 0.0);

        let create_rotamer = |cb_pos: Point3<f64>| -> Rotamer {
            let placeholder_residue_id = ResidueId::default();
            let mut atoms = vec![
                Atom::new("N", placeholder_residue_id, Point3::new(0.0, 1.0, 0.0)),
                Atom::new("CA", placeholder_residue_id, Point3::new(0.0, 0.0, 0.0)),
                Atom::new("C", placeholder_residue_id, Point3::new(1.0, 0.0, 0.0)),
                Atom::new("CB", placeholder_residue_id, cb_pos),
            ];
            atoms.iter_mut().for_each(|a| {
                a.force_field_type = if a.name == "CB" {
                    "C_SC".to_string()
                } else {
                    "BB".to_string()
                }
            });
            let mut rotamer = Rotamer {
                atoms,
                bonds: vec![(0, 1), (1, 2), (1, 3)],
            };
            let topo = topology_registry.get("ALA").unwrap();
            parameterizer
                .parameterize_rotamer(&mut rotamer, "ALA", topo)
                .unwrap();
            rotamer
        };

        let rotamer0 = create_rotamer(Point3::new(0.0, 0.0, 0.0));
        let rotamer1 = create_rotamer(Point3::new(-0.5, -0.8, 0.0));

        let mut rot_lib_map = HashMap::new();
        rot_lib_map.insert(ResidueType::Alanine, vec![rotamer0, rotamer1]);
        let rotamer_library = RotamerLibrary {
            rotamers: rot_lib_map,
        };

        parameterizer.parameterize_system(&mut system).unwrap();

        let res_c_cb_id = system
            .residue(res_c_id)
            .unwrap()
            .get_first_atom_id_by_name("CB")
            .unwrap();
        system.atom_mut(res_c_cb_id).unwrap().position = Point3::new(5.05, 0.0, 0.0);

        let mut el_cache = ELCache::new();
        el_cache.insert(res_a_id, ResidueType::Alanine, 0, Default::default());
        el_cache.insert(res_a_id, ResidueType::Alanine, 1, Default::default());
        el_cache.insert(res_b_id, ResidueType::Alanine, 0, Default::default());
        el_cache.insert(res_b_id, ResidueType::Alanine, 1, Default::default());
        el_cache.insert(res_c_id, ResidueType::Alanine, 0, Default::default());

        TestSetupWithEnv {
            system,
            res_a_id,
            res_b_id,
            res_c_id,
            el_cache,
            forcefield,
            rotamer_library,
            topology_registry,
            _temp_dir: temp_dir,
        }
    }

    #[test]
    fn run_finds_optimal_pair_with_less_clash() {
        let setup = setup_basic_test_data();

        let config = PlacementConfigBuilder::new()
            .forcefield_path("dummy.ff")
            .delta_params_path("dummy.delta")
            .s_factor(0.0)
            .rotamer_library_path("dummy.rotlib")
            .topology_registry_path("dummy.topo")
            .max_iterations(1)
            .final_refinement_iterations(0)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.1,
                patience_iterations: 1,
            })
            .num_solutions(1)
            .residues_to_optimize(ResidueSelection::All)
            .build()
            .unwrap();

        let reporter = ProgressReporter::new();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let active_residues: HashSet<ResidueId> =
            [setup.res_a_id, setup.res_b_id].into_iter().collect();
        let result = run(
            setup.res_a_id,
            setup.res_b_id,
            &setup.system,
            &setup.el_cache,
            &context,
            &active_residues,
        )
        .unwrap();

        assert_eq!(result.rotamer_idx_a, 0);
        assert_eq!(
            result.rotamer_idx_b, 0,
            "Expected LEU rotamer 0 to be selected as it results in less steric clash"
        );
    }

    #[test]
    fn run_handles_empty_rotamer_list() {
        let mut setup = setup_basic_test_data();
        setup
            .rotamer_library
            .rotamers
            .get_mut(&ResidueType::Alanine)
            .unwrap()
            .clear();

        let config = PlacementConfigBuilder::new()
            .forcefield_path("dummy.ff")
            .delta_params_path("dummy.delta")
            .s_factor(0.0)
            .rotamer_library_path("dummy.rotlib")
            .topology_registry_path("dummy.topo")
            .max_iterations(1)
            .final_refinement_iterations(0)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.1,
                patience_iterations: 1,
            })
            .num_solutions(1)
            .residues_to_optimize(ResidueSelection::All)
            .build()
            .unwrap();

        let reporter = ProgressReporter::default();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let active_residues = HashSet::new();
        let result = run(
            setup.res_a_id,
            setup.res_b_id,
            &setup.system,
            &setup.el_cache,
            &context,
            &active_residues,
        );

        assert!(matches!(result, Err(EngineError::RotamerLibrary { .. })));
    }

    #[test]
    fn run_considers_other_active_residues_in_decision() {
        let setup = setup_env_test_data();

        let config = PlacementConfigBuilder::new()
            .forcefield_path("dummy.ff")
            .delta_params_path("dummy.delta")
            .s_factor(0.0)
            .rotamer_library_path("dummy.rotlib")
            .topology_registry_path("dummy.topo")
            .max_iterations(1)
            .final_refinement_iterations(0)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.1,
                patience_iterations: 1,
            })
            .num_solutions(1)
            .residues_to_optimize(ResidueSelection::All)
            .build()
            .unwrap();

        let reporter = ProgressReporter::new();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let active_residues: HashSet<ResidueId> = [setup.res_a_id, setup.res_b_id, setup.res_c_id]
            .into_iter()
            .collect();

        let result = run(
            setup.res_a_id,
            setup.res_b_id,
            &setup.system,
            &setup.el_cache,
            &context,
            &active_residues,
        )
        .unwrap();

        assert!(
            result.rotamer_idx_b == 1,
            "Expected residue B to choose rotamer 1 to minimize clash with residue C (got {})",
            result.rotamer_idx_b
        );
    }
}
