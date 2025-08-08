use crate::core::forcefield::scoring::Scorer;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::core::models::system::MolecularSystem;
use crate::engine::cache::ELCache;
use crate::engine::context::{OptimizationContext, ProvidesResidueSelections};
use crate::engine::error::EngineError;
use crate::engine::placement;
use tracing::{debug, instrument, trace};

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
) -> Result<DoubletResult, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    let residue_a = system.residue(res_a_id).unwrap();
    let residue_b = system.residue(res_b_id).unwrap();
    let res_type_a = residue_a.residue_type.unwrap();
    let res_type_b = residue_b.residue_type.unwrap();

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

    let res_name_a = res_type_a.to_three_letter();
    let topology_a =
        context
            .topology_registry
            .get(res_name_a)
            .ok_or_else(|| EngineError::TopologyNotFound {
                residue_name: res_name_a.to_string(),
            })?;

    let res_name_b = res_type_b.to_three_letter();
    let topology_b =
        context
            .topology_registry
            .get(res_name_b)
            .ok_or_else(|| EngineError::TopologyNotFound {
                residue_name: res_name_b.to_string(),
            })?;

    #[cfg(not(feature = "parallel"))]
    let iterator = index_pairs.iter();

    #[cfg(feature = "parallel")]
    let iterator = index_pairs.par_iter();

    let best_pair = iterator
        .map(|&(idx_a, idx_b)| {
            let rot_a = &rotamers_a[idx_a];
            let rot_b = &rotamers_b[idx_b];

            let el_a_term = el_cache
                .get(res_a_id, res_type_a, idx_a)
                .copied()
                .unwrap_or_default();
            let el_b_term = el_cache
                .get(res_b_id, res_type_b, idx_b)
                .copied()
                .unwrap_or_default();
            let el_a_total = el_a_term.total();
            let el_b_total = el_b_term.total();

            let mut temp_system = system.clone();
            placement::place_rotamer_on_system(&mut temp_system, res_a_id, rot_a, topology_a)?;
            placement::place_rotamer_on_system(&mut temp_system, res_b_id, rot_b, topology_b)?;

            let atoms_a_all_ids: Vec<AtomId> =
                temp_system.residue(res_a_id).unwrap().atoms().to_vec();
            let atoms_b_all_ids: Vec<AtomId> =
                temp_system.residue(res_b_id).unwrap().atoms().to_vec();

            let scorer = Scorer::new(&temp_system, context.forcefield);
            let interaction_term = scorer.score_interaction(&atoms_a_all_ids, &atoms_b_all_ids)?;
            let interaction_total = interaction_term.total();
            let local_energy = el_a_total + el_b_total + interaction_total;

            trace!(
                "Pair (A:{}, B:{}): E_local = {:.2} (EL_A={:.2}, EL_B={:.2}, Int={:.2})",
                idx_a, idx_b, local_energy, el_a_total, el_b_total, interaction_total
            );

            Ok((local_energy, (idx_a, idx_b)))
        })
        .collect::<Result<Vec<_>, EngineError>>()?
        .into_iter()
        .min_by(|(energy_1, _), (energy_2, _)| {
            energy_1
                .partial_cmp(energy_2)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

    match best_pair {
        Some((best_energy, (best_idx_a, best_idx_b))) => Ok(DoubletResult {
            rotamer_idx_a: best_idx_a,
            rotamer_idx_b: best_idx_b,
            best_local_energy: best_energy,
        }),
        None => Err(EngineError::PhaseFailed {
            phase: "Doublet Optimization",
            reason: "No valid rotamer pairs could be evaluated.".to_string(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::atom::CachedVdwParam;
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
    use std::{collections::HashMap, fs, path::Path};
    use tempfile::TempDir;

    struct TestSetup {
        system: MolecularSystem,
        res_a_id: ResidueId,
        res_b_id: ResidueId,
        el_cache: ELCache,
        forcefield: Forcefield,
        rotamer_library: RotamerLibrary,
        topology_registry: TopologyRegistry,
        _temp_dir: TempDir,
    }

    fn write_file(path: &Path, content: &str) {
        fs::write(path, content).expect("Failed to write temporary file for test");
    }

    fn setup_test_data() -> TestSetup {
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
                atom.vdw_param = CachedVdwParam::LennardJones {
                    radius: 1.0,
                    well_depth: 0.0,
                };
                system.add_atom_to_residue(res_id, atom).unwrap();
            }
            let mut cb_atom = Atom::new("CB", res_id, Point3::new(offset, -0.5, 1.2));
            cb_atom.force_field_type = "C_SC".to_string();
            cb_atom.vdw_param = CachedVdwParam::LennardJones {
                radius: 3.8,
                well_depth: 0.1,
            };
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
            },
            deltas: HashMap::new(),
        };

        let topology_registry = TopologyRegistry::load(&topology_path).unwrap();

        let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 0.0);

        let create_rotamer = |residue_id, cb_pos: Point3<f64>| -> Rotamer {
            let mut atoms = Vec::new();
            let mut n = Atom::new("N", residue_id, Point3::new(0.0, 1.0, 0.0));
            n.force_field_type = "BB".to_string();
            n.vdw_param = CachedVdwParam::LennardJones {
                radius: 1.0,
                well_depth: 0.0,
            };
            atoms.push(n);
            let mut ca = Atom::new("CA", residue_id, Point3::new(0.0, 0.0, 0.0));
            ca.force_field_type = "BB".to_string();
            ca.vdw_param = CachedVdwParam::LennardJones {
                radius: 1.0,
                well_depth: 0.0,
            };
            atoms.push(ca);
            let mut c = Atom::new("C", residue_id, Point3::new(1.0, 0.0, 0.0));
            c.force_field_type = "BB".to_string();
            c.vdw_param = CachedVdwParam::LennardJones {
                radius: 1.0,
                well_depth: 0.0,
            };
            atoms.push(c);
            let mut cb = Atom::new("CB", residue_id, cb_pos);
            cb.force_field_type = "C_SC".to_string();
            cb.vdw_param = CachedVdwParam::LennardJones {
                radius: 3.8,
                well_depth: 0.1,
            };
            atoms.push(cb);

            let bonds = vec![(0, 1), (1, 2), (1, 3)];

            let mut rotamer = Rotamer { atoms, bonds };

            let res_name = system.residue(residue_id).unwrap().name.as_str();
            let topo = topology_registry.get(res_name).unwrap();
            for atom in &mut rotamer.atoms {
                parameterizer
                    .parameterize_protein_atom(atom, res_name, topo)
                    .unwrap();
            }
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

        TestSetup {
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

    #[test]
    fn run_finds_optimal_pair_with_less_clash() {
        let setup = setup_test_data();

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

        let result = run(
            setup.res_a_id,
            setup.res_b_id,
            &setup.system,
            &setup.el_cache,
            &context,
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
        let mut setup = setup_test_data();
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

        let result = run(
            setup.res_a_id,
            setup.res_b_id,
            &setup.system,
            &setup.el_cache,
            &context,
        );

        assert!(matches!(result, Err(EngineError::RotamerLibrary { .. })));
    }
}
