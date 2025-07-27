use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
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

#[instrument(skip_all, name = "el_energy_task")]
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
        return Ok(ELCache::new());
    }

    context.reporter.report(Progress::TaskStart {
        total_steps: work_list.len() as u64,
    });

    #[cfg(not(feature = "parallel"))]
    let iterator = work_list.iter();

    #[cfg(feature = "parallel")]
    let iterator = work_list.par_iter();

    let results: Vec<WorkResult> = iterator
        .map(|unit| compute_energies_for_unit(unit, context))
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
        cached_combinations = cache.len(),
        "EL pre-computation finished."
    );
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
                for &residue_type in allowed_types {
                    work_list.push(WorkUnit {
                        residue_id,
                        residue_type,
                    });
                }
                is_design_site = true;
            }
        }

        if !is_design_site {
            if let Some(native_type) = residue.res_type {
                work_list.push(WorkUnit {
                    residue_id,
                    residue_type: native_type,
                });
            }
        }
    }
    Ok(work_list)
}

#[instrument(skip_all, fields(residue_id = ?unit.residue_id, residue_type = %unit.residue_type))]
fn compute_energies_for_unit<C>(unit: &WorkUnit, context: &OptimizationContext<C>) -> WorkResult
where
    C: ProvidesResidueSelections + Sync,
{
    let rotamers = context
        .rotamer_library
        .get_rotamers_for(unit.residue_type)
        .ok_or_else(|| EngineError::RotamerLibrary {
            residue_type: unit.residue_type.to_string(),
            message: "No rotamers found for this residue type.".to_string(),
        })?;

    let placement_info = context
        .rotamer_library
        .get_placement_info_for(unit.residue_type)
        .ok_or_else(|| EngineError::RotamerLibrary {
            residue_type: unit.residue_type.to_string(),
            message: "No placement info found for this residue type.".to_string(),
        })?;

    let active_residue_ids = context.resolve_all_active_residues()?;

    let environment_atoms: Vec<AtomId> = context
        .system
        .atoms_iter()
        .filter_map(|(atom_id, atom)| {
            if !active_residue_ids.contains(&atom.residue_id) {
                Some(atom_id)
            } else {
                None
            }
        })
        .collect();

    let mut energy_map = HashMap::with_capacity(rotamers.len());

    for (rotamer_idx, rotamer) in rotamers.iter().enumerate() {
        let mut temp_system = context.system.clone();

        place_rotamer_on_system(&mut temp_system, unit.residue_id, rotamer, placement_info)?;

        let query_atoms: Vec<AtomId> = temp_system
            .residue(unit.residue_id)
            .unwrap()
            .atoms()
            .to_vec();

        let scorer = Scorer::new(&temp_system, context.forcefield);
        let energy = scorer.score_interaction(&query_atoms, &environment_atoms)?;

        energy_map.insert(rotamer_idx, energy);
    }

    context.reporter.report(Progress::TaskIncrement);

    Ok(((unit.residue_id, unit.residue_type), energy_map))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::params::{Forcefield, GlobalParams, NonBondedParams, VdwParam};
    use crate::core::models::atom::Atom;
    use crate::core::models::chain::ChainType;
    use crate::core::models::system::MolecularSystem;
    use crate::core::rotamers::library::RotamerLibrary;
    use crate::core::rotamers::placement::PlacementInfo;
    use crate::core::rotamers::rotamer::Rotamer;
    use crate::engine::config::ConvergenceConfig;
    use crate::engine::config::{
        DesignConfig, DesignConfigBuilder, DesignSpec, PlacementConfig, PlacementConfigBuilder,
        ResidueSelection, ResidueSpecifier,
    };
    use crate::engine::context::OptimizationContext;
    use crate::engine::progress::ProgressReporter;
    use nalgebra::Point3;
    use std::collections::HashMap;

    fn create_test_system() -> MolecularSystem {
        let mut system = MolecularSystem::new();
        let chain_a = system.add_chain('A', ChainType::Protein);
        system
            .add_residue(chain_a, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        system
            .add_residue(chain_a, 2, "GLY", Some(ResidueType::Glycine))
            .unwrap();
        system
            .add_residue(chain_a, 3, "LEU", Some(ResidueType::Leucine))
            .unwrap();
        system
    }

    fn create_test_forcefield() -> Forcefield {
        let mut vdw = HashMap::new();
        vdw.insert(
            "C_BB".to_string(),
            VdwParam::LennardJones {
                radius: 1.0,
                well_depth: 1.0,
            },
        );
        vdw.insert(
            "N_BB".to_string(),
            VdwParam::LennardJones {
                radius: 1.0,
                well_depth: 1.0,
            },
        );
        vdw.insert(
            "C_SC".to_string(),
            VdwParam::LennardJones {
                radius: 1.0,
                well_depth: 1.0,
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

    fn create_test_rotamer_library() -> RotamerLibrary {
        let mut atoms = Vec::new();
        let mut n = Atom::new(1, "N", ResidueId::default(), Point3::new(0.0, 1.0, 0.0));
        n.force_field_type = "N_BB".to_string();
        atoms.push(n);
        let mut ca = Atom::new(2, "CA", ResidueId::default(), Point3::new(0.0, 0.0, 0.0));
        ca.force_field_type = "C_BB".to_string();
        atoms.push(ca);
        let mut c = Atom::new(3, "C", ResidueId::default(), Point3::new(1.0, 0.0, 0.0));
        c.force_field_type = "C_BB".to_string();
        atoms.push(c);
        let mut cb = Atom::new(4, "CB", ResidueId::default(), Point3::new(0.0, -1.0, -1.0));
        cb.force_field_type = "C_SC".to_string();
        atoms.push(cb);

        let ala_rotamer = Rotamer { atoms };

        let mut rotamers = HashMap::new();
        rotamers.insert(ResidueType::Alanine, vec![ala_rotamer.clone()]);
        rotamers.insert(ResidueType::Leucine, vec![ala_rotamer]);

        let ala_placement = PlacementInfo {
            anchor_atoms: vec!["N".to_string(), "CA".to_string(), "C".to_string()],
            sidechain_atoms: vec!["CB".to_string()],
            exact_match_atoms: vec![],
            connection_points: vec![],
        };
        let mut placement_info = HashMap::new();
        placement_info.insert(ResidueType::Alanine, ala_placement.clone());
        placement_info.insert(ResidueType::Leucine, ala_placement);

        RotamerLibrary {
            rotamers,
            placement_info,
        }
    }

    fn create_test_placement_config(selection: ResidueSelection) -> PlacementConfig {
        PlacementConfigBuilder::new()
            .forcefield_path("")
            .delta_params_path("")
            .s_factor(0.0)
            .rotamer_library_path("")
            .placement_registry_path("")
            .max_iterations(1)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.1,
                patience_iterations: 1,
            })
            .num_solutions(1)
            .residues_to_optimize(selection)
            .build()
            .unwrap()
    }

    fn create_test_design_config(
        design_spec: DesignSpec,
        repack_selection: ResidueSelection,
    ) -> DesignConfig {
        DesignConfigBuilder::new()
            .forcefield_path("")
            .delta_params_path("")
            .s_factor(0.0)
            .rotamer_library_path("")
            .placement_registry_path("")
            .max_iterations(1)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.1,
                patience_iterations: 1,
            })
            .num_solutions(1)
            .design_spec(design_spec)
            .neighbors_to_repack(repack_selection)
            .build()
            .unwrap()
    }

    #[test]
    fn build_work_list_for_placement_works() {
        let system = create_test_system();
        let ff = create_test_forcefield();
        let rot_lib = create_test_rotamer_library();
        let reporter = ProgressReporter::default();
        let selection = ResidueSelection::List {
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
        };
        let config = create_test_placement_config(selection);
        let context = OptimizationContext::new(&system, &ff, &reporter, &config, &rot_lib);

        let work_list = build_work_list(&context).unwrap();

        assert_eq!(work_list.len(), 1);
        assert_eq!(work_list[0].residue_type, ResidueType::Alanine);
    }

    #[test]
    fn build_work_list_for_design_works() {
        let system = create_test_system();
        let ff = create_test_forcefield();
        let rot_lib = create_test_rotamer_library();
        let reporter = ProgressReporter::default();

        let mut design_spec = DesignSpec::new();
        design_spec.insert(
            ResidueSpecifier {
                chain_id: 'A',
                residue_number: 1,
            },
            vec![ResidueType::Leucine],
        );

        let config = create_test_design_config(design_spec, ResidueSelection::All);
        let context = OptimizationContext::new(&system, &ff, &reporter, &config, &rot_lib);

        let work_list = build_work_list(&context).unwrap();

        assert_eq!(work_list.len(), 2);
        let has_design_site = work_list.iter().any(|w| {
            w.residue_type == ResidueType::Leucine
                && w.residue_id
                    == system
                        .find_residue_by_id(system.find_chain_by_id('A').unwrap(), 1)
                        .unwrap()
        });
        let has_repack_site = work_list.iter().any(|w| {
            w.residue_type == ResidueType::Leucine
                && w.residue_id
                    == system
                        .find_residue_by_id(system.find_chain_by_id('A').unwrap(), 3)
                        .unwrap()
        });

        assert!(has_design_site);
        assert!(has_repack_site);
    }

    #[test]
    fn run_with_simple_config_succeeds() {
        let mut system = create_test_system();
        let chain_a = system.find_chain_by_id('A').unwrap();

        let res1_id = system.find_residue_by_id(chain_a, 1).unwrap();
        let mut n1 = Atom::new(1, "N", res1_id, Point3::new(0.0, 1.0, 0.0));
        n1.force_field_type = "N_BB".to_string();
        system.add_atom_to_residue(res1_id, n1).unwrap();
        let mut ca1 = Atom::new(2, "CA", res1_id, Point3::new(0.0, 0.0, 0.0));
        ca1.force_field_type = "C_BB".to_string();
        system.add_atom_to_residue(res1_id, ca1).unwrap();
        let mut c1 = Atom::new(3, "C", res1_id, Point3::new(1.0, 0.0, 0.0));
        c1.force_field_type = "C_BB".to_string();
        system.add_atom_to_residue(res1_id, c1).unwrap();

        let res3_id = system.find_residue_by_id(chain_a, 3).unwrap();
        let mut n3 = Atom::new(4, "N", res3_id, Point3::new(2.0, 1.0, 0.0));
        n3.force_field_type = "N_BB".to_string();
        system.add_atom_to_residue(res3_id, n3).unwrap();
        let mut ca3 = Atom::new(5, "CA", res3_id, Point3::new(2.0, 0.0, 0.0));
        ca3.force_field_type = "C_BB".to_string();
        system.add_atom_to_residue(res3_id, ca3).unwrap();
        let mut c3 = Atom::new(6, "C", res3_id, Point3::new(3.0, 0.0, 0.0));
        c3.force_field_type = "C_BB".to_string();
        system.add_atom_to_residue(res3_id, c3).unwrap();

        let ff = create_test_forcefield();
        let rot_lib = create_test_rotamer_library();
        let reporter = ProgressReporter::default();
        let config = create_test_placement_config(ResidueSelection::All);
        let context = OptimizationContext::new(&system, &ff, &reporter, &config, &rot_lib);

        let result = run(&context);
        let cache = result.unwrap();

        assert!(!cache.is_empty());
        assert_eq!(cache.len(), 2);
    }

    #[test]
    fn run_with_empty_work_list_returns_empty_cache() {
        let system = create_test_system();
        let ff = create_test_forcefield();
        let rot_lib = create_test_rotamer_library();
        let reporter = ProgressReporter::default();
        let config = create_test_placement_config(ResidueSelection::List {
            include: vec![],
            exclude: vec![],
        });
        let context = OptimizationContext::new(&system, &ff, &reporter, &config, &rot_lib);

        let result = run(&context).unwrap();
        assert!(result.is_empty());
    }
}
