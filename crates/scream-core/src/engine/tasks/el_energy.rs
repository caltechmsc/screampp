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
        total: work_list.len() as u64,
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
            if let Some(native_type) = residue.residue_type {
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

    let residue_name = unit.residue_type.to_three_letter();
    let topology = context.topology_registry.get(residue_name).ok_or_else(|| {
        EngineError::TopologyNotFound {
            residue_name: residue_name.to_string(),
        }
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

        place_rotamer_on_system(&mut temp_system, unit.residue_id, rotamer, topology)?;

        let query_atoms: Vec<AtomId> = temp_system
            .residue(unit.residue_id)
            .unwrap()
            .atoms()
            .to_vec();

        let scorer = Scorer::new(&temp_system, context.forcefield);
        let energy = scorer.score_interaction(&query_atoms, &environment_atoms)?;

        energy_map.insert(rotamer_idx, energy);
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
        forcefield::{parameterization::Parameterizer, params::Forcefield},
        models::{atom::Atom, chain::ChainType, system::MolecularSystem},
        rotamers::library::RotamerLibrary,
        topology::registry::TopologyRegistry,
    };
    use crate::engine::{
        config::{
            ConvergenceConfig, DesignConfig, DesignConfigBuilder, DesignSpec, PlacementConfig,
            PlacementConfigBuilder, ResidueSelection, ResidueSpecifier,
        },
        context::OptimizationContext,
        progress::ProgressReporter,
    };
    use nalgebra::Point3;
    use std::{fs, path::Path};
    use tempfile::TempDir;

    struct TestSetup {
        system: MolecularSystem,
        forcefield: Forcefield,
        rotamer_library: RotamerLibrary,
        topology_registry: TopologyRegistry,
        _temp_dir: TempDir,
    }

    fn write_file(path: &Path, content: &str) {
        fs::write(path, content).expect("Failed to write temporary file for test");
    }

    fn setup() -> TestSetup {
        let temp_dir = tempfile::tempdir().expect("Failed to create temp dir");
        let dir = temp_dir.path();

        let ff_path = dir.join("ff.toml");
        write_file(
            &ff_path,
            r#"
[globals]
dielectric_constant = 1.0
potential_function = "lennard-jones-12-6"
[vdw]
C_BB = { radius = 1.0, well_depth = 1.0 }
N_BB = { radius = 1.0, well_depth = 1.0 }
C_SC = { radius = 1.0, well_depth = 1.0 }
"#,
        );

        let delta_path = dir.join("delta.csv");
        write_file(&delta_path, "residue_type,atom_name,mu,sigma\n");

        let topo_path = dir.join("topo.toml");
        write_file(
            &topo_path,
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
"#,
        );

        let rot_path = dir.join("rot.toml");
        write_file(
            &rot_path,
            r#"
# Minimal rotamer library for testing. Atoms are pre-parameterized.
[[ALA]]
atoms = [ { serial = 1, atom_name = "CB", position = [0.0, -1.0, -1.0], partial_charge = 0.0, force_field_type = "C_SC" } ]
bonds = []
[[LEU]]
atoms = [ { serial = 1, atom_name = "CB", position = [0.0, -1.0, -1.0], partial_charge = 0.0, force_field_type = "C_SC" } ]
bonds = []
"#,
        );

        let forcefield = Forcefield::load(&ff_path, &delta_path).unwrap();
        let topology_registry = TopologyRegistry::load(&topo_path).unwrap();

        let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 0.0);
        let rotamer_library =
            RotamerLibrary::load(&rot_path, &topology_registry, &forcefield, 0.0).unwrap();

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
        parameterizer.parameterize_system(&mut system).unwrap();

        TestSetup {
            system,
            forcefield,
            rotamer_library,
            topology_registry,
            _temp_dir: temp_dir,
        }
    }

    fn create_test_placement_config(selection: ResidueSelection) -> PlacementConfig {
        PlacementConfigBuilder::new()
            .forcefield_path("")
            .delta_params_path("")
            .rotamer_library_path("")
            .topology_registry_path("")
            .s_factor(0.0)
            .max_iterations(1)
            .num_solutions(1)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.1,
                patience_iterations: 1,
            })
            .residues_to_optimize(selection)
            .build()
            .unwrap()
    }

    fn create_test_design_config(spec: DesignSpec, repack: ResidueSelection) -> DesignConfig {
        DesignConfigBuilder::new()
            .forcefield_path("")
            .delta_params_path("")
            .rotamer_library_path("")
            .topology_registry_path("")
            .s_factor(0.0)
            .max_iterations(1)
            .num_solutions(1)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.1,
                patience_iterations: 1,
            })
            .design_spec(spec)
            .neighbors_to_repack(repack)
            .build()
            .unwrap()
    }

    #[test]
    fn build_work_list_for_placement_works() {
        let setup = setup();
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
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let work_list = build_work_list(&context).unwrap();

        assert_eq!(work_list.len(), 1);
        assert_eq!(work_list[0].residue_type, ResidueType::Alanine);
    }

    #[test]
    fn build_work_list_for_design_works() {
        let setup = setup();
        let reporter = ProgressReporter::default();
        let system = &setup.system;

        let mut design_spec = DesignSpec::new();
        design_spec.insert(
            ResidueSpecifier {
                chain_id: 'A',
                residue_number: 2,
            },
            vec![ResidueType::Alanine],
        );

        let config = create_test_design_config(
            design_spec,
            ResidueSelection::List {
                include: vec![ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: 3,
                }],
                exclude: vec![],
            },
        );
        let context = OptimizationContext::new(
            system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let work_list = build_work_list(&context).unwrap();

        assert_eq!(work_list.len(), 2);
        let has_design_site = work_list.iter().any(|w| {
            w.residue_type == ResidueType::Alanine
                && w.residue_id
                    == system
                        .find_residue_by_id(system.find_chain_by_id('A').unwrap(), 2)
                        .unwrap()
        });
        let has_repack_site = work_list.iter().any(|w| {
            w.residue_type == ResidueType::Leucine
                && w.residue_id
                    == system
                        .find_residue_by_id(system.find_chain_by_id('A').unwrap(), 3)
                        .unwrap()
        });

        assert!(has_design_site, "Work list should contain the design site");
        assert!(has_repack_site, "Work list should contain the repack site");
    }

    #[test]
    fn run_with_simple_config_succeeds() {
        let mut setup = setup();
        let chain_a = setup.system.find_chain_by_id('A').unwrap();
        let res1_id = setup.system.find_residue_by_id(chain_a, 1).unwrap();
        let mut n1 = Atom::new("N", res1_id, Point3::new(10.0, 0.0, 0.0));
        n1.force_field_type = "N_BB".to_string();
        setup.system.add_atom_to_residue(res1_id, n1).unwrap();

        let reporter = ProgressReporter::default();
        let config = create_test_placement_config(ResidueSelection::All);

        let parameterizer = Parameterizer::new(&setup.forcefield, &setup.topology_registry, 0.0);
        parameterizer
            .parameterize_system(&mut setup.system)
            .unwrap();

        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let cache = run(&context).unwrap();

        assert!(!cache.is_empty());
        assert_eq!(
            cache.len(),
            2,
            "Should have EL energies for ALA and LEU, but not GLY"
        );
    }

    #[test]
    fn run_with_empty_work_list_returns_empty_cache() {
        let setup = setup();
        let reporter = ProgressReporter::default();
        let config = create_test_placement_config(ResidueSelection::List {
            include: vec![],
            exclude: vec![],
        });
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let result = run(&context).unwrap();
        assert!(result.is_empty());
    }
}
