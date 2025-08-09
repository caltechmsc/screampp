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
            .to_vec();

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
H_ = { radius = 1.0, well_depth = 1.0 }
[hbond]
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
[[ALA]]
atoms = [
    { serial = 1, atom_name = "N", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "N_BB" },
    { serial = 2, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_BB" },
    { serial = 3, atom_name = "C", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_BB" },
    { serial = 4, atom_name = "CB", position = [0.0, -1.0, -1.0], partial_charge = 0.0, force_field_type = "C_SC" }
]
bonds = []
[[LEU]]
atoms = [
    { serial = 1, atom_name = "N", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "N_BB" },
    { serial = 2, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_BB" },
    { serial = 3, atom_name = "C", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_BB" },
    { serial = 4, atom_name = "CB", position = [0.0, -1.0, -1.0], partial_charge = 0.0, force_field_type = "C_SC" },
    { serial = 5, atom_name = "CG", position = [0.0, -2.0, -2.0], partial_charge = 0.0, force_field_type = "C_SC" }
]
bonds = []
[[GLY]]
atoms = [
    { serial = 1, atom_name = "N", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "N_BB" },
    { serial = 2, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_BB" },
    { serial = 3, atom_name = "C", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_BB" },
    { serial = 4, atom_name = "HA1", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "H_" },
    { serial = 5, atom_name = "HA2", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "H_" }
]
bonds = []
"#,
        );

        let forcefield = Forcefield::load(&ff_path, &delta_path).unwrap();
        let topology_registry = TopologyRegistry::load(&topo_path).unwrap();

        let rotamer_library =
            RotamerLibrary::load(&rot_path, &topology_registry, &forcefield, 0.0).unwrap();

        let mut system = MolecularSystem::new();
        let chain_a = system.add_chain('A', ChainType::Protein);

        let ala_id = system
            .add_residue(chain_a, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        for name in ["N", "CA", "C"].iter() {
            let mut atom = Atom::new(name, ala_id, Point3::origin());
            atom.force_field_type = "C_BB".to_string();
            system.add_atom_to_residue(ala_id, atom).unwrap();
        }
        let mut cb_atom = Atom::new("CB", ala_id, Point3::origin());
        cb_atom.force_field_type = "C_SC".to_string();
        system.add_atom_to_residue(ala_id, cb_atom).unwrap();

        let gly_id = system
            .add_residue(chain_a, 2, "GLY", Some(ResidueType::Glycine))
            .unwrap();
        for name in ["N", "CA", "C"].iter() {
            let mut atom = Atom::new(name, gly_id, Point3::origin());
            atom.force_field_type = "C_BB".to_string();
            system.add_atom_to_residue(gly_id, atom).unwrap();
        }

        let leu_id = system
            .add_residue(chain_a, 3, "LEU", Some(ResidueType::Leucine))
            .unwrap();
        for name in ["N", "CA", "C", "CB", "CG"].iter() {
            let mut atom = Atom::new(name, leu_id, Point3::origin());
            atom.force_field_type = if name.len() == 1 {
                "C_BB".to_string()
            } else {
                "C_SC".to_string()
            };
            system.add_atom_to_residue(leu_id, atom).unwrap();
        }

        let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 0.0);
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
            .final_refinement_iterations(0)
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
            .include_input_conformation(false)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.1,
                patience_iterations: 1,
            })
            .design_spec(spec)
            .neighbors_to_repack(repack)
            .final_refinement_iterations(0)
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

        assert_eq!(work_list.len(), 2);

        let has_ala = work_list
            .iter()
            .any(|w| w.residue_type == ResidueType::Alanine);
        let has_gly = work_list
            .iter()
            .any(|w| w.residue_type == ResidueType::Glycine);

        assert!(has_ala, "Work list should include Alanine");
        assert!(has_gly, "Work list should include Glycine");
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
                    residue_number: 3,
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

        let cache = run(&context).unwrap();

        assert!(!cache.is_empty());
        assert_eq!(
            cache.len(),
            2,
            "Should have EL energies for ALA and LEU, which were active."
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
