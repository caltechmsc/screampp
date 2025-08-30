use super::config::{AnalyzeConfig, DesignConfig, DesignSpec, PlacementConfig, ResidueSelection};
use super::error::EngineError;
use super::progress::ProgressReporter;
use super::utils::query;
use crate::core::{
    forcefield::params::Forcefield,
    models::{ids::ResidueId, system::MolecularSystem},
    rotamers::library::RotamerLibrary,
    topology::registry::TopologyRegistry,
};
use std::collections::HashSet;

#[derive(Clone, Copy)]
pub struct OptimizationContext<'a, C>
where
    C: ProvidesResidueSelections + Sync,
{
    pub system: &'a MolecularSystem,
    pub forcefield: &'a Forcefield,
    pub reporter: &'a ProgressReporter<'a>,
    pub config: &'a C,
    pub rotamer_library: &'a RotamerLibrary,
    pub topology_registry: &'a TopologyRegistry,
}

impl<'a, C> OptimizationContext<'a, C>
where
    C: ProvidesResidueSelections + Sync,
{
    pub fn new(
        system: &'a MolecularSystem,
        forcefield: &'a Forcefield,
        reporter: &'a ProgressReporter<'a>,
        config: &'a C,
        rotamer_library: &'a RotamerLibrary,
        topology_registry: &'a TopologyRegistry,
    ) -> Self {
        Self {
            system,
            forcefield,
            reporter,
            config,
            rotamer_library,
            topology_registry,
        }
    }

    pub fn resolve_repack_residues(&self) -> Result<HashSet<ResidueId>, EngineError> {
        let selection = self.config.repack_selection();
        query::resolve_selection_to_ids(self.system, selection, self.rotamer_library)
    }

    pub fn resolve_design_residues(&self) -> Result<HashSet<ResidueId>, EngineError> {
        match self.config.design_spec() {
            Some(spec) => {
                let specifiers: Vec<_> = spec.keys().cloned().collect();
                let selection = ResidueSelection::List {
                    include: specifiers,
                    exclude: vec![],
                };
                query::resolve_selection_to_ids(self.system, &selection, self.rotamer_library)
            }
            None => Ok(HashSet::new()),
        }
    }

    pub fn resolve_all_active_residues(&self) -> Result<HashSet<ResidueId>, EngineError> {
        let repack_ids = self.resolve_repack_residues()?;
        let design_ids = self.resolve_design_residues()?;
        Ok(repack_ids.union(&design_ids).cloned().collect())
    }
}

#[derive(Clone, Copy)]
pub struct AnalysisContext<'a> {
    pub system: &'a MolecularSystem,
    pub forcefield: &'a Forcefield,
    pub reporter: &'a ProgressReporter<'a>,
    pub config: &'a AnalyzeConfig,
    pub topology_registry: &'a TopologyRegistry,
}

impl<'a> AnalysisContext<'a> {
    pub fn new(
        system: &'a MolecularSystem,
        forcefield: &'a Forcefield,
        reporter: &'a ProgressReporter<'a>,
        config: &'a AnalyzeConfig,
        topology_registry: &'a TopologyRegistry,
    ) -> Self {
        Self {
            system,
            forcefield,
            reporter,
            config,
            topology_registry,
        }
    }
}

pub trait ProvidesResidueSelections {
    fn repack_selection(&self) -> &ResidueSelection;
    fn design_spec(&self) -> Option<&DesignSpec> {
        None
    }
}
impl ProvidesResidueSelections for PlacementConfig {
    fn repack_selection(&self) -> &ResidueSelection {
        &self.residues_to_optimize
    }
}
impl ProvidesResidueSelections for DesignConfig {
    fn repack_selection(&self) -> &ResidueSelection {
        &self.neighbors_to_repack
    }
    fn design_spec(&self) -> Option<&DesignSpec> {
        Some(&self.design_spec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{
        forcefield::params::{Forcefield, GlobalParams, NonBondedParams},
        models::{chain::ChainType, residue::ResidueType, system::MolecularSystem},
        rotamers::{library::RotamerLibrary, rotamer::Rotamer},
        topology::registry::TopologyRegistry,
    };
    use crate::engine::config::ConvergenceConfig;
    use crate::engine::{
        config::{PlacementConfig, PlacementConfigBuilder, ResidueSelection, ResidueSpecifier},
        progress::ProgressReporter,
    };
    use nalgebra::Point3;
    use std::collections::HashMap;

    struct TestSetup {
        system: MolecularSystem,
        rotamer_library: RotamerLibrary,
        topology_registry: TopologyRegistry,
        residue_ids: HashMap<String, ResidueId>,
    }

    fn setup() -> TestSetup {
        let mut system = MolecularSystem::new();
        let mut residue_ids = HashMap::new();

        let chain_a_id = system.add_chain('A', ChainType::Protein);
        let gly_id = system
            .add_residue(chain_a_id, 1, "GLY", Some(ResidueType::Glycine))
            .unwrap();
        residue_ids.insert("GLY".to_string(), gly_id);

        let ala_id = system
            .add_residue(chain_a_id, 2, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let atom_ala =
            crate::core::models::atom::Atom::new("CA", ala_id, Point3::new(5.0, 0.0, 0.0));
        system.add_atom_to_residue(ala_id, atom_ala).unwrap();
        residue_ids.insert("ALA".to_string(), ala_id);

        let lys_id = system
            .add_residue(chain_a_id, 3, "LYS", Some(ResidueType::Lysine))
            .unwrap();
        let atom_lys =
            crate::core::models::atom::Atom::new("CA", lys_id, Point3::new(10.0, 0.0, 0.0));
        system.add_atom_to_residue(lys_id, atom_lys).unwrap();
        residue_ids.insert("LYS".to_string(), lys_id);

        let chain_b_id = system.add_chain('B', ChainType::Ligand);
        let lig_id = system.add_residue(chain_b_id, 10, "LIG", None).unwrap();
        let atom_lig =
            crate::core::models::atom::Atom::new("C1", lig_id, Point3::new(0.0, 0.0, 0.0));
        system.add_atom_to_residue(lig_id, atom_lig).unwrap();
        residue_ids.insert("LIG".to_string(), lig_id);

        let mut rotamer_library = RotamerLibrary::default();
        let empty_rotamer = Rotamer {
            atoms: vec![],
            bonds: vec![],
        };
        rotamer_library
            .rotamers
            .insert(ResidueType::Alanine, vec![empty_rotamer.clone()]);
        rotamer_library
            .rotamers
            .insert(ResidueType::Lysine, vec![empty_rotamer]);

        let topology_registry = TopologyRegistry::default();

        TestSetup {
            system,
            rotamer_library,
            topology_registry,
            residue_ids,
        }
    }

    fn create_mock_config(selection: ResidueSelection) -> PlacementConfig {
        PlacementConfigBuilder::new()
            .forcefield_path("")
            .delta_params_path("")
            .rotamer_library_path("")
            .topology_registry_path("")
            .s_factor(0.0)
            .max_iterations(1)
            .num_solutions(1)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.01,
                patience_iterations: 10,
            })
            .residues_to_optimize(selection)
            .final_refinement_iterations(0)
            .build()
            .unwrap()
    }

    #[test]
    fn context_resolve_repack_residues_works() {
        let setup = setup();
        let selection = ResidueSelection::List {
            include: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 2,
            }],
            exclude: vec![],
        };
        let config = create_mock_config(selection);
        let reporter = ProgressReporter::default();
        let forcefield = Forcefield {
            non_bonded: NonBondedParams {
                globals: GlobalParams {
                    dielectric_constant: 1.0,
                    potential_function: "".to_string(),
                },
                vdw: HashMap::new(),
                hbond: HashMap::new(),
                hbond_donors: HashSet::new(),
                hbond_acceptors: HashSet::new(),
            },
            deltas: HashMap::new(),
            weight_map: HashMap::new(),
        };

        let context = OptimizationContext::new(
            &setup.system,
            &forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );

        let active_residues = context.resolve_repack_residues().unwrap();

        assert_eq!(active_residues.len(), 1);
        assert!(active_residues.contains(setup.residue_ids.get("ALA").unwrap()));
    }
}
