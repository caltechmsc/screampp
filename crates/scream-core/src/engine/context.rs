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

/// Provides a unified context for optimization operations in SCREAM++.
///
/// This struct encapsulates all the data and configuration needed to perform
/// side-chain placement and protein design operations. It provides convenient
/// access to the molecular system, force field, rotamer library, and other
/// essential components while ensuring thread safety through the Sync bound.
#[derive(Clone, Copy)]
pub struct OptimizationContext<'a, C>
where
    C: ProvidesResidueSelections + Sync,
{
    /// Reference to the molecular system being optimized.
    pub system: &'a MolecularSystem,
    /// Reference to the force field parameters for energy calculations.
    pub forcefield: &'a Forcefield,
    /// Progress reporter for tracking optimization progress.
    pub reporter: &'a ProgressReporter<'a>,
    /// Configuration object providing residue selection criteria.
    pub config: &'a C,
    /// Reference to the rotamer library for sampling conformations.
    pub rotamer_library: &'a RotamerLibrary,
    /// Reference to the topology registry for molecular structure definitions.
    pub topology_registry: &'a TopologyRegistry,
}

impl<'a, C> OptimizationContext<'a, C>
where
    C: ProvidesResidueSelections + Sync,
{
    /// Creates a new optimization context with the provided components.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system to optimize.
    /// * `forcefield` - The force field for energy calculations.
    /// * `reporter` - Progress reporter for tracking operations.
    /// * `config` - Configuration providing residue selections.
    /// * `rotamer_library` - Library of rotamer conformations.
    /// * `topology_registry` - Registry of molecular topologies.
    ///
    /// # Return
    ///
    /// Returns a new `OptimizationContext` instance.
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

    /// Resolves the repack residue selection to a set of residue IDs.
    ///
    /// This method converts the residue selection criteria from the configuration
    /// into actual residue IDs present in the molecular system, filtering for
    /// residues that have available rotamers.
    ///
    /// # Return
    ///
    /// Returns a `HashSet` of `ResidueId`s representing residues to repack.
    ///
    /// # Errors
    ///
    /// Returns `EngineError` if the selection cannot be resolved or if there
    /// are issues with the molecular system structure.
    pub fn resolve_repack_residues(&self) -> Result<HashSet<ResidueId>, EngineError> {
        let selection = self.config.repack_selection();
        query::resolve_selection_to_ids(self.system, selection, self.rotamer_library)
    }

    /// Resolves the design residue selection to a set of residue IDs.
    ///
    /// For design operations, this method extracts residue IDs from the design
    /// specification. If no design spec is provided, returns an empty set.
    ///
    /// # Return
    ///
    /// Returns a `HashSet` of `ResidueId`s representing residues to design.
    ///
    /// # Errors
    ///
    /// Returns `EngineError` if the selection cannot be resolved or if there
    /// are issues with the molecular system structure.
    pub fn resolve_design_residues(&self) -> Result<HashSet<ResidueId>, EngineError> {
        match self.config.design_spec() {
            Some(spec) => {
                // Extract residue specifiers from the design specification keys
                // to create a selection for residues that will be mutated
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

    /// Resolves all active residues for the optimization operation.
    ///
    /// This combines both repack and design residues into a single set,
    /// representing all residues that will be actively optimized.
    ///
    /// # Return
    ///
    /// Returns a `HashSet` of `ResidueId`s representing all active residues.
    ///
    /// # Errors
    ///
    /// Returns `EngineError` if any residue selection cannot be resolved.
    pub fn resolve_all_active_residues(&self) -> Result<HashSet<ResidueId>, EngineError> {
        let repack_ids = self.resolve_repack_residues()?;
        let design_ids = self.resolve_design_residues()?;
        Ok(repack_ids.union(&design_ids).cloned().collect())
    }
}

/// Provides a unified context for analysis operations in SCREAM++.
///
/// This struct encapsulates the necessary components for performing molecular
/// analysis operations such as interaction energy calculations and steric
/// clash detection. It provides access to the molecular system, force field,
/// and analysis configuration.
#[derive(Clone, Copy)]
pub struct AnalysisContext<'a> {
    /// Reference to the molecular system being analyzed.
    pub system: &'a MolecularSystem,
    /// Reference to the force field parameters for energy calculations.
    pub forcefield: &'a Forcefield,
    /// Progress reporter for tracking analysis progress.
    pub reporter: &'a ProgressReporter<'a>,
    /// Configuration object specifying the analysis type and parameters.
    pub config: &'a AnalyzeConfig,
    /// Reference to the topology registry for molecular structure definitions.
    pub topology_registry: &'a TopologyRegistry,
}

impl<'a> AnalysisContext<'a> {
    /// Creates a new analysis context with the provided components.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system to analyze.
    /// * `forcefield` - The force field for energy calculations.
    /// * `reporter` - Progress reporter for tracking operations.
    /// * `config` - Configuration specifying the analysis type.
    /// * `topology_registry` - Registry of molecular topologies.
    ///
    /// # Return
    ///
    /// Returns a new `AnalysisContext` instance.
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

/// Defines the interface for configuration objects that provide residue selections.
///
/// This trait is implemented by configuration structs that need to specify
/// which residues should be included in optimization operations. It provides
/// methods to access repack selections and optional design specifications.
pub trait ProvidesResidueSelections {
    /// Returns the residue selection for repacking operations.
    ///
    /// # Return
    ///
    /// Returns a reference to the `ResidueSelection` specifying which residues
    /// should be repacked during optimization.
    fn repack_selection(&self) -> &ResidueSelection;
    /// Returns the design specification if available.
    ///
    /// # Return
    ///
    /// Returns `Some(&DesignSpec)` if the configuration includes design operations,
    /// or `None` if only repacking is performed.
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
