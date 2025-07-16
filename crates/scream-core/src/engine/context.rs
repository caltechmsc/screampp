use super::config::{DesignConfig, DesignSpec, PlacementConfig, ResidueSelection};
use super::error::EngineError;
use crate::core::forcefield::params::Forcefield;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::core::rotamers::library::RotamerLibrary;
use std::collections::HashSet;

#[derive(Clone, Copy)]
pub struct Context<'a, C> {
    pub system: &'a MolecularSystem,
    pub config: &'a C,
    pub forcefield: &'a Forcefield,
    pub rotamer_library: &'a RotamerLibrary,
}

impl<'a, C> Context<'a, C> {
    pub fn new(
        system: &'a MolecularSystem,
        config: &'a C,
        forcefield: &'a Forcefield,
        rotamer_library: &'a RotamerLibrary,
    ) -> Self {
        Self {
            system,
            config,
            forcefield,
            rotamer_library,
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

impl<'a, C: ProvidesResidueSelections> Context<'a, C> {
    pub fn resolve_repack_residues(&self) -> Result<HashSet<ResidueId>, EngineError> {
        let selection = self.config.repack_selection();
        resolve_selection_to_ids(self.system, selection, self.rotamer_library)
    }

    pub fn resolve_design_residues(&self) -> Result<HashSet<ResidueId>, EngineError> {
        match self.config.design_spec() {
            Some(spec) => {
                let specifiers: Vec<_> = spec.keys().cloned().collect();
                let selection = ResidueSelection::List {
                    include: specifiers,
                    exclude: vec![],
                };
                resolve_selection_to_ids(self.system, &selection, self.rotamer_library)
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

fn resolve_selection_to_ids(
    system: &MolecularSystem,
    selection: &ResidueSelection,
    library: &RotamerLibrary,
) -> Result<HashSet<ResidueId>, EngineError> {
    let candidate_ids: HashSet<ResidueId> = match selection {
        ResidueSelection::All => system.residues_iter().map(|(id, _)| id).collect(),
        ResidueSelection::List { include, exclude } => {
            let mut initial_set = HashSet::new();

            if include.is_empty() && !exclude.is_empty() {
                initial_set = system.residues_iter().map(|(id, _)| id).collect();
            } else {
                for spec in include {
                    let chain_id = system
                        .find_chain_by_id(spec.chain_id)
                        .ok_or_else(|| EngineError::ResidueNotFound { spec: spec.clone() })?;
                    let residue_id = system
                        .find_residue_by_id(chain_id, spec.residue_number)
                        .ok_or_else(|| EngineError::ResidueNotFound { spec: spec.clone() })?;
                    initial_set.insert(residue_id);
                }
            }

            for spec in exclude {
                if let Some(chain_id) = system.find_chain_by_id(spec.chain_id) {
                    if let Some(residue_id) =
                        system.find_residue_by_id(chain_id, spec.residue_number)
                    {
                        initial_set.remove(&residue_id);
                    }
                }
            }
            initial_set
        }
        ResidueSelection::LigandBindingSite { .. } => {
            // TODO: Implement ligand binding site selection
            unimplemented!("Ligand binding site selection is not yet implemented.");
        }
    };

    let final_active_residues = candidate_ids
        .into_iter()
        .filter(|&residue_id| {
            system
                .residue(residue_id)
                .and_then(|res| res.res_type)
                .map_or(false, |res_type| {
                    library.get_rotamers_for(res_type).is_some()
                })
        })
        .collect();

    Ok(final_active_residues)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::params::Forcefield;
    use crate::core::models::chain::ChainType;
    use crate::core::models::residue::ResidueType;
    use crate::engine::config::{DesignSpec, OptimizationConfig, ResidueSpecifier, ScoringConfig};

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
        let chain_b = system.add_chain('B', ChainType::Protein);
        system
            .add_residue(chain_b, 10, "TRP", Some(ResidueType::Tryptophan))
            .unwrap();
        system
            .add_residue(chain_b, 11, "XXX", None) // Non-standard residue
            .unwrap();
        system
    }

    fn create_test_scoring_config() -> ScoringConfig {
        ScoringConfig {
            forcefield_path: "".into(),
            rotamer_library_path: "".into(),
            delta_params_path: "".into(),
            s_factor: 1.0,
        }
    }

    fn create_test_optimization_config() -> OptimizationConfig {
        OptimizationConfig {
            max_iterations: 100,
            convergence_threshold: 0.01,
            num_solutions: 1,
        }
    }

    fn create_test_forcefield() -> Forcefield {
        use crate::core::forcefield::params::{GlobalParams, NonBondedParams};
        use std::collections::HashMap;

        let globals = GlobalParams {
            dielectric_constant: 1.0,
            potential_function: "mixed".to_string(),
        };

        let non_bonded = NonBondedParams {
            globals,
            vdw: HashMap::new(),
            hbond: HashMap::new(),
        };

        Forcefield {
            non_bonded,
            deltas: HashMap::new(),
        }
    }

    fn create_test_rotamer_library() -> RotamerLibrary {
        let mut lib = RotamerLibrary::default();
        lib.rotamers
            .insert(ResidueType::Alanine, Default::default());
        lib.rotamers
            .insert(ResidueType::Leucine, Default::default());
        lib
    }

    #[test]
    fn context_creation_is_correct() {
        let system = MolecularSystem::new();
        let ff = create_test_forcefield();
        let rot_lib = create_test_rotamer_library();
        let config = ();
        let context = Context::new(&system, &config, &ff, &rot_lib);
        assert!(std::ptr::eq(context.system, &system));
        assert!(std::ptr::eq(context.forcefield, &ff));
        assert!(std::ptr::eq(context.config, &config));
        assert!(std::ptr::eq(context.rotamer_library, &rot_lib));
    }

    #[test]
    fn resolve_all_residues_selects_all_supported_except_glycine() {
        let system = create_test_system();
        let library = create_test_rotamer_library();
        let selection = ResidueSelection::All;
        let resolved = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        assert_eq!(resolved.len(), 2);

        let chain_a = system.find_chain_by_id('A').unwrap();
        let ala_id = system.find_residue_by_id(chain_a, 1).unwrap();
        let gly_id = system.find_residue_by_id(chain_a, 2).unwrap();
        let leu_id = system.find_residue_by_id(chain_a, 3).unwrap();
        let chain_b = system.find_chain_by_id('B').unwrap();
        let trp_id = system.find_residue_by_id(chain_b, 10).unwrap();
        let xxx_id = system.find_residue_by_id(chain_b, 11).unwrap();

        assert!(resolved.contains(&ala_id));
        assert!(resolved.contains(&leu_id));
        assert!(!resolved.contains(&gly_id));
        assert!(!resolved.contains(&trp_id));
        assert!(!resolved.contains(&xxx_id));
    }

    #[test]
    fn resolve_list_with_include_only_selects_specified_and_supported_residues() {
        let system = create_test_system();
        let library = create_test_rotamer_library();
        let selection = ResidueSelection::List {
            include: vec![
                ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: 1,
                },
                ResidueSpecifier {
                    chain_id: 'B',
                    residue_number: 10,
                },
            ],
            exclude: vec![],
        };
        let resolved = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        assert_eq!(resolved.len(), 1);
        let chain_a = system.find_chain_by_id('A').unwrap();
        let ala_id = system.find_residue_by_id(chain_a, 1).unwrap();
        assert!(resolved.contains(&ala_id));
    }

    #[test]
    fn resolve_list_with_exclude_only_selects_all_supported_except_excluded() {
        let system = create_test_system();
        let library = create_test_rotamer_library();
        let selection = ResidueSelection::List {
            include: vec![],
            exclude: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 1,
            }],
        };
        let resolved = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        assert_eq!(resolved.len(), 1);
        let chain_a = system.find_chain_by_id('A').unwrap();
        let ala_id = system.find_residue_by_id(chain_a, 1).unwrap();
        let leu_id = system.find_residue_by_id(chain_a, 3).unwrap();

        assert!(!resolved.contains(&ala_id));
        assert!(resolved.contains(&leu_id));
    }

    #[test]
    fn resolve_list_with_include_and_exclude_works_correctly() {
        let system = create_test_system();
        let library = create_test_rotamer_library();
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
            exclude: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 1,
            }],
        };
        let resolved = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        assert_eq!(resolved.len(), 1);
        let chain_a = system.find_chain_by_id('A').unwrap();
        let leu_id = system.find_residue_by_id(chain_a, 3).unwrap();
        assert!(resolved.contains(&leu_id));
    }

    #[test]
    fn resolve_list_with_nonexistent_include_residue_returns_error() {
        let system = create_test_system();
        let library = create_test_rotamer_library();
        let spec = ResidueSpecifier {
            chain_id: 'Z',
            residue_number: 99,
        };
        let selection = ResidueSelection::List {
            include: vec![spec.clone()],
            exclude: vec![],
        };
        let result = resolve_selection_to_ids(&system, &selection, &library);

        assert!(matches!(
            result,
            Err(EngineError::ResidueNotFound { spec: e_spec }) if e_spec == spec
        ));
    }

    #[test]
    fn resolve_list_with_nonexistent_exclude_residue_is_ignored() {
        let system = create_test_system();
        let library = create_test_rotamer_library();
        let selection = ResidueSelection::List {
            include: vec![],
            exclude: vec![ResidueSpecifier {
                chain_id: 'Z',
                residue_number: 99,
            }],
        };
        let resolved = resolve_selection_to_ids(&system, &selection, &library).unwrap();
        assert_eq!(resolved.len(), 2);
    }

    #[test]
    fn resolve_list_with_empty_include_and_exclude_returns_empty_set() {
        let system = create_test_system();
        let library = create_test_rotamer_library();
        let selection = ResidueSelection::List {
            include: vec![],
            exclude: vec![],
        };
        let resolved = resolve_selection_to_ids(&system, &selection, &library).unwrap();
        assert!(resolved.is_empty());
    }

    // TODO: Implement ligand binding site resolution and update this test
    #[test]
    #[should_panic(expected = "Ligand binding site selection is not yet implemented.")]
    fn resolve_ligand_binding_site_panics() {
        let system = create_test_system();
        let library = create_test_rotamer_library();
        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue_name: "LIG".to_string(),
            radius_angstroms: 5.0,
        };
        let _ = resolve_selection_to_ids(&system, &selection, &library);
    }

    #[test]
    fn resolve_design_residues_works_for_design_config() {
        let system = create_test_system();
        let library = create_test_rotamer_library();
        let ff = create_test_forcefield();

        let spec = ResidueSpecifier {
            chain_id: 'A',
            residue_number: 1,
        };
        let mut design_spec = DesignSpec::new();
        design_spec.insert(spec.clone(), vec![ResidueType::Leucine]);

        let config = DesignConfig {
            design_spec,
            neighbors_to_repack: ResidueSelection::List {
                include: vec![],
                exclude: vec![],
            },
            scoring: create_test_scoring_config(),
            optimization: create_test_optimization_config(),
        };

        let context = Context::new(&system, &config, &ff, &library);
        let design_residues = context.resolve_design_residues().unwrap();

        assert_eq!(design_residues.len(), 1);
        let chain_a = system.find_chain_by_id('A').unwrap();
        let ala_id = system.find_residue_by_id(chain_a, 1).unwrap();
        assert!(design_residues.contains(&ala_id));
    }

    #[test]
    fn resolve_design_residues_returns_empty_for_placement_config() {
        let system = create_test_system();
        let library = create_test_rotamer_library();
        let ff = create_test_forcefield();

        let config = PlacementConfig {
            residues_to_optimize: ResidueSelection::All,
            scoring: create_test_scoring_config(),
            optimization: create_test_optimization_config(),
        };

        let context = Context::new(&system, &config, &ff, &library);
        let design_residues = context.resolve_design_residues().unwrap();

        assert!(design_residues.is_empty());
    }

    #[test]
    fn resolve_all_active_residues_combines_repack_and_design() {
        let system = create_test_system();
        let library = create_test_rotamer_library();
        let ff = create_test_forcefield();

        let design_specifier = ResidueSpecifier {
            chain_id: 'A',
            residue_number: 1,
        };
        let mut design_spec = DesignSpec::new();
        design_spec.insert(design_specifier, vec![ResidueType::Leucine]);

        let repack_specifier = ResidueSpecifier {
            chain_id: 'A',
            residue_number: 3,
        };
        let repack_selection = ResidueSelection::List {
            include: vec![repack_specifier],
            exclude: vec![],
        };

        let config = DesignConfig {
            design_spec,
            neighbors_to_repack: repack_selection,
            scoring: create_test_scoring_config(),
            optimization: create_test_optimization_config(),
        };

        let context = Context::new(&system, &config, &ff, &library);
        let all_active = context.resolve_all_active_residues().unwrap();

        assert_eq!(all_active.len(), 2);
        let chain_a = system.find_chain_by_id('A').unwrap();
        let ala_id = system.find_residue_by_id(chain_a, 1).unwrap();
        let leu_id = system.find_residue_by_id(chain_a, 3).unwrap();
        assert!(all_active.contains(&ala_id));
        assert!(all_active.contains(&leu_id));
    }
}
