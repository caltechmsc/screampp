use super::config::ResidueSelection;
use super::error::EngineError;
use crate::core::forcefield::params::Forcefield;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::core::rotamers::library::RotamerLibrary;
use std::collections::HashSet;

#[derive(Clone)]
pub struct Context<'a, C> {
    pub system: &'a MolecularSystem,
    pub forcefield: &'a Forcefield,
    pub config: &'a C,
}

impl<'a, C> Context<'a, C> {
    pub fn new(system: &'a MolecularSystem, forcefield: &'a Forcefield, config: &'a C) -> Self {
        Self {
            system,
            forcefield,
            config,
        }
    }
}

#[derive(Clone)]
pub struct OptimizationContext<'a, C> {
    pub universal: Context<'a, C>,
    pub rotamer_library: &'a RotamerLibrary,
    pub active_residues: &'a HashSet<ResidueId>,
}

impl<'a, C> OptimizationContext<'a, C> {
    pub fn new(
        universal: Context<'a, C>,
        rotamer_library: &'a RotamerLibrary,
        active_residues: &'a HashSet<ResidueId>,
    ) -> Self {
        Self {
            universal,
            rotamer_library,
            active_residues,
        }
    }

    pub fn with_rotamer_library(&self, new_library: &'a RotamerLibrary) -> Self {
        Self {
            universal: Context {
                system: self.universal.system,
                forcefield: self.universal.forcefield,
                config: self.universal.config,
            },
            rotamer_library: new_library,
            active_residues: self.active_residues,
        }
    }
}

pub fn resolve_residue_selection(
    system: &MolecularSystem,
    selection: &ResidueSelection,
) -> Result<HashSet<ResidueId>, EngineError> {
    match selection {
        ResidueSelection::All => {
            let placeable_residues = system
                .residues_iter()
                .filter_map(|(id, res)| {
                    if res.res_type != Some(crate::core::models::residue::ResidueType::Glycine) {
                        Some(id)
                    } else {
                        None
                    }
                })
                .collect();
            Ok(placeable_residues)
        }
        ResidueSelection::List { include, exclude } => {
            let mut final_ids = if include.is_empty() && !exclude.is_empty() {
                resolve_residue_selection(system, &ResidueSelection::All)?
            } else {
                let mut included_ids = HashSet::with_capacity(include.len());
                for spec in include {
                    let chain_id = system
                        .find_chain_by_id(spec.chain_id)
                        .ok_or_else(|| EngineError::ResidueNotFound { spec: spec.clone() })?;
                    let residue_id = system
                        .find_residue_by_id(chain_id, spec.residue_number)
                        .ok_or_else(|| EngineError::ResidueNotFound { spec: spec.clone() })?;
                    included_ids.insert(residue_id);
                }
                included_ids
            };

            for spec in exclude {
                if let Some(chain_id) = system.find_chain_by_id(spec.chain_id) {
                    if let Some(residue_id) =
                        system.find_residue_by_id(chain_id, spec.residue_number)
                    {
                        final_ids.remove(&residue_id);
                    }
                }
            }
            Ok(final_ids)
        }
        ResidueSelection::LigandBindingSite { .. } => {
            // TODO: Implement resolution of ligand binding site residues
            unimplemented!("Ligand binding site selection is not yet implemented.");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::params::Forcefield;
    use crate::core::models::chain::ChainType;
    use crate::core::models::residue::ResidueType;
    use crate::engine::config::ResidueSpecifier;

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

    #[test]
    fn context_creation_is_correct() {
        let system = MolecularSystem::new();
        let ff = create_test_forcefield();
        let config = ();
        let context = Context::new(&system, &ff, &config);
        assert!(std::ptr::eq(context.system, &system));
        assert!(std::ptr::eq(context.forcefield, &ff));
        assert!(std::ptr::eq(context.config, &config));
    }

    #[test]
    fn optimization_context_creation_is_correct() {
        let system = MolecularSystem::new();
        let ff = create_test_forcefield();
        let config = ();
        let rot_lib = RotamerLibrary::default();
        let active_res = HashSet::new();
        let context = Context::new(&system, &ff, &config);
        let opt_context = OptimizationContext::new(context.clone(), &rot_lib, &active_res);

        assert!(std::ptr::eq(opt_context.universal.system, &system));
        assert!(std::ptr::eq(opt_context.rotamer_library, &rot_lib));
        assert!(std::ptr::eq(opt_context.active_residues, &active_res));
    }

    #[test]
    fn optimization_context_with_rotamer_library_updates_library() {
        let system = MolecularSystem::new();
        let ff = create_test_forcefield();
        let config = ();
        let rot_lib1 = RotamerLibrary::default();
        let rot_lib2 = RotamerLibrary::default();
        let active_res = HashSet::new();
        let context = Context::new(&system, &ff, &config);
        let opt_context1 = OptimizationContext::new(context.clone(), &rot_lib1, &active_res);

        let opt_context2 = opt_context1.with_rotamer_library(&rot_lib2);

        assert!(std::ptr::eq(opt_context1.rotamer_library, &rot_lib1));
        assert!(std::ptr::eq(opt_context2.rotamer_library, &rot_lib2));
        assert!(std::ptr::eq(
            opt_context1.universal.system,
            opt_context2.universal.system
        ));
    }

    #[test]
    fn resolve_all_residues_selects_all_except_glycine() {
        let system = create_test_system();
        let selection = ResidueSelection::All;
        let resolved = resolve_residue_selection(&system, &selection).unwrap();

        assert_eq!(resolved.len(), 3);

        let chain_a = system.find_chain_by_id('A').unwrap();
        let ala_id = system.find_residue_by_id(chain_a, 1).unwrap();
        let gly_id = system.find_residue_by_id(chain_a, 2).unwrap();
        let leu_id = system.find_residue_by_id(chain_a, 3).unwrap();
        let chain_b = system.find_chain_by_id('B').unwrap();
        let trp_id = system.find_residue_by_id(chain_b, 10).unwrap();

        assert!(resolved.contains(&ala_id));
        assert!(resolved.contains(&leu_id));
        assert!(resolved.contains(&trp_id));
        assert!(!resolved.contains(&gly_id));
    }

    #[test]
    fn resolve_list_with_include_only_selects_specified_residues() {
        let system = create_test_system();
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
        let resolved = resolve_residue_selection(&system, &selection).unwrap();

        assert_eq!(resolved.len(), 2);
        let chain_a = system.find_chain_by_id('A').unwrap();
        let ala_id = system.find_residue_by_id(chain_a, 1).unwrap();
        let chain_b = system.find_chain_by_id('B').unwrap();
        let trp_id = system.find_residue_by_id(chain_b, 10).unwrap();
        assert!(resolved.contains(&ala_id));
        assert!(resolved.contains(&trp_id));
    }

    #[test]
    fn resolve_list_with_exclude_only_selects_all_except_excluded() {
        let system = create_test_system();
        let selection = ResidueSelection::List {
            include: vec![],
            exclude: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 1,
            }],
        };
        let resolved = resolve_residue_selection(&system, &selection).unwrap();

        assert_eq!(resolved.len(), 2);
        let chain_a = system.find_chain_by_id('A').unwrap();
        let ala_id = system.find_residue_by_id(chain_a, 1).unwrap();
        let leu_id = system.find_residue_by_id(chain_a, 3).unwrap();
        let chain_b = system.find_chain_by_id('B').unwrap();
        let trp_id = system.find_residue_by_id(chain_b, 10).unwrap();

        assert!(!resolved.contains(&ala_id));
        assert!(resolved.contains(&leu_id));
        assert!(resolved.contains(&trp_id));
    }

    #[test]
    fn resolve_list_with_include_and_exclude_works_correctly() {
        let system = create_test_system();
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
        let resolved = resolve_residue_selection(&system, &selection).unwrap();

        assert_eq!(resolved.len(), 1);
        let chain_a = system.find_chain_by_id('A').unwrap();
        let leu_id = system.find_residue_by_id(chain_a, 3).unwrap();
        assert!(resolved.contains(&leu_id));
    }

    #[test]
    fn resolve_list_with_nonexistent_include_residue_returns_error() {
        let system = create_test_system();
        let spec = ResidueSpecifier {
            chain_id: 'Z',
            residue_number: 99,
        };
        let selection = ResidueSelection::List {
            include: vec![spec.clone()],
            exclude: vec![],
        };
        let result = resolve_residue_selection(&system, &selection);

        assert!(matches!(
            result,
            Err(EngineError::ResidueNotFound { spec: e_spec }) if e_spec == spec
        ));
    }

    #[test]
    fn resolve_list_with_nonexistent_exclude_residue_is_ignored() {
        let system = create_test_system();
        let selection = ResidueSelection::List {
            include: vec![],
            exclude: vec![ResidueSpecifier {
                chain_id: 'Z',
                residue_number: 99,
            }],
        };
        let resolved = resolve_residue_selection(&system, &selection).unwrap();
        assert_eq!(resolved.len(), 3);
    }

    #[test]
    fn resolve_list_with_empty_include_and_exclude_returns_empty_set() {
        let system = create_test_system();
        let selection = ResidueSelection::List {
            include: vec![],
            exclude: vec![],
        };
        let resolved = resolve_residue_selection(&system, &selection).unwrap();
        assert!(resolved.is_empty());
    }

    // TODO: Implement ligand binding site resolution and update this test
    #[test]
    #[should_panic(expected = "Ligand binding site selection is not yet implemented.")]
    fn resolve_ligand_binding_site_panics() {
        let system = create_test_system();
        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue_name: "LIG".to_string(),
            radius_angstroms: 5.0,
        };
        let _ = resolve_residue_selection(&system, &selection);
    }
}
