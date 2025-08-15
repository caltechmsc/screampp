use super::config::{AnalyzeConfig, DesignConfig, DesignSpec, PlacementConfig, ResidueSelection};
use super::error::EngineError;
use super::progress::ProgressReporter;
use crate::core::{
    forcefield::params::Forcefield,
    models::{chain::ChainType, ids::ResidueId, system::MolecularSystem},
    rotamers::library::RotamerLibrary,
    topology::registry::TopologyRegistry,
};
use kiddo::{KdTree, SquaredEuclidean};
use std::collections::HashSet;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

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

pub fn resolve_selection_to_ids(
    system: &MolecularSystem,
    selection: &ResidueSelection,
    library: &RotamerLibrary,
) -> Result<HashSet<ResidueId>, EngineError> {
    let mut candidate_ids: HashSet<ResidueId> = HashSet::new();

    match selection {
        ResidueSelection::All => {
            candidate_ids = system.residues_iter().map(|(id, _)| id).collect();
        }
        ResidueSelection::List { include, exclude } => {
            if include.is_empty() && !exclude.is_empty() {
                candidate_ids = system.residues_iter().map(|(id, _)| id).collect();
            } else {
                for spec in include {
                    let chain_id = system
                        .find_chain_by_id(spec.chain_id)
                        .ok_or_else(|| EngineError::ResidueNotFound { spec: spec.clone() })?;
                    let residue_id = system
                        .find_residue_by_id(chain_id, spec.residue_number)
                        .ok_or_else(|| EngineError::ResidueNotFound { spec: spec.clone() })?;
                    candidate_ids.insert(residue_id);
                }
            }

            for spec in exclude {
                if let Some(chain_id) = system.find_chain_by_id(spec.chain_id) {
                    if let Some(residue_id) =
                        system.find_residue_by_id(chain_id, spec.residue_number)
                    {
                        candidate_ids.remove(&residue_id);
                    }
                }
            }
        }
        ResidueSelection::LigandBindingSite {
            ligand_residue,
            radius_angstroms,
        } => {
            let ligand_chain_id = system
                .find_chain_by_id(ligand_residue.chain_id)
                .ok_or_else(|| EngineError::ResidueNotFound {
                    spec: ligand_residue.clone(),
                })?;
            let ligand_res_id = system
                .find_residue_by_id(ligand_chain_id, ligand_residue.residue_number)
                .ok_or_else(|| EngineError::ResidueNotFound {
                    spec: ligand_residue.clone(),
                })?;

            let mut ligand_heavy_atom_positions: Vec<[f64; 3]> = Vec::new();
            if let Some(ligand_res) = system.residue(ligand_res_id) {
                for atom_id in ligand_res.atoms() {
                    if let Some(atom) = system.atom(*atom_id) {
                        if is_heavy_atom(&atom.name) {
                            ligand_heavy_atom_positions.push([
                                atom.position.x,
                                atom.position.y,
                                atom.position.z,
                            ]);
                        }
                    }
                }
            }

            if ligand_heavy_atom_positions.is_empty() {
                return Ok(HashSet::new());
            }

            let kdtree: KdTree<f64, 3> = (&ligand_heavy_atom_positions).into();
            let radius_sq = radius_angstroms * radius_angstroms;

            #[cfg(not(feature = "parallel"))]
            {
                for (res_id, residue) in system.residues_iter() {
                    if res_id == ligand_res_id {
                        continue;
                    }
                    if let Some(chain) = system.chain(residue.chain_id) {
                        if chain.chain_type != ChainType::Protein {
                            continue;
                        }
                    } else {
                        continue;
                    }

                    let is_in_binding_site = residue.atoms().iter().any(|protein_atom_id| {
                        if let Some(protein_atom) = system.atom(*protein_atom_id) {
                            if is_heavy_atom(&protein_atom.name) {
                                let protein_pos = [
                                    protein_atom.position.x,
                                    protein_atom.position.y,
                                    protein_atom.position.z,
                                ];
                                let nearest = kdtree.nearest_one::<SquaredEuclidean>(&protein_pos);
                                return nearest.distance <= radius_sq;
                            }
                        }
                        false
                    });

                    if is_in_binding_site {
                        candidate_ids.insert(res_id);
                    }
                }
            }

            #[cfg(feature = "parallel")]
            {
                let binding_site_ids: HashSet<ResidueId> = system
                    .residues_iter()
                    .par_bridge()
                    .filter_map(|(res_id, residue)| {
                        if res_id == ligand_res_id {
                            return None;
                        }
                        if let Some(chain) = system.chain(residue.chain_id) {
                            if chain.chain_type != ChainType::Protein {
                                return None;
                            }
                        } else {
                            return None;
                        }

                        let is_in_binding_site = residue.atoms().iter().any(|protein_atom_id| {
                            if let Some(protein_atom) = system.atom(*protein_atom_id) {
                                if is_heavy_atom(&protein_atom.name) {
                                    let protein_pos = [
                                        protein_atom.position.x,
                                        protein_atom.position.y,
                                        protein_atom.position.z,
                                    ];
                                    let nearest =
                                        kdtree.nearest_one::<SquaredEuclidean>(&protein_pos);
                                    return nearest.distance <= radius_sq;
                                }
                            }
                            false
                        });

                        if is_in_binding_site {
                            Some(res_id)
                        } else {
                            None
                        }
                    })
                    .collect();

                candidate_ids.extend(binding_site_ids);
            }
        }
    };

    let final_active_residues = candidate_ids
        .into_iter()
        .filter(|&residue_id| {
            system
                .residue(residue_id)
                .and_then(|res| res.residue_type)
                .map_or(false, |residue_type| {
                    library.get_rotamers_for(residue_type).is_some()
                })
        })
        .collect();

    Ok(final_active_residues)
}

fn is_heavy_atom(atom_name: &str) -> bool {
    let first_char = atom_name
        .trim()
        .chars()
        .next()
        .map(|c| c.to_ascii_uppercase());
    !matches!(first_char, Some('H') | Some('D'))
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

    mod resolve_selection_to_ids_tests {
        use super::*;

        #[test]
        fn resolve_all_selection_returns_residues_with_rotamers() {
            let setup = setup();
            let selection = ResidueSelection::All;

            let result =
                resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library)
                    .unwrap();

            assert_eq!(result.len(), 2);
            assert!(result.contains(setup.residue_ids.get("ALA").unwrap()));
            assert!(result.contains(setup.residue_ids.get("LYS").unwrap()));
        }

        #[test]
        fn resolve_list_selection_with_include_only() {
            let setup = setup();
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

            let result =
                resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library)
                    .unwrap();

            assert_eq!(result.len(), 1);
            assert!(result.contains(setup.residue_ids.get("ALA").unwrap()));
        }

        #[test]
        fn resolve_list_selection_with_include_and_exclude() {
            let setup = setup();
            let selection = ResidueSelection::List {
                include: vec![
                    ResidueSpecifier {
                        chain_id: 'A',
                        residue_number: 2,
                    },
                    ResidueSpecifier {
                        chain_id: 'A',
                        residue_number: 3,
                    },
                ],
                exclude: vec![ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: 2,
                }],
            };

            let result =
                resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library)
                    .unwrap();

            assert_eq!(result.len(), 1);
            assert!(result.contains(setup.residue_ids.get("LYS").unwrap()));
        }

        #[test]
        fn resolve_list_selection_with_exclude_only() {
            let setup = setup();
            let selection = ResidueSelection::List {
                include: vec![],
                exclude: vec![
                    ResidueSpecifier {
                        chain_id: 'A',
                        residue_number: 1,
                    },
                    ResidueSpecifier {
                        chain_id: 'A',
                        residue_number: 2,
                    },
                ],
            };

            let result =
                resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library)
                    .unwrap();

            assert_eq!(result.len(), 1);
            assert!(result.contains(setup.residue_ids.get("LYS").unwrap()));
        }

        #[test]
        fn resolve_list_selection_fails_for_nonexistent_residue() {
            let setup = setup();
            let selection = ResidueSelection::List {
                include: vec![ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: 99,
                }],
                exclude: vec![],
            };

            let result =
                resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library);
            assert!(matches!(result, Err(EngineError::ResidueNotFound { .. })));
        }

        #[test]
        fn resolve_ligand_binding_site_selection() {
            let setup = setup();
            let selection = ResidueSelection::LigandBindingSite {
                ligand_residue: ResidueSpecifier {
                    chain_id: 'B',
                    residue_number: 10,
                },
                radius_angstroms: 6.0,
            };

            let result =
                resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library)
                    .unwrap();

            assert_eq!(result.len(), 1);
            assert!(result.contains(setup.residue_ids.get("ALA").unwrap()));
        }

        #[test]
        fn resolve_ligand_binding_site_selection_with_larger_radius() {
            let setup = setup();
            let selection = ResidueSelection::LigandBindingSite {
                ligand_residue: ResidueSpecifier {
                    chain_id: 'B',
                    residue_number: 10,
                },
                radius_angstroms: 11.0,
            };

            let result =
                resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library)
                    .unwrap();

            assert_eq!(result.len(), 2);
            assert!(result.contains(setup.residue_ids.get("ALA").unwrap()));
            assert!(result.contains(setup.residue_ids.get("LYS").unwrap()));
        }

        #[test]
        fn resolve_ligand_binding_site_selection_empty_when_no_protein_nearby() {
            let mut setup = setup();
            let lig_res_id = setup
                .system
                .find_residue_by_id(setup.system.find_chain_by_id('B').unwrap(), 10)
                .unwrap();
            let lig_atom_id = setup
                .system
                .residue(lig_res_id)
                .unwrap()
                .get_first_atom_id_by_name("C1")
                .expect("Ligand residue should have a 'C1' atom");

            setup.system.atom_mut(lig_atom_id).unwrap().position = Point3::new(100.0, 100.0, 100.0);

            let selection = ResidueSelection::LigandBindingSite {
                ligand_residue: ResidueSpecifier {
                    chain_id: 'B',
                    residue_number: 10,
                },
                radius_angstroms: 10.0,
            };

            let result =
                resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library)
                    .unwrap();
            assert!(result.is_empty());
        }

        #[test]
        fn resolve_selection_filters_residues_without_rotamers() {
            let mut setup = setup();
            setup.rotamer_library.rotamers.remove(&ResidueType::Alanine);

            let selection = ResidueSelection::All;
            let result =
                resolve_selection_to_ids(&setup.system, &selection, &setup.rotamer_library)
                    .unwrap();

            assert_eq!(result.len(), 1);
            assert!(result.contains(setup.residue_ids.get("LYS").unwrap()));
            assert!(!result.contains(setup.residue_ids.get("ALA").unwrap()));
        }
    }

    mod optimization_context_tests {
        use super::*;

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
}
