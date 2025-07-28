use super::config::{AnalyzeConfig, DesignConfig, DesignSpec, PlacementConfig, ResidueSelection};
use super::error::EngineError;
use super::progress::ProgressReporter;
use crate::core::forcefield::params::Forcefield;
use crate::core::models::chain::ChainType;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::core::rotamers::library::RotamerLibrary;
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
    ) -> Self {
        Self {
            system,
            forcefield,
            reporter,
            config,
            rotamer_library,
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
}

impl<'a> AnalysisContext<'a> {
    pub fn new(
        system: &'a MolecularSystem,
        forcefield: &'a Forcefield,
        reporter: &'a ProgressReporter<'a>,
        config: &'a AnalyzeConfig,
    ) -> Self {
        Self {
            system,
            forcefield,
            reporter,
            config,
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
            let residue_iter = system
                .residues_iter()
                .filter_map(|(res_id, res)| Some((res_id, res)));

            #[cfg(feature = "parallel")]
            let residue_iter = system
                .residues_iter()
                .par_bridge()
                .filter_map(|(res_id, res)| Some((res_id, res)));

            for (res_id, residue) in residue_iter {
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

                let mut is_in_binding_site = false;
                for protein_atom_id in residue.atoms() {
                    if let Some(protein_atom) = system.atom(*protein_atom_id) {
                        if is_heavy_atom(&protein_atom.name) {
                            let protein_pos = [
                                protein_atom.position.x,
                                protein_atom.position.y,
                                protein_atom.position.z,
                            ];

                            let nearest = kdtree.nearest_one::<SquaredEuclidean>(&protein_pos);

                            if nearest.distance <= radius_sq {
                                is_in_binding_site = true;
                                break;
                            }
                        }
                    }
                }

                if is_in_binding_site {
                    candidate_ids.insert(res_id);
                }
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

pub fn is_heavy_atom(atom_name: &str) -> bool {
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
    use crate::core::models::chain::ChainType;
    use crate::core::models::residue::ResidueType;
    use crate::core::rotamers::library::RotamerLibrary;
    use crate::core::rotamers::rotamer::Rotamer;
    use crate::engine::config::{ResidueSelection, ResidueSpecifier};
    use nalgebra::Point3;
    use std::collections::HashMap;

    fn create_test_system_and_library()
    -> (MolecularSystem, RotamerLibrary, HashMap<String, ResidueId>) {
        let mut system = MolecularSystem::new();
        let chain_a_id = system.add_chain('A', ChainType::Protein);
        let res_gly_id = system
            .add_residue(chain_a_id, 1, "GLY", Some(ResidueType::Glycine))
            .unwrap();
        let res_ala_id = system
            .add_residue(chain_a_id, 2, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let res_lys_id = system
            .add_residue(chain_a_id, 3, "LYS", Some(ResidueType::Lysine))
            .unwrap();

        let atom_ala =
            crate::core::models::atom::Atom::new("CA", res_ala_id, Point3::new(5.0, 0.0, 0.0));
        system.add_atom_to_residue(res_ala_id, atom_ala).unwrap();

        let atom_lys =
            crate::core::models::atom::Atom::new("CA", res_lys_id, Point3::new(10.0, 0.0, 0.0));
        system.add_atom_to_residue(res_lys_id, atom_lys).unwrap();

        let chain_b_id = system.add_chain('B', ChainType::Ligand);
        let res_lig_id = system.add_residue(chain_b_id, 10, "LIG", None).unwrap();
        let atom_lig =
            crate::core::models::atom::Atom::new("C1", res_lig_id, Point3::new(0.0, 0.0, 0.0));
        system.add_atom_to_residue(res_lig_id, atom_lig).unwrap();

        let mut library = RotamerLibrary::default();
        library
            .rotamers
            .insert(ResidueType::Alanine, vec![Rotamer { atoms: vec![] }]);
        library
            .rotamers
            .insert(ResidueType::Lysine, vec![Rotamer { atoms: vec![] }]);

        let mut ids = HashMap::new();
        ids.insert("GLY".to_string(), res_gly_id);
        ids.insert("ALA".to_string(), res_ala_id);
        ids.insert("LYS".to_string(), res_lys_id);
        ids.insert("LIG".to_string(), res_lig_id);

        (system, library, ids)
    }

    #[test]
    fn resolve_all_selection_returns_residues_with_rotamers() {
        let (system, library, ids) = create_test_system_and_library();
        let selection = ResidueSelection::All;

        let result = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        assert_eq!(result.len(), 2);
        assert!(result.contains(ids.get("ALA").unwrap()));
        assert!(result.contains(ids.get("LYS").unwrap()));
    }

    #[test]
    fn resolve_list_selection_with_include_only() {
        let (system, library, ids) = create_test_system_and_library();
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

        let result = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        assert_eq!(result.len(), 1);
        assert!(result.contains(ids.get("ALA").unwrap()));
    }

    #[test]
    fn resolve_list_selection_with_include_and_exclude() {
        let (system, library, ids) = create_test_system_and_library();
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

        let result = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        assert_eq!(result.len(), 1);
        assert!(result.contains(ids.get("LYS").unwrap()));
    }

    #[test]
    fn resolve_list_selection_with_exclude_only() {
        let (system, library, ids) = create_test_system_and_library();
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

        let result = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        assert_eq!(result.len(), 1);
        assert!(result.contains(ids.get("LYS").unwrap()));
    }

    #[test]
    fn resolve_list_selection_fails_for_nonexistent_residue() {
        let (system, library, _) = create_test_system_and_library();
        let selection = ResidueSelection::List {
            include: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 99,
            }],
            exclude: vec![],
        };

        let result = resolve_selection_to_ids(&system, &selection, &library);
        assert!(matches!(result, Err(EngineError::ResidueNotFound { .. })));
    }

    #[test]
    fn resolve_ligand_binding_site_selection() {
        let (system, library, ids) = create_test_system_and_library();
        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue: ResidueSpecifier {
                chain_id: 'B',
                residue_number: 10,
            },
            radius_angstroms: 6.0,
        };

        let result = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        assert_eq!(result.len(), 1);
        assert!(result.contains(ids.get("ALA").unwrap()));
    }

    #[test]
    fn resolve_ligand_binding_site_selection_with_larger_radius() {
        let (system, library, ids) = create_test_system_and_library();
        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue: ResidueSpecifier {
                chain_id: 'B',
                residue_number: 10,
            },
            radius_angstroms: 11.0,
        };

        let result = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        assert_eq!(result.len(), 2);
        assert!(result.contains(ids.get("ALA").unwrap()));
        assert!(result.contains(ids.get("LYS").unwrap()));
    }

    #[test]
    fn resolve_ligand_binding_site_selection_empty_when_no_protein_nearby() {
        let (mut system, library, _) = create_test_system_and_library();
        let lig_res_id = system
            .find_residue_by_id(system.find_chain_by_id('B').unwrap(), 10)
            .unwrap();
        let lig_atom_id = system
            .residue(lig_res_id)
            .unwrap()
            .get_atom_id_by_name("C1")
            .unwrap();
        system.atom_mut(lig_atom_id).unwrap().position = Point3::new(100.0, 100.0, 100.0);

        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue: ResidueSpecifier {
                chain_id: 'B',
                residue_number: 10,
            },
            radius_angstroms: 10.0,
        };

        let result = resolve_selection_to_ids(&system, &selection, &library).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn resolve_selection_filters_residues_without_rotamers() {
        let (system, mut library, ids) = create_test_system_and_library();
        library.rotamers.remove(&ResidueType::Alanine);

        let selection = ResidueSelection::All;
        let result = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        assert_eq!(result.len(), 1);
        assert!(result.contains(ids.get("LYS").unwrap()));
        assert!(!result.contains(ids.get("ALA").unwrap()));
    }
}
