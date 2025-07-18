use super::config::{DesignConfig, DesignSpec, PlacementConfig, ResidueSelection};
use super::error::EngineError;
use super::progress::ProgressReporter;
use crate::core::forcefield::params::Forcefield;
use crate::core::models::chain::ChainType;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::core::rotamers::library::RotamerLibrary;
use kiddo::KdTree;
use kiddo::SquaredEuclidean;
use std::collections::HashSet;

#[cfg(feature = "parallel")]
use rayon::prelude::*;
#[cfg(feature = "parallel")]
use std::sync::Arc;

#[derive(Clone, Copy)]
pub struct Context<'a, C> {
    pub system: &'a MolecularSystem,
    pub config: &'a C,
    pub forcefield: &'a Forcefield,
    pub rotamer_library: &'a RotamerLibrary,
    pub reporter: &'a ProgressReporter<'a>,
}

impl<'a, C> Context<'a, C> {
    pub fn new(
        system: &'a MolecularSystem,
        config: &'a C,
        forcefield: &'a Forcefield,
        rotamer_library: &'a RotamerLibrary,
        reporter: &'a ProgressReporter<'a>,
    ) -> Self {
        Self {
            system,
            config,
            forcefield,
            rotamer_library,
            reporter,
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
                        if crate::core::utils::identifiers::is_heavy_atom(&atom.name) {
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
                        if crate::core::utils::identifiers::is_heavy_atom(&protein_atom.name) {
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
    use crate::core::models::atom::Atom;
    use crate::core::models::chain::ChainType;
    use crate::core::models::residue::ResidueType;
    use crate::core::models::system::MolecularSystem;
    use crate::core::rotamers::library::RotamerLibrary;
    use crate::engine::config::ResidueSpecifier;
    use nalgebra::Point3;

    fn create_simple_test_system() -> (MolecularSystem, RotamerLibrary) {
        let mut system = MolecularSystem::new();
        let chain_a_id = system.add_chain('A', ChainType::Protein);
        system
            .add_residue(chain_a_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        system
            .add_residue(chain_a_id, 2, "GLY", Some(ResidueType::Glycine))
            .unwrap();

        let mut library = RotamerLibrary::default();
        library
            .rotamers
            .insert(ResidueType::Alanine, Default::default());
        library
            .rotamers
            .insert(ResidueType::Glycine, Default::default());

        (system, library)
    }

    fn create_test_system_with_ligand_and_protein() -> MolecularSystem {
        let mut system = MolecularSystem::new();

        let chain_a_id = system.add_chain('A', ChainType::Protein);
        let res_ala_id = system
            .add_residue(chain_a_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let res_gly_id = system
            .add_residue(chain_a_id, 2, "GLY", Some(ResidueType::Glycine))
            .unwrap();
        let res_leu_id = system
            .add_residue(chain_a_id, 3, "LEU", Some(ResidueType::Leucine))
            .unwrap();
        let res_phe_id = system
            .add_residue(chain_a_id, 4, "PHE", Some(ResidueType::Phenylalanine))
            .unwrap();

        let mut ala_ca = Atom::new(10, "CA", res_ala_id, Point3::new(1.0, 0.0, 0.0));
        ala_ca.force_field_type = "C_sp3".to_string();
        let mut ala_cb = Atom::new(11, "CB", res_ala_id, Point3::new(0.5, 0.5, 0.0));
        ala_cb.force_field_type = "C_sp3".to_string();
        system.add_atom_to_residue(res_ala_id, ala_ca).unwrap();
        system.add_atom_to_residue(res_ala_id, ala_cb).unwrap();

        let mut gly_ca = Atom::new(20, "CA", res_gly_id, Point3::new(2.5, 0.0, 0.0));
        gly_ca.force_field_type = "C_sp3".to_string();
        system.add_atom_to_residue(res_gly_id, gly_ca).unwrap();

        let mut leu_ca = Atom::new(30, "CA", res_leu_id, Point3::new(10.0, 0.0, 0.0));
        leu_ca.force_field_type = "C_sp3".to_string();
        system.add_atom_to_residue(res_leu_id, leu_ca).unwrap();

        let mut phe_ca = Atom::new(40, "CA", res_phe_id, Point3::new(12.0, 0.0, 0.0));
        phe_ca.force_field_type = "C_sp3".to_string();
        system.add_atom_to_residue(res_phe_id, phe_ca).unwrap();

        let chain_l_id = system.add_chain('L', ChainType::Ligand);
        let res_lig_id = system.add_residue(chain_l_id, 1, "LIG", None).unwrap();

        let mut lig_c1 = Atom::new(1, "C1", res_lig_id, Point3::new(0.0, 0.0, 0.0));
        lig_c1.force_field_type = "C_sp3".to_string();
        let mut lig_h1 = Atom::new(2, "H1", res_lig_id, Point3::new(-0.5, 0.0, 0.0));
        lig_h1.force_field_type = "H_".to_string();
        let mut lig_c2 = Atom::new(3, "C2", res_lig_id, Point3::new(2.0, 0.0, 0.0));
        lig_c2.force_field_type = "C_sp3".to_string();
        let mut lig_n1 = Atom::new(4, "N1", res_lig_id, Point3::new(0.0, 1.0, 0.0));
        lig_n1.force_field_type = "N_ar".to_string();

        system.add_atom_to_residue(res_lig_id, lig_c1).unwrap();
        system.add_atom_to_residue(res_lig_id, lig_h1).unwrap();
        system.add_atom_to_residue(res_lig_id, lig_c2).unwrap();
        system.add_atom_to_residue(res_lig_id, lig_n1).unwrap();

        system
    }

    fn create_test_rotamer_library() -> RotamerLibrary {
        let mut lib = RotamerLibrary::default();
        lib.rotamers.insert(ResidueType::Alanine, vec![]);
        lib
    }

    #[test]
    fn all_selection_selects_all_valid_residues() {
        let (system, library) = create_simple_test_system();
        let selection = ResidueSelection::All;
        let active_residues = resolve_selection_to_ids(&system, &selection, &library).unwrap();
        assert_eq!(active_residues.len(), 2);
    }

    #[test]
    fn list_selection_with_include() {
        let (system, library) = create_simple_test_system();
        let selection = ResidueSelection::List {
            include: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 1,
            }],
            exclude: vec![],
        };
        let active_residues = resolve_selection_to_ids(&system, &selection, &library).unwrap();
        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let ala_id = system.find_residue_by_id(chain_a_id, 1).unwrap();
        assert_eq!(active_residues.len(), 1);
        assert!(active_residues.contains(&ala_id));
    }

    #[test]
    fn list_selection_with_exclude() {
        let (system, library) = create_simple_test_system();
        let selection = ResidueSelection::List {
            include: vec![],
            exclude: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 1,
            }],
        };
        let active_residues = resolve_selection_to_ids(&system, &selection, &library).unwrap();
        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let gly_id = system.find_residue_by_id(chain_a_id, 2).unwrap();
        assert_eq!(active_residues.len(), 1);
        assert!(active_residues.contains(&gly_id));
    }

    #[test]
    fn list_selection_with_include_and_exclude() {
        let (system, library) = create_simple_test_system();
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
            exclude: vec![ResidueSpecifier {
                chain_id: 'A',
                residue_number: 1,
            }],
        };
        let active_residues = resolve_selection_to_ids(&system, &selection, &library).unwrap();
        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let gly_id = system.find_residue_by_id(chain_a_id, 2).unwrap();
        assert_eq!(active_residues.len(), 1);
        assert!(active_residues.contains(&gly_id));
    }

    #[test]
    fn selection_filters_residues_not_in_rotamer_library() {
        let (system, mut library) = create_simple_test_system();
        library.rotamers.remove(&ResidueType::Glycine);

        let selection = ResidueSelection::All;
        let active_residues = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let ala_id = system.find_residue_by_id(chain_a_id, 1).unwrap();

        assert_eq!(active_residues.len(), 1);
        assert!(active_residues.contains(&ala_id));
    }

    #[test]
    fn ligand_binding_site_selects_close_protein_residues() {
        let system = create_test_system_with_ligand_and_protein();
        let library = create_test_rotamer_library();

        let ligand_specifier = ResidueSpecifier {
            chain_id: 'L',
            residue_number: 1,
        };
        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue: ligand_specifier,
            radius_angstroms: 2.0,
        };

        let active_residues = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let ala_id = system.find_residue_by_id(chain_a_id, 1).unwrap();
        let gly_id = system.find_residue_by_id(chain_a_id, 2).unwrap();
        let leu_id = system.find_residue_by_id(chain_a_id, 3).unwrap();
        let phe_id = system.find_residue_by_id(chain_a_id, 4).unwrap();

        assert!(active_residues.contains(&ala_id));

        assert!(!active_residues.contains(&gly_id));

        assert!(!active_residues.contains(&leu_id));
        assert!(!active_residues.contains(&phe_id));

        assert_eq!(active_residues.len(), 1);
    }

    #[test]
    fn ligand_binding_site_filters_by_rotamer_library() {
        let system = create_test_system_with_ligand_and_protein();
        let mut library_only_ala = RotamerLibrary::default();
        library_only_ala
            .rotamers
            .insert(ResidueType::Alanine, Default::default());

        let ligand_specifier = ResidueSpecifier {
            chain_id: 'L',
            residue_number: 1,
        };
        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue: ligand_specifier,
            radius_angstroms: 3.0,
        };

        let active_residues =
            resolve_selection_to_ids(&system, &selection, &library_only_ala).unwrap();

        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let ala_id = system.find_residue_by_id(chain_a_id, 1).unwrap();
        let gly_id = system.find_residue_by_id(chain_a_id, 2).unwrap();

        assert!(active_residues.contains(&ala_id));
        assert!(!active_residues.contains(&gly_id));
        assert_eq!(active_residues.len(), 1);
    }

    #[test]
    fn ligand_binding_site_handles_nonexistent_ligand() {
        let system = create_test_system_with_ligand_and_protein();
        let library = create_test_rotamer_library();

        let nonexistent_ligand = ResidueSpecifier {
            chain_id: 'X',
            residue_number: 99,
        };
        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue: nonexistent_ligand.clone(),
            radius_angstroms: 5.0,
        };

        let result = resolve_selection_to_ids(&system, &selection, &library);

        assert!(
            matches!(result, Err(EngineError::ResidueNotFound { spec }) if spec == nonexistent_ligand)
        );
    }

    #[test]
    fn ligand_binding_site_handles_ligand_with_no_heavy_atoms() {
        let mut system = MolecularSystem::new();
        let library = create_test_rotamer_library();

        let chain_a_id = system.add_chain('A', ChainType::Protein);
        let res_ala_id = system
            .add_residue(chain_a_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let mut ala_ca = Atom::new(10, "CA", res_ala_id, Point3::new(1.0, 0.0, 0.0));
        ala_ca.force_field_type = "C_sp3".to_string();
        system.add_atom_to_residue(res_ala_id, ala_ca).unwrap();

        let chain_l_id = system.add_chain('L', ChainType::Ligand);
        let res_lig_id = system.add_residue(chain_l_id, 1, "LIG", None).unwrap();
        let mut lig_h1 = Atom::new(1, "H1", res_lig_id, Point3::new(0.0, 0.0, 0.0));
        lig_h1.force_field_type = "H_".to_string();
        system.add_atom_to_residue(res_lig_id, lig_h1).unwrap();

        let ligand_specifier = ResidueSpecifier {
            chain_id: 'L',
            residue_number: 1,
        };
        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue: ligand_specifier,
            radius_angstroms: 2.0,
        };

        let active_residues = resolve_selection_to_ids(&system, &selection, &library).unwrap();
        assert!(active_residues.is_empty());
    }

    #[test]
    fn ligand_binding_site_excludes_ligand_itself() {
        let mut system = create_test_system_with_ligand_and_protein();
        let library = create_test_rotamer_library();

        let ligand_specifier = ResidueSpecifier {
            chain_id: 'L',
            residue_number: 1,
        };
        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue: ligand_specifier.clone(),
            radius_angstroms: 5.0,
        };

        let active_residues = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        let chain_l_id = system.find_chain_by_id('L').unwrap();
        let lig_res_id = system.find_residue_by_id(chain_l_id, 1).unwrap();

        assert!(!active_residues.contains(&lig_res_id));
    }

    #[test]
    fn ligand_binding_site_works_with_large_radius() {
        let system = create_test_system_with_ligand_and_protein();
        let mut library = create_test_rotamer_library();
        library
            .rotamers
            .insert(ResidueType::Leucine, Default::default());
        library
            .rotamers
            .insert(ResidueType::Phenylalanine, Default::default());

        let ligand_specifier = ResidueSpecifier {
            chain_id: 'L',
            residue_number: 1,
        };
        let selection = ResidueSelection::LigandBindingSite {
            ligand_residue: ligand_specifier,
            radius_angstroms: 20.0,
        };

        let active_residues = resolve_selection_to_ids(&system, &selection, &library).unwrap();

        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let ala_id = system.find_residue_by_id(chain_a_id, 1).unwrap();
        let gly_id = system.find_residue_by_id(chain_a_id, 2).unwrap();
        let leu_id = system.find_residue_by_id(chain_a_id, 3).unwrap();
        let phe_id = system.find_residue_by_id(chain_a_id, 4).unwrap();

        assert!(active_residues.contains(&ala_id));
        assert!(!active_residues.contains(&gly_id));
        assert!(active_residues.contains(&leu_id));
        assert!(active_residues.contains(&phe_id));
        assert_eq!(active_residues.len(), 3);
    }
}
