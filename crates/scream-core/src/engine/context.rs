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
