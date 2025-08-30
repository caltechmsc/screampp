use crate::core::models::chain::ChainType;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::core::models::system::MolecularSystem;
use crate::core::rotamers::library::RotamerLibrary;
use crate::engine::config::ResidueSelection;
use crate::engine::error::EngineError;
use kiddo::{KdTree, SquaredEuclidean};
use std::collections::{HashMap, HashSet};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub fn collect_active_sidechain_atoms(
    system: &MolecularSystem,
    active_residues: &HashSet<ResidueId>,
) -> HashMap<ResidueId, Vec<AtomId>> {
    let mut map = HashMap::with_capacity(active_residues.len());
    for &residue_id in active_residues {
        if let Some(residue) = system.residue(residue_id) {
            let sidechain_ids: Vec<AtomId> = residue
                .atoms()
                .iter()
                .filter_map(|&atom_id| {
                    system.atom(atom_id).and_then(|atom| {
                        if atom.role == crate::core::models::atom::AtomRole::Sidechain {
                            Some(atom_id)
                        } else {
                            None
                        }
                    })
                })
                .collect();
            map.insert(residue_id, sidechain_ids);
        }
    }
    map
}

pub fn precompute_environment_atoms(
    system: &MolecularSystem,
    active_residues: &HashSet<ResidueId>,
) -> Vec<AtomId> {
    system
        .atoms_iter()
        .filter_map(|(atom_id, atom)| {
            if !active_residues.contains(&atom.residue_id)
                || atom.role == crate::core::models::atom::AtomRole::Backbone
            {
                Some(atom_id)
            } else {
                None
            }
        })
        .collect()
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
            let iterator = system.residues_iter();

            #[cfg(feature = "parallel")]
            let iterator = system.residues_iter().par_bridge();

            let binding_site_ids: HashSet<ResidueId> = iterator
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
                                let nearest = kdtree.nearest_one::<SquaredEuclidean>(&protein_pos);
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
