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
