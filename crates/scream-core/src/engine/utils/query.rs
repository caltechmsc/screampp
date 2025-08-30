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
