use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::ids::ResidueId;
use crate::core::models::residue::ResidueType;
use std::collections::HashMap;

#[derive(Debug, Default, Clone)]
pub struct ELCache {
    data: HashMap<(ResidueId, ResidueType), HashMap<usize, EnergyTerm>>,
}

impl ELCache {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn insert(
        &mut self,
        residue_id: ResidueId,
        residue_type: ResidueType,
        rotamer_idx: usize,
        energy: EnergyTerm,
    ) {
        self.data
            .entry((residue_id, residue_type))
            .or_default()
            .insert(rotamer_idx, energy);
    }

    pub fn get(
        &self,
        residue_id: ResidueId,
        residue_type: ResidueType,
        rotamer_idx: usize,
    ) -> Option<&EnergyTerm> {
        self.data
            .get(&(residue_id, residue_type))
            .and_then(|inner_map| inner_map.get(&rotamer_idx))
    }

    pub fn get_energies_for(
        &self,
        residue_id: ResidueId,
        residue_type: ResidueType,
    ) -> Option<&HashMap<usize, EnergyTerm>> {
        self.data.get(&(residue_id, residue_type))
    }

    pub fn find_ground_state_for(
        &self,
        residue_id: ResidueId,
        residue_type: ResidueType,
    ) -> Option<(usize, &EnergyTerm)> {
        self.get_energies_for(residue_id, residue_type)
            .and_then(|energies| {
                energies
                    .iter()
                    .min_by(|(_, term_a), (_, term_b)| {
                        term_a
                            .total()
                            .partial_cmp(&term_b.total())
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .map(|(idx, term)| (*idx, term))
            })
    }
}
