use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::ids::ResidueId;
use crate::core::models::residue::ResidueType;
use std::collections::HashMap;

/// Caches energy terms for residue-rotamer combinations to avoid recomputation.
///
/// This struct provides an efficient caching mechanism for storing and retrieving
/// energy calculations for different rotamer conformations of protein residues.
/// It uses a nested HashMap structure to organize energies by residue identity
/// and type, enabling fast lookups during optimization processes.
#[derive(Debug, Default, Clone)]
pub struct ELCache {
    /// Internal storage mapping (residue_id, residue_type) to rotamer energies.
    data: HashMap<(ResidueId, ResidueType), HashMap<usize, EnergyTerm>>,
}

impl ELCache {
    /// Creates a new empty energy cache.
    ///
    /// # Return
    ///
    /// Returns a new `ELCache` instance with no cached energies.
    pub fn new() -> Self {
        Self::default()
    }

    /// Inserts an energy term for a specific residue-rotamer combination.
    ///
    /// If an energy term already exists for the given residue and rotamer,
    /// it will be overwritten with the new value.
    ///
    /// # Arguments
    ///
    /// * `residue_id` - The unique identifier of the residue.
    /// * `residue_type` - The type of the residue (e.g., Alanine, Glycine).
    /// * `rotamer_idx` - The index of the rotamer conformation.
    /// * `energy` - The energy term to cache for this combination.
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

    /// Retrieves the cached energy term for a specific residue-rotamer combination.
    ///
    /// # Arguments
    ///
    /// * `residue_id` - The unique identifier of the residue.
    /// * `residue_type` - The type of the residue.
    /// * `rotamer_idx` - The index of the rotamer conformation.
    ///
    /// # Return
    ///
    /// Returns `Some(&EnergyTerm)` if the combination exists in the cache,
    /// or `None` if no energy has been cached for this combination.
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

    /// Retrieves all cached energy terms for a specific residue.
    ///
    /// This method returns the complete set of rotamer energies that have
    /// been cached for the given residue, allowing iteration over all
    /// available conformations.
    ///
    /// # Arguments
    ///
    /// * `residue_id` - The unique identifier of the residue.
    /// * `residue_type` - The type of the residue.
    ///
    /// # Return
    ///
    /// Returns `Some(&HashMap<usize, EnergyTerm>)` containing all cached
    /// rotamer energies for the residue, or `None` if no energies are cached.
    pub fn get_energies_for(
        &self,
        residue_id: ResidueId,
        residue_type: ResidueType,
    ) -> Option<&HashMap<usize, EnergyTerm>> {
        self.data.get(&(residue_id, residue_type))
    }

    /// Finds the rotamer with the lowest total energy for a specific residue.
    ///
    /// This method searches through all cached rotamers for the given residue
    /// and returns the one with the minimum total energy value. This is useful
    /// for identifying the most favorable conformation.
    ///
    /// # Arguments
    ///
    /// * `residue_id` - The unique identifier of the residue.
    /// * `residue_type` - The type of the residue.
    ///
    /// # Return
    ///
    /// Returns `Some((usize, &EnergyTerm))` containing the rotamer index and
    /// energy of the ground state, or `None` if no energies are cached for
    /// the residue.
    pub fn find_ground_state_for(
        &self,
        residue_id: ResidueId,
        residue_type: ResidueType,
    ) -> Option<(usize, &EnergyTerm)> {
        self.get_energies_for(residue_id, residue_type)
            .and_then(|energies| {
                // Find the rotamer with the lowest total energy by comparing energy values
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

    /// Returns the number of residues for which energies are cached.
    ///
    /// # Return
    ///
    /// Returns the count of unique (residue_id, residue_type) combinations
    /// that have at least one cached energy term.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Checks whether the cache contains any cached energies.
    ///
    /// # Return
    ///
    /// Returns `true` if no energies are cached, `false` otherwise.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::term::EnergyTerm;
    use crate::core::models::ids::ResidueId;
    use crate::core::models::residue::ResidueType;
    use slotmap::KeyData;

    fn dummy_residue_id(n: u64) -> ResidueId {
        ResidueId::from(KeyData::from_ffi(n))
    }

    #[test]
    fn insert_and_get_returns_correct_energy_term() {
        let mut cache = ELCache::new();
        let rid = dummy_residue_id(1);
        let rtype = ResidueType::Alanine;
        let energy = EnergyTerm::new(1.0, 2.0, 3.0);

        cache.insert(rid, rtype, 0, energy);
        let retrieved = cache.get(rid, rtype, 0).unwrap();
        assert_eq!(*retrieved, energy);
    }

    #[test]
    fn get_returns_none_for_nonexistent_entry() {
        let cache = ELCache::new();
        let rid = dummy_residue_id(2);
        let rtype = ResidueType::Glycine;
        assert!(cache.get(rid, rtype, 0).is_none());
    }

    #[test]
    fn get_energies_for_returns_all_rotamers_for_residue() {
        let mut cache = ELCache::new();
        let rid = dummy_residue_id(3);
        let rtype = ResidueType::Leucine;
        cache.insert(rid, rtype, 0, EnergyTerm::new(1.0, 0.0, 0.0));
        cache.insert(rid, rtype, 1, EnergyTerm::new(2.0, 0.0, 0.0));

        let energies = cache.get_energies_for(rid, rtype).unwrap();
        assert_eq!(energies.len(), 2);
        assert!(energies.contains_key(&0));
        assert!(energies.contains_key(&1));
    }

    #[test]
    fn get_energies_for_returns_none_for_unknown_residue() {
        let cache = ELCache::new();
        let rid = dummy_residue_id(4);
        let rtype = ResidueType::Proline;
        assert!(cache.get_energies_for(rid, rtype).is_none());
    }

    #[test]
    fn find_ground_state_for_returns_rotamer_with_lowest_total_energy() {
        let mut cache = ELCache::new();
        let rid = dummy_residue_id(5);
        let rtype = ResidueType::Valine;
        cache.insert(rid, rtype, 0, EnergyTerm::new(1.0, 2.0, 3.0));
        cache.insert(rid, rtype, 1, EnergyTerm::new(-1.0, 0.0, 0.0));
        cache.insert(rid, rtype, 2, EnergyTerm::new(0.0, 0.0, 0.0));

        let (idx, energy) = cache.find_ground_state_for(rid, rtype).unwrap();
        assert_eq!(idx, 1);
        assert_eq!(*energy, EnergyTerm::new(-1.0, 0.0, 0.0));
    }

    #[test]
    fn find_ground_state_for_returns_none_when_no_rotamers() {
        let cache = ELCache::new();
        let rid = dummy_residue_id(6);
        let rtype = ResidueType::Phenylalanine;
        assert!(cache.find_ground_state_for(rid, rtype).is_none());
    }

    #[test]
    fn insert_overwrites_existing_rotamer_energy() {
        let mut cache = ELCache::new();
        let rid = dummy_residue_id(7);
        let rtype = ResidueType::Serine;
        cache.insert(rid, rtype, 0, EnergyTerm::new(1.0, 1.0, 1.0));
        cache.insert(rid, rtype, 0, EnergyTerm::new(2.0, 2.0, 2.0));
        let retrieved = cache.get(rid, rtype, 0).unwrap();
        assert_eq!(*retrieved, EnergyTerm::new(2.0, 2.0, 2.0));
    }

    #[test]
    fn cache_is_empty_on_creation() {
        let cache = ELCache::new();
        assert!(cache.data.is_empty());
    }

    #[test]
    fn len_returns_zero_for_empty_cache() {
        let cache = ELCache::new();
        assert_eq!(cache.len(), 0);
    }

    #[test]
    fn len_returns_number_of_residue_entries() {
        let mut cache = ELCache::new();
        let rid1 = dummy_residue_id(10);
        let rid2 = dummy_residue_id(20);
        cache.insert(rid1, ResidueType::Alanine, 0, EnergyTerm::default());
        cache.insert(rid2, ResidueType::Glycine, 0, EnergyTerm::default());
        assert_eq!(cache.len(), 2);
    }

    #[test]
    fn is_empty_returns_true_for_empty_cache() {
        let cache = ELCache::new();
        assert!(cache.is_empty());
    }

    #[test]
    fn is_empty_returns_false_when_cache_has_entries() {
        let mut cache = ELCache::new();
        let rid = dummy_residue_id(30);
        cache.insert(rid, ResidueType::Leucine, 0, EnergyTerm::default());
        assert!(!cache.is_empty());
    }
}
