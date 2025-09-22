use super::atom::{Atom, AtomRole};
use super::chain::{Chain, ChainType};
use super::ids::{AtomId, ChainId, ResidueId};
use super::residue::{Residue, ResidueType};
use super::topology::{Bond, BondOrder};
use slotmap::{SecondaryMap, SlotMap};
use std::collections::{HashMap, HashSet};

const CYSTEINE_SULFUR_GAMMA_ATOM_NAME: &str = "SG";

/// Represents a complete molecular system with atoms, residues, chains, and bonds.
///
/// This struct serves as the central data structure for molecular modeling,
/// providing efficient storage and access to all molecular components.
/// It maintains internal caches and lookup maps for performance optimization.
#[derive(Debug, Clone, Default)]
pub struct MolecularSystem {
    /// Primary storage for atoms using a slot map for efficient ID management.
    atoms: SlotMap<AtomId, Atom>,
    /// Primary storage for residues using a slot map for efficient ID management.
    residues: SlotMap<ResidueId, Residue>,
    /// Primary storage for chains using a slot map for efficient ID management.
    chains: SlotMap<ChainId, Chain>,
    /// List of all bonds in the system.
    bonds: Vec<Bond>,
    /// Lookup map for finding residues by chain ID and residue number.
    residue_id_map: HashMap<(ChainId, isize), ResidueId>,
    /// Lookup map for finding chains by their single-character identifier.
    chain_id_map: HashMap<char, ChainId>,
    /// Cached adjacency list for bond connectivity, indexed by atom ID.
    bond_adjacency: SecondaryMap<AtomId, Vec<AtomId>>,
}

impl MolecularSystem {
    /// Creates a new, empty molecular system.
    ///
    /// This constructor initializes all internal data structures
    /// and is ready for adding chains, residues, and atoms.
    pub fn new() -> Self {
        Self::default()
    }

    /// Retrieves an immutable reference to an atom by its ID.
    ///
    /// # Arguments
    ///
    /// * `id` - The atom ID to look up.
    ///
    /// # Return
    ///
    /// Returns `Some(&Atom)` if the atom exists, otherwise `None`.
    pub fn atom(&self, id: AtomId) -> Option<&Atom> {
        self.atoms.get(id)
    }

    /// Retrieves a mutable reference to an atom by its ID.
    ///
    /// # Arguments
    ///
    /// * `id` - The atom ID to look up.
    ///
    /// # Return
    ///
    /// Returns `Some(&mut Atom)` if the atom exists, otherwise `None`.
    pub fn atom_mut(&mut self, id: AtomId) -> Option<&mut Atom> {
        self.atoms.get_mut(id)
    }

    /// Returns an iterator over all atoms in the system.
    ///
    /// # Return
    ///
    /// An iterator yielding `(AtomId, &Atom)` pairs.
    pub fn atoms_iter(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.atoms.iter()
    }

    /// Returns a mutable iterator over all atoms in the system.
    ///
    /// # Return
    ///
    /// An iterator yielding `(AtomId, &mut Atom)` pairs.
    pub fn atoms_iter_mut(&mut self) -> impl Iterator<Item = (AtomId, &mut Atom)> {
        self.atoms.iter_mut()
    }

    /// Retrieves an immutable reference to a residue by its ID.
    ///
    /// # Arguments
    ///
    /// * `id` - The residue ID to look up.
    ///
    /// # Return
    ///
    /// Returns `Some(&Residue)` if the residue exists, otherwise `None`.
    pub fn residue(&self, id: ResidueId) -> Option<&Residue> {
        self.residues.get(id)
    }

    /// Retrieves a mutable reference to a residue by its ID.
    ///
    /// # Arguments
    ///
    /// * `id` - The residue ID to look up.
    ///
    /// # Return
    ///
    /// Returns `Some(&mut Residue)` if the residue exists, otherwise `None`.
    pub fn residue_mut(&mut self, id: ResidueId) -> Option<&mut Residue> {
        self.residues.get_mut(id)
    }

    /// Returns an iterator over all residues in the system.
    ///
    /// # Return
    ///
    /// An iterator yielding `(ResidueId, &Residue)` pairs.
    pub fn residues_iter(&self) -> impl Iterator<Item = (ResidueId, &Residue)> {
        self.residues.iter()
    }

    /// Retrieves an immutable reference to a chain by its ID.
    ///
    /// # Arguments
    ///
    /// * `id` - The chain ID to look up.
    ///
    /// # Return
    ///
    /// Returns `Some(&Chain)` if the chain exists, otherwise `None`.
    pub fn chain(&self, id: ChainId) -> Option<&Chain> {
        self.chains.get(id)
    }

    /// Retrieves a mutable reference to a chain by its ID.
    ///
    /// # Arguments
    ///
    /// * `id` - The chain ID to look up.
    ///
    /// # Return
    ///
    /// Returns `Some(&mut Chain)` if the chain exists, otherwise `None`.
    pub fn chain_mut(&mut self, id: ChainId) -> Option<&mut Chain> {
        self.chains.get_mut(id)
    }

    /// Returns an iterator over all chains in the system.
    ///
    /// # Return
    ///
    /// An iterator yielding `(ChainId, &Chain)` pairs.
    pub fn chains_iter(&self) -> impl Iterator<Item = (ChainId, &Chain)> {
        self.chains.iter()
    }

    /// Returns a slice of all bonds in the system.
    ///
    /// # Return
    ///
    /// A slice containing all bonds.
    pub fn bonds(&self) -> &[Bond] {
        &self.bonds
    }

    /// Finds a chain ID by its single-character identifier.
    ///
    /// # Arguments
    ///
    /// * `id` - The character identifier of the chain.
    ///
    /// # Return
    ///
    /// Returns `Some(ChainId)` if the chain exists, otherwise `None`.
    pub fn find_chain_by_id(&self, id: char) -> Option<ChainId> {
        self.chain_id_map.get(&id).copied()
    }

    /// Finds a residue ID by its chain ID and residue number.
    ///
    /// # Arguments
    ///
    /// * `chain_id` - The ID of the chain containing the residue.
    /// * `residue_number` - The sequential number of the residue.
    ///
    /// # Return
    ///
    /// Returns `Some(ResidueId)` if the residue exists, otherwise `None`.
    pub fn find_residue_by_id(
        &self,
        chain_id: ChainId,
        residue_number: isize,
    ) -> Option<ResidueId> {
        self.residue_id_map
            .get(&(chain_id, residue_number))
            .copied()
    }

    /// Adds a new chain to the system or returns the existing one.
    ///
    /// This method is idempotent; if a chain with the given ID already exists,
    /// it returns the existing chain ID without creating a duplicate.
    ///
    /// # Arguments
    ///
    /// * `id` - The single-character identifier for the chain.
    /// * `chain_type` - The type of the chain.
    ///
    /// # Return
    ///
    /// The ID of the chain (new or existing).
    pub fn add_chain(&mut self, id: char, chain_type: ChainType) -> ChainId {
        *self.chain_id_map.entry(id).or_insert_with(|| {
            let chain = Chain::new(id, chain_type);
            self.chains.insert(chain)
        })
    }

    /// Adds a new residue to the system or returns the existing one.
    ///
    /// This method is idempotent; if a residue with the given chain ID and
    /// residue number already exists, it returns the existing residue ID.
    ///
    /// # Arguments
    ///
    /// * `chain_id` - The ID of the chain to add the residue to.
    /// * `residue_number` - The sequential number of the residue.
    /// * `name` - The name of the residue.
    /// * `residue_type` - The type of the residue, if known.
    ///
    /// # Return
    ///
    /// Returns `Some(ResidueId)` if successful, otherwise `None` (e.g., if chain doesn't exist).
    pub fn add_residue(
        &mut self,
        chain_id: ChainId,
        residue_number: isize,
        name: &str,
        residue_type: Option<ResidueType>,
    ) -> Option<ResidueId> {
        let chain = self.chains.get_mut(chain_id)?;
        let key = (chain_id, residue_number);

        let residue_id = *self.residue_id_map.entry(key).or_insert_with(|| {
            let residue = Residue::new(residue_number, name, residue_type, chain_id);
            self.residues.insert(residue)
        });

        if !chain.residues.contains(&residue_id) {
            chain.residues.push(residue_id);
        }

        Some(residue_id)
    }

    /// Adds an atom to a specific residue.
    ///
    /// This method inserts the atom into the system and registers it with the given residue.
    /// It also initializes the bond adjacency list for the new atom.
    ///
    /// # Arguments
    ///
    /// * `residue_id` - The ID of the residue to add the atom to.
    /// * `atom` - The atom to add.
    ///
    /// # Return
    ///
    /// Returns `Some(AtomId)` if successful, otherwise `None` (e.g., if residue doesn't exist).
    pub fn add_atom_to_residue(&mut self, residue_id: ResidueId, atom: Atom) -> Option<AtomId> {
        if !self.residues.contains_key(residue_id) {
            return None;
        }

        let name = atom.name.clone();

        let atom_id = self.atoms.insert(atom);
        self.bond_adjacency.insert(atom_id, Vec::new());

        let residue = self.residues.get_mut(residue_id).unwrap();
        residue.add_atom(&name, atom_id);

        Some(atom_id)
    }

    /// Adds a bond between two atoms.
    ///
    /// This method creates a bond between the specified atoms and updates
    /// the adjacency cache. It is idempotent; adding an existing bond
    /// succeeds without creating duplicates.
    ///
    /// # Arguments
    ///
    /// * `atom1_id` - ID of the first atom.
    /// * `atom2_id` - ID of the second atom.
    /// * `order` - The order of the bond.
    ///
    /// # Return
    ///
    /// Returns `Some(())` if successful, otherwise `None` (e.g., if atoms don't exist).
    pub fn add_bond(&mut self, atom1_id: AtomId, atom2_id: AtomId, order: BondOrder) -> Option<()> {
        if !self.atoms.contains_key(atom1_id) || !self.atoms.contains_key(atom2_id) {
            return None;
        }

        if let Some(neighbors) = self.bond_adjacency.get(atom1_id) {
            if neighbors.contains(&atom2_id) {
                // Bond already exists, operation is successful (idempotent)
                return Some(());
            }
        }

        self.bonds.push(Bond::new(atom1_id, atom2_id, order));
        self.bond_adjacency[atom1_id].push(atom2_id);
        self.bond_adjacency[atom2_id].push(atom1_id);
        Some(())
    }

    /// Removes an atom from the system.
    ///
    /// This method removes the atom and all associated data, including
    /// bonds and adjacency information. It updates the parent residue
    /// and cleans up all references.
    ///
    /// # Arguments
    ///
    /// * `atom_id` - The ID of the atom to remove.
    ///
    /// # Return
    ///
    /// Returns `Some(Atom)` if the atom existed and was removed, otherwise `None`.
    pub fn remove_atom(&mut self, atom_id: AtomId) -> Option<Atom> {
        let atom = self.atoms.remove(atom_id)?;

        // 1. Remove from parent residue
        if let Some(residue) = self.residues.get_mut(atom.residue_id) {
            residue.remove_atom(&atom.name, atom_id);
        }

        // 2. Remove all bonds connected to this atom
        let original_bonds = std::mem::take(&mut self.bonds);
        self.bonds = original_bonds
            .into_iter()
            .filter(|bond| !bond.contains(atom_id))
            .collect();

        // 3. Clean up adjacency list
        let neighbors = self.bond_adjacency.remove(atom_id).unwrap_or_default();
        for neighbor_id in neighbors {
            if let Some(adjacency) = self.bond_adjacency.get_mut(neighbor_id) {
                adjacency.retain(|&id| id != atom_id);
            }
        }

        Some(atom)
    }

    /// Removes a residue from the system.
    ///
    /// This method removes the residue and all its atoms, updating
    /// the parent chain and cleaning up all references and bonds.
    ///
    /// # Arguments
    ///
    /// * `residue_id` - The ID of the residue to remove.
    ///
    /// # Return
    ///
    /// Returns `Some(Residue)` if the residue existed and was removed, otherwise `None`.
    pub fn remove_residue(&mut self, residue_id: ResidueId) -> Option<Residue> {
        let residue = self.residues.get(residue_id)?.clone(); // Clone to avoid borrow checker issues

        // 1. Remove all atoms within the residue
        for atom_id in residue.atoms().to_vec() {
            // to_vec to avoid borrow issues
            self.remove_atom(atom_id);
        }

        // 2. Remove from parent chain
        if let Some(chain) = self.chains.get_mut(residue.chain_id) {
            chain.residues.retain(|&id| id != residue_id);
        }

        // 3. Remove from residue maps
        self.residue_id_map
            .remove(&(residue.chain_id, residue.residue_number));

        // 4. Finally, remove the residue itself
        self.residues.remove(residue_id)
    }

    /// Retrieves the bonded neighbors of an atom.
    ///
    /// This method returns the list of atoms directly bonded to the given atom,
    /// using the cached adjacency information.
    ///
    /// # Arguments
    ///
    /// * `atom_id` - The ID of the atom to query.
    ///
    /// # Return
    ///
    /// Returns `Some(&[AtomId])` if the atom exists, otherwise `None`.
    pub fn get_bonded_neighbors(&self, atom_id: AtomId) -> Option<&[AtomId]> {
        self.bond_adjacency.get(atom_id).map(|v| v.as_slice())
    }

    /// Returns an iterator over atoms with a specific role.
    ///
    /// # Arguments
    ///
    /// * `role` - The atom role to filter by.
    ///
    /// # Return
    ///
    /// An iterator yielding `(AtomId, &Atom)` pairs for matching atoms.
    pub fn atoms_by_role(&self, role: AtomRole) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.atoms.iter().filter(move |(_, atom)| atom.role == role)
    }

    /// Returns a vector of atom IDs with a specific role.
    ///
    /// # Arguments
    ///
    /// * `role` - The atom role to filter by.
    ///
    /// # Return
    ///
    /// A vector containing the IDs of matching atoms.
    pub fn atom_ids_by_role(&self, role: AtomRole) -> Vec<AtomId> {
        self.atoms_by_role(role).map(|(id, _)| id).collect()
    }

    /// Returns an iterator over protein atoms (backbone and sidechain).
    ///
    /// # Return
    ///
    /// An iterator yielding `(AtomId, &Atom)` pairs for protein atoms.
    pub fn protein_atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.atoms
            .iter()
            .filter(|(_, atom)| matches!(atom.role, AtomRole::Backbone | AtomRole::Sidechain))
    }

    /// Returns a vector of protein atom IDs.
    ///
    /// # Return
    ///
    /// A vector containing the IDs of protein atoms.
    pub fn protein_atom_ids(&self) -> Vec<AtomId> {
        self.protein_atoms().map(|(id, _)| id).collect()
    }

    /// Returns an iterator over background atoms (ligands, water, other).
    ///
    /// # Return
    ///
    /// An iterator yielding `(AtomId, &Atom)` pairs for background atoms.
    pub fn background_atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.atoms.iter().filter(|(_, atom)| {
            matches!(
                atom.role,
                AtomRole::Ligand | AtomRole::Water | AtomRole::Other
            )
        })
    }

    /// Returns a vector of background atom IDs.
    ///
    /// # Return
    ///
    /// A vector containing the IDs of background atoms.
    pub fn background_atom_ids(&self) -> Vec<AtomId> {
        self.background_atoms().map(|(id, _)| id).collect()
    }

    /// Detects and returns the residue IDs of all Cysteine residues involved in disulfide bonds.
    ///
    /// A disulfide bond is identified by a covalent bond between the Sulfur-Gamma (SG)
    /// atoms of two different Cysteine residues. This method correctly handles both
    /// `CYS` and `CYX` residue names by relying on the `ResidueType`.
    ///
    /// # Returns
    ///
    /// A `HashSet<ResidueId>` containing the IDs of all residues participating in
    /// any disulfide bond within the system.
    pub fn find_disulfide_bonded_residues(&self) -> HashSet<ResidueId> {
        let mut bonded_residue_ids = HashSet::new();

        // 1. Collect all Cysteine residues and their SG atoms
        let cysteine_sg_atoms: HashMap<ResidueId, AtomId> = self
            .residues_iter()
            .filter_map(|(res_id, residue)| {
                if matches!(residue.residue_type, Some(ResidueType::Cysteine)) {
                    residue
                        .get_first_atom_id_by_name(CYSTEINE_SULFUR_GAMMA_ATOM_NAME)
                        .map(|sg_id| (res_id, sg_id))
                } else {
                    None
                }
            })
            .collect();

        // If there are fewer than two Cysteines, no disulfide bonds are possible
        if cysteine_sg_atoms.len() < 2 {
            return bonded_residue_ids;
        }

        // Create a reverse map from SG AtomId to ResidueId for quick lookups
        let sg_atom_to_residue: HashMap<AtomId, ResidueId> = cysteine_sg_atoms
            .iter()
            .map(|(&res_id, &atom_id)| (atom_id, res_id))
            .collect();

        // 2. Iterate through the collected SG atoms and check their bonds
        for (res_id_a, sg_atom_id_a) in &cysteine_sg_atoms {
            if bonded_residue_ids.contains(res_id_a) {
                continue;
            }

            if let Some(neighbors) = self.get_bonded_neighbors(*sg_atom_id_a) {
                for &neighbor_atom_id in neighbors {
                    // 3. Check if the neighbor is also an SG atom of another Cysteine
                    if let Some(res_id_b) = sg_atom_to_residue.get(&neighbor_atom_id) {
                        if res_id_a != res_id_b {
                            // Found a disulfide bond!
                            bonded_residue_ids.insert(*res_id_a);
                            bonded_residue_ids.insert(*res_id_b);
                            break;
                        }
                    }
                }
            }
        }

        bonded_residue_ids
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::residue::ResidueType;
    use nalgebra::Point3;

    mod core_functionality {
        use super::*;
        use crate::core::models::residue::ResidueType;
        use nalgebra::Point3;

        struct TestRefs {
            chain_a_id: ChainId,
            gly_id: ResidueId,
            gly_n_id: AtomId,
            gly_ca_id: AtomId,
            ala_id: ResidueId,
            ala_ca_id: AtomId,
        }

        fn create_standard_test_system() -> (MolecularSystem, TestRefs) {
            let mut system = MolecularSystem::new();

            let chain_a_id = system.add_chain('A', ChainType::Protein);

            let gly_id = system
                .add_residue(chain_a_id, 1, "GLY", Some(ResidueType::Glycine))
                .unwrap();
            let gly_n_atom = Atom::new("N", gly_id, Point3::new(0.0, 0.0, 0.0));
            let gly_ca_atom = Atom::new("CA", gly_id, Point3::new(1.4, 0.0, 0.0));

            let gly_n_id = system.add_atom_to_residue(gly_id, gly_n_atom).unwrap();
            let gly_ca_id = system.add_atom_to_residue(gly_id, gly_ca_atom).unwrap();
            system
                .add_bond(gly_n_id, gly_ca_id, BondOrder::Single)
                .unwrap();

            let ala_id = system
                .add_residue(chain_a_id, 2, "ALA", Some(ResidueType::Alanine))
                .unwrap();
            let ala_ca_atom = Atom::new("CA", ala_id, Point3::new(2.0, 1.0, 0.0));
            let ala_ca_id = system.add_atom_to_residue(ala_id, ala_ca_atom).unwrap();
            system
                .add_bond(gly_ca_id, ala_ca_id, BondOrder::Single)
                .unwrap();

            let refs = TestRefs {
                chain_a_id,
                gly_id,
                gly_n_id,
                gly_ca_id,
                ala_id,
                ala_ca_id,
            };

            (system, refs)
        }

        fn create_system_with_duplicate_atoms() -> MolecularSystem {
            let mut system = MolecularSystem::new();
            let chain_id = system.add_chain('A', ChainType::Protein);
            let gly_id = system
                .add_residue(chain_id, 1, "GLY", Some(ResidueType::Glycine))
                .unwrap();

            let n_atom = Atom::new("N", gly_id, Point3::origin());
            let hca1_atom = Atom::new("HCA", gly_id, Point3::origin());
            let hca2_atom = Atom::new("HCA", gly_id, Point3::origin());

            system.add_atom_to_residue(gly_id, n_atom).unwrap();
            system.add_atom_to_residue(gly_id, hca1_atom).unwrap();
            system.add_atom_to_residue(gly_id, hca2_atom).unwrap();
            system
        }

        struct DisulfideTestRefs {
            cys1_a_sg_id: AtomId,
            cys3_a_sg_id: AtomId,
            cys1_b_id: ResidueId,
            cys1_b_sg_id: AtomId,
        }

        fn create_disulfide_test_system() -> (MolecularSystem, DisulfideTestRefs) {
            let mut system = MolecularSystem::new();

            let chain_a_id = system.add_chain('A', ChainType::Protein);
            let cys1_a_id = system
                .add_residue(chain_a_id, 1, "CYS", Some(ResidueType::Cysteine))
                .unwrap();
            system
                .add_residue(chain_a_id, 2, "GLY", Some(ResidueType::Glycine))
                .unwrap();
            let cys3_a_id = system
                .add_residue(chain_a_id, 3, "CYS", Some(ResidueType::Cysteine))
                .unwrap();

            let cys1_a_sg = Atom::new("SG", cys1_a_id, Point3::new(0.0, 0.0, 0.0));
            let cys3_a_sg = Atom::new("SG", cys3_a_id, Point3::new(0.0, 2.0, 0.0));

            let cys1_a_sg_id = system.add_atom_to_residue(cys1_a_id, cys1_a_sg).unwrap();
            let cys3_a_sg_id = system.add_atom_to_residue(cys3_a_id, cys3_a_sg).unwrap();

            let chain_b_id = system.add_chain('B', ChainType::Protein);
            let cys1_b_id = system
                .add_residue(chain_b_id, 1, "CYS", Some(ResidueType::Cysteine))
                .unwrap();
            let cys1_b_sg = Atom::new("SG", cys1_b_id, Point3::new(10.0, 0.0, 0.0));
            let cys1_b_sg_id = system.add_atom_to_residue(cys1_b_id, cys1_b_sg).unwrap();

            let refs = DisulfideTestRefs {
                cys1_a_sg_id,
                cys3_a_sg_id,
                cys1_b_id,
                cys1_b_sg_id,
            };
            (system, refs)
        }

        #[test]
        fn system_creation_and_access() {
            let (system, refs) = create_standard_test_system();

            assert_eq!(system.atoms_iter().count(), 3);
            assert_eq!(system.residues_iter().count(), 2);
            assert_eq!(system.chains_iter().count(), 1);
            assert_eq!(system.bonds().len(), 2);
            assert!(system.find_chain_by_id('B').is_none());

            let found_gly = system.find_residue_by_id(refs.chain_a_id, 1).unwrap();
            let found_ala = system.find_residue_by_id(refs.chain_a_id, 2).unwrap();
            assert_eq!(found_gly, refs.gly_id);
            assert_eq!(found_ala, refs.ala_id);

            assert_eq!(system.residue(refs.gly_id).unwrap().name, "GLY");
            assert_eq!(system.atom(refs.gly_n_id).unwrap().name, "N");
        }

        #[test]
        fn atom_removal_updates_system_correctly() {
            let (mut system, refs) = create_standard_test_system();

            assert_eq!(system.bonds().len(), 2);
            assert!(system.bonds().iter().any(|b| b.contains(refs.gly_n_id)));
            assert_eq!(system.residue(refs.gly_id).unwrap().atoms().len(), 2);
            assert!(
                system
                    .get_bonded_neighbors(refs.gly_ca_id)
                    .unwrap()
                    .contains(&refs.gly_n_id)
            );

            let removed_atom = system.remove_atom(refs.gly_n_id).unwrap();

            assert_eq!(removed_atom.name, "N");
            assert_eq!(system.atoms_iter().count(), 2);
            assert!(system.atom(refs.gly_n_id).is_none());
            assert_eq!(system.bonds().len(), 1);
            assert!(!system.bond_adjacency.contains_key(refs.gly_n_id));
            assert!(
                !system
                    .get_bonded_neighbors(refs.gly_ca_id)
                    .unwrap()
                    .contains(&refs.gly_n_id)
            );
            assert_eq!(system.residue(refs.gly_id).unwrap().atoms().len(), 1);
        }

        #[test]
        fn atom_removal_handles_duplicates_correctly() {
            let mut system = create_system_with_duplicate_atoms();
            let gly_id = system
                .find_residue_by_id(system.find_chain_by_id('A').unwrap(), 1)
                .unwrap();
            let gly_res = system.residue(gly_id).unwrap();

            let hca_ids = gly_res.get_atom_ids_by_name("HCA").unwrap().to_vec();
            assert_eq!(hca_ids.len(), 2);
            let id_to_remove = hca_ids[0];
            let id_to_keep = hca_ids[1];

            system.remove_atom(id_to_remove);

            let remaining_hca_ids = system
                .residue(gly_id)
                .unwrap()
                .get_atom_ids_by_name("HCA")
                .unwrap();
            assert_eq!(remaining_hca_ids.len(), 1);
            assert_eq!(remaining_hca_ids[0], id_to_keep);
            assert_eq!(system.residue(gly_id).unwrap().atoms().len(), 2);
        }

        #[test]
        fn residue_removal_updates_system_correctly() {
            let (mut system, refs) = create_standard_test_system();

            assert_eq!(system.atoms_iter().count(), 3);
            assert_eq!(system.residues_iter().count(), 2);
            assert!(system.find_residue_by_id(refs.chain_a_id, 1).is_some());

            let removed_residue = system.remove_residue(refs.gly_id).unwrap();

            assert_eq!(removed_residue.name, "GLY");
            assert_eq!(system.residues_iter().count(), 1);
            assert!(system.residue(refs.gly_id).is_none());
            assert!(system.find_residue_by_id(refs.chain_a_id, 1).is_none());
            assert_eq!(system.atoms_iter().count(), 1);
            assert!(system.atom(refs.gly_n_id).is_none());
            assert!(system.atom(refs.gly_ca_id).is_none());
            assert!(system.atom(refs.ala_ca_id).is_some());
            assert!(
                system.bonds().is_empty(),
                "Bonds to removed atoms should be gone"
            );
            assert_eq!(system.chain(refs.chain_a_id).unwrap().residues().len(), 1);
        }

        #[test]
        fn get_bonded_neighbors_returns_correct_neighbors() {
            let (system, refs) = create_standard_test_system();
            let _gly_residue = system.residue(refs.gly_id).unwrap();

            let n_neighbors = system.get_bonded_neighbors(refs.gly_n_id).unwrap();
            assert_eq!(n_neighbors, &[refs.gly_ca_id]);

            let ca_neighbors = system.get_bonded_neighbors(refs.gly_ca_id).unwrap();
            assert_eq!(ca_neighbors.len(), 2);
            assert!(ca_neighbors.contains(&refs.gly_n_id));
            assert!(ca_neighbors.contains(&refs.ala_ca_id));

            let ala_ca_neighbors = system.get_bonded_neighbors(refs.ala_ca_id).unwrap();
            assert_eq!(ala_ca_neighbors, &[refs.gly_ca_id]);
        }

        #[test]
        fn handles_intrachain_disulfide_bond_correctly() {
            let (mut system, refs) = create_disulfide_test_system();
            system
                .add_bond(refs.cys1_a_sg_id, refs.cys3_a_sg_id, BondOrder::Single)
                .unwrap();

            assert_eq!(
                system.get_bonded_neighbors(refs.cys1_a_sg_id).unwrap(),
                &[refs.cys3_a_sg_id]
            );
            assert_eq!(
                system.get_bonded_neighbors(refs.cys3_a_sg_id).unwrap(),
                &[refs.cys1_a_sg_id]
            );

            system.remove_atom(refs.cys1_a_sg_id);

            assert!(system.bonds().is_empty());
            assert!(system.get_bonded_neighbors(refs.cys1_a_sg_id).is_none());
            assert!(
                system
                    .get_bonded_neighbors(refs.cys3_a_sg_id)
                    .unwrap()
                    .is_empty()
            );
        }

        #[test]
        fn handles_interchain_disulfide_bond_correctly() {
            let (mut system, refs) = create_disulfide_test_system();
            system
                .add_bond(refs.cys1_a_sg_id, refs.cys1_b_sg_id, BondOrder::Single)
                .unwrap();

            assert_eq!(
                system.get_bonded_neighbors(refs.cys1_a_sg_id).unwrap(),
                &[refs.cys1_b_sg_id]
            );
            assert_eq!(
                system.get_bonded_neighbors(refs.cys1_b_sg_id).unwrap(),
                &[refs.cys1_a_sg_id]
            );

            system.remove_residue(refs.cys1_b_id);

            assert!(system.bonds().is_empty());
            assert!(system.atom(refs.cys1_b_sg_id).is_none());
            assert!(system.residue(refs.cys1_b_id).is_none());
            assert!(
                system
                    .get_bonded_neighbors(refs.cys1_a_sg_id)
                    .unwrap()
                    .is_empty()
            );
        }

        #[test]
        fn idempotent_add_bond_does_not_create_duplicates() {
            let (mut system, refs) = create_disulfide_test_system();
            system
                .add_bond(refs.cys1_a_sg_id, refs.cys3_a_sg_id, BondOrder::Single)
                .unwrap();
            system
                .add_bond(refs.cys3_a_sg_id, refs.cys1_a_sg_id, BondOrder::Single)
                .unwrap();

            assert_eq!(
                system.bonds().len(),
                1,
                "Adding an existing bond should be idempotent"
            );
            let neighbors = system.get_bonded_neighbors(refs.cys1_a_sg_id).unwrap();
            assert_eq!(
                neighbors.len(),
                1,
                "Adjacency list should not contain duplicates"
            );
        }
    }

    mod role_based_queries {
        use super::*;
        use std::collections::HashMap;

        fn create_system_with_roles() -> (MolecularSystem, HashMap<&'static str, AtomId>) {
            let mut system = MolecularSystem::new();
            let mut id_map = HashMap::new();

            let chain_a_id = system.add_chain('A', ChainType::Protein);
            let ala_id = system
                .add_residue(chain_a_id, 1, "ALA", Some(ResidueType::Alanine))
                .unwrap();

            let mut ala_n = Atom::new("N", ala_id, Point3::origin());
            ala_n.role = AtomRole::Backbone;
            id_map.insert("ALA_N", system.add_atom_to_residue(ala_id, ala_n).unwrap());

            let mut ala_ca = Atom::new("CA", ala_id, Point3::origin());
            ala_ca.role = AtomRole::Backbone;
            id_map.insert(
                "ALA_CA",
                system.add_atom_to_residue(ala_id, ala_ca).unwrap(),
            );

            let mut ala_cb = Atom::new("CB", ala_id, Point3::origin());
            ala_cb.role = AtomRole::Sidechain;
            id_map.insert(
                "ALA_CB",
                system.add_atom_to_residue(ala_id, ala_cb).unwrap(),
            );

            let chain_l_id = system.add_chain('L', ChainType::Ligand);
            let lig_id = system.add_residue(chain_l_id, 101, "LIG", None).unwrap();

            let mut lig_c1 = Atom::new("C1", lig_id, Point3::origin());
            lig_c1.role = AtomRole::Ligand;
            id_map.insert(
                "LIG_C1",
                system.add_atom_to_residue(lig_id, lig_c1).unwrap(),
            );

            let chain_w_id = system.add_chain('W', ChainType::Water);
            let hoh_id = system.add_residue(chain_w_id, 201, "HOH", None).unwrap();

            let mut hoh_o = Atom::new("O", hoh_id, Point3::origin());
            hoh_o.role = AtomRole::Water;
            id_map.insert("HOH_O", system.add_atom_to_residue(hoh_id, hoh_o).unwrap());

            let unknown_atom = Atom::new("X", ala_id, Point3::origin());
            id_map.insert(
                "UNKNOWN",
                system.add_atom_to_residue(ala_id, unknown_atom).unwrap(),
            );

            (system, id_map)
        }

        #[test]
        fn atoms_by_role_returns_correct_atoms() {
            let (system, id_map) = create_system_with_roles();

            let backbone_atoms: Vec<AtomId> = system
                .atoms_by_role(AtomRole::Backbone)
                .map(|(id, _)| id)
                .collect();
            assert_eq!(backbone_atoms.len(), 2);
            assert!(backbone_atoms.contains(id_map.get("ALA_N").unwrap()));
            assert!(backbone_atoms.contains(id_map.get("ALA_CA").unwrap()));

            let sidechain_atoms: Vec<AtomId> = system
                .atoms_by_role(AtomRole::Sidechain)
                .map(|(id, _)| id)
                .collect();
            assert_eq!(sidechain_atoms, vec![*id_map.get("ALA_CB").unwrap()]);

            let ligand_atoms: Vec<AtomId> = system
                .atoms_by_role(AtomRole::Ligand)
                .map(|(id, _)| id)
                .collect();
            assert_eq!(ligand_atoms, vec![*id_map.get("LIG_C1").unwrap()]);

            let unknown_atoms: Vec<AtomId> = system
                .atoms_by_role(AtomRole::Other)
                .map(|(id, _)| id)
                .collect();
            assert_eq!(unknown_atoms, vec![*id_map.get("UNKNOWN").unwrap()]);
        }

        #[test]
        fn protein_atoms_returns_backbone_and_sidechain() {
            let (system, id_map) = create_system_with_roles();
            let protein_ids: Vec<AtomId> = system.protein_atom_ids();

            assert_eq!(protein_ids.len(), 3);
            assert!(protein_ids.contains(id_map.get("ALA_N").unwrap()));
            assert!(protein_ids.contains(id_map.get("ALA_CA").unwrap()));
            assert!(protein_ids.contains(id_map.get("ALA_CB").unwrap()));
        }

        #[test]
        fn background_atoms_returns_all_non_protein_atoms() {
            let (system, id_map) = create_system_with_roles();
            let background_ids: Vec<AtomId> = system.background_atom_ids();

            assert_eq!(background_ids.len(), 3);
            assert!(background_ids.contains(id_map.get("LIG_C1").unwrap()));
            assert!(background_ids.contains(id_map.get("HOH_O").unwrap()));
            assert!(background_ids.contains(id_map.get("UNKNOWN").unwrap()));
        }
    }
}
