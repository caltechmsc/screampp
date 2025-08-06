use super::atom::{Atom, AtomRole};
use super::chain::{Chain, ChainType};
use super::ids::{AtomId, ChainId, ResidueId};
use super::residue::{Residue, ResidueType};
use super::topology::{Bond, BondOrder};
use slotmap::{SecondaryMap, SlotMap};
use std::collections::HashMap;

#[derive(Debug, Clone, Default)]
pub struct MolecularSystem {
    // --- Primary Data Stores (Source of Truth) ---
    atoms: SlotMap<AtomId, Atom>,
    residues: SlotMap<ResidueId, Residue>,
    chains: SlotMap<ChainId, Chain>,
    bonds: Vec<Bond>,

    // --- Lookup Maps for Fast Access ---
    residue_id_map: HashMap<(ChainId, isize), ResidueId>,
    chain_id_map: HashMap<char, ChainId>,

    // --- Adjacency information (cache for performance) ---
    bond_adjacency: SecondaryMap<AtomId, Vec<AtomId>>,
}

impl MolecularSystem {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn atom(&self, id: AtomId) -> Option<&Atom> {
        self.atoms.get(id)
    }
    pub fn atom_mut(&mut self, id: AtomId) -> Option<&mut Atom> {
        self.atoms.get_mut(id)
    }
    pub fn atoms_iter(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.atoms.iter()
    }

    pub fn residue(&self, id: ResidueId) -> Option<&Residue> {
        self.residues.get(id)
    }
    pub fn residue_mut(&mut self, id: ResidueId) -> Option<&mut Residue> {
        self.residues.get_mut(id)
    }
    pub fn residues_iter(&self) -> impl Iterator<Item = (ResidueId, &Residue)> {
        self.residues.iter()
    }

    pub fn chain(&self, id: ChainId) -> Option<&Chain> {
        self.chains.get(id)
    }
    pub fn chain_mut(&mut self, id: ChainId) -> Option<&mut Chain> {
        self.chains.get_mut(id)
    }
    pub fn chains_iter(&self) -> impl Iterator<Item = (ChainId, &Chain)> {
        self.chains.iter()
    }

    pub fn bonds(&self) -> &[Bond] {
        &self.bonds
    }

    pub fn find_chain_by_id(&self, id: char) -> Option<ChainId> {
        self.chain_id_map.get(&id).copied()
    }
    pub fn find_residue_by_id(
        &self,
        chain_id: ChainId,
        residue_number: isize,
    ) -> Option<ResidueId> {
        self.residue_id_map
            .get(&(chain_id, residue_number))
            .copied()
    }

    pub fn add_chain(&mut self, id: char, chain_type: ChainType) -> ChainId {
        *self.chain_id_map.entry(id).or_insert_with(|| {
            let chain = Chain::new(id, chain_type);
            self.chains.insert(chain)
        })
    }

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

    pub fn get_bonded_neighbors(&self, atom_id: AtomId) -> Option<&[AtomId]> {
        self.bond_adjacency.get(atom_id).map(|v| v.as_slice())
    }

    pub fn atoms_by_role(&self, role: AtomRole) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.atoms.iter().filter(move |(_, atom)| atom.role == role)
    }

    pub fn atom_ids_by_role(&self, role: AtomRole) -> Vec<AtomId> {
        self.atoms_by_role(role).map(|(id, _)| id).collect()
    }

    pub fn protein_atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.atoms
            .iter()
            .filter(|(_, atom)| matches!(atom.role, AtomRole::Backbone | AtomRole::Sidechain))
    }

    pub fn protein_atom_ids(&self) -> Vec<AtomId> {
        self.protein_atoms().map(|(id, _)| id).collect()
    }

    pub fn background_atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.atoms.iter().filter(|(_, atom)| {
            matches!(
                atom.role,
                AtomRole::Ligand | AtomRole::Water | AtomRole::Ion | AtomRole::Unknown
            )
        })
    }

    pub fn background_atom_ids(&self) -> Vec<AtomId> {
        self.background_atoms().map(|(id, _)| id).collect()
    }
}

#[cfg(test)]
mod tests {
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
        assert_eq!(system.residue(gly_id).unwrap().atoms.len(), 2);
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
