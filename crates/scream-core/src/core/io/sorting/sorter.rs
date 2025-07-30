use crate::core::io::sorting::rules::{ATOM_NAME_ALIASES, ATOM_ORDER_WEIGHTS};
use crate::core::models::atom::Atom;
use crate::core::models::ids::AtomId;
use crate::core::models::system::MolecularSystem;
use std::cmp::Ordering;

#[derive(Debug)]
pub struct CanonicalAtom<'a> {
    pub id: AtomId,
    pub source: &'a Atom,
    pub chain_char: char,
    pub residue_number: isize,
}

pub fn sort_system_atoms(system: &MolecularSystem) -> Vec<CanonicalAtom> {
    let mut atoms_to_sort: Vec<CanonicalAtom> = system
        .atoms_iter()
        .map(|(atom_id, atom)| {
            let residue = system
                .residue(atom.residue_id)
                .expect("Atom must belong to a residue");
            let chain = system
                .chain(residue.chain_id)
                .expect("Residue must belong to a chain");
            CanonicalAtom {
                id: atom_id,
                source: atom,
                chain_char: chain.id,
                residue_number: residue.residue_number,
            }
        })
        .collect();

    atoms_to_sort.sort_unstable_by(|a, b| {
        // Level 1: Compare by chain character.
        a.chain_char
            .cmp(&b.chain_char)
            // Level 2: If chains are the same, compare by residue sequence number.
            .then_with(|| a.residue_number.cmp(&b.residue_number))
            // Level 3: If residues are the same, compare by canonical atom name order.
            .then_with(|| compare_atom_names(&a.source.name, &b.source.name))
    });

    atoms_to_sort
}

fn compare_atom_names(name_a: &str, name_b: &str) -> Ordering {
    let trimmed_a = name_a.trim();
    let trimmed_b = name_b.trim();

    let canonical_a = ATOM_NAME_ALIASES
        .get(trimmed_a)
        .copied()
        .unwrap_or(trimmed_a);
    let canonical_b = ATOM_NAME_ALIASES
        .get(trimmed_b)
        .copied()
        .unwrap_or(trimmed_b);

    let weight_a = ATOM_ORDER_WEIGHTS
        .get(canonical_a)
        .copied()
        .unwrap_or(i32::MAX);
    let weight_b = ATOM_ORDER_WEIGHTS
        .get(canonical_b)
        .copied()
        .unwrap_or(i32::MAX);

    match weight_a.cmp(&weight_b) {
        Ordering::Equal => trimmed_a.cmp(trimmed_b),
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::atom::Atom;
    use crate::core::models::chain::ChainType;
    use crate::core::models::residue::ResidueType;
    use nalgebra::Point3;

    fn create_disordered_system() -> MolecularSystem {
        let mut system = MolecularSystem::new();

        let chain_b_id = system.add_chain('B', ChainType::Protein);
        let chain_a_id = system.add_chain('A', ChainType::Protein);

        let res_a2_id = system
            .add_residue(chain_a_id, 2, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let res_a1_id = system
            .add_residue(chain_a_id, 1, "GLY", Some(ResidueType::Glycine))
            .unwrap();
        let res_b1_id = system
            .add_residue(chain_b_id, 1, "SER", Some(ResidueType::Serine))
            .unwrap();

        system
            .add_atom_to_residue(res_a1_id, Atom::new("CA", res_a1_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(res_a1_id, Atom::new("N", res_a1_id, Point3::origin()))
            .unwrap();

        system
            .add_atom_to_residue(res_a2_id, Atom::new("CB", res_a2_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(res_a2_id, Atom::new("CA", res_a2_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(res_a2_id, Atom::new("HCA", res_a2_id, Point3::origin()))
            .unwrap();

        system
            .add_atom_to_residue(res_b1_id, Atom::new("OG", res_b1_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(res_b1_id, Atom::new("CA", res_b1_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(res_b1_id, Atom::new("XXX", res_b1_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(res_b1_id, Atom::new("YYY", res_b1_id, Point3::origin()))
            .unwrap();

        system
    }

    #[test]
    fn test_sort_system_atoms_produces_canonical_order() {
        let system = create_disordered_system();
        let sorted_atoms = sort_system_atoms(&system);

        let sorted_names: Vec<_> = sorted_atoms
            .iter()
            .map(|ca| ca.source.name.as_str())
            .collect();

        let res_a1_id = system
            .find_residue_by_id(system.find_chain_by_id('A').unwrap(), 1)
            .unwrap();
        let res_a2_id = system
            .find_residue_by_id(system.find_chain_by_id('A').unwrap(), 2)
            .unwrap();
        let res_b1_id = system
            .find_residue_by_id(system.find_chain_by_id('B').unwrap(), 1)
            .unwrap();

        let res_a1 = system.residue(res_a1_id).unwrap();
        let res_a2 = system.residue(res_a2_id).unwrap();
        let res_b1 = system.residue(res_b1_id).unwrap();

        let mut expected_a1_names: Vec<_> = res_a1
            .atoms()
            .iter()
            .map(|id| system.atom(*id).unwrap().name.as_str())
            .collect();
        expected_a1_names.sort_by(|a, b| compare_atom_names(a, b));

        let mut expected_a2_names: Vec<_> = res_a2
            .atoms()
            .iter()
            .map(|id| system.atom(*id).unwrap().name.as_str())
            .collect();
        expected_a2_names.sort_by(|a, b| compare_atom_names(a, b));

        let mut expected_b1_names: Vec<_> = res_b1
            .atoms()
            .iter()
            .map(|id| system.atom(*id).unwrap().name.as_str())
            .collect();
        expected_b1_names.sort_by(|a, b| compare_atom_names(a, b));

        let final_expected_order: Vec<&str> = expected_a1_names
            .into_iter()
            .chain(expected_a2_names.into_iter())
            .chain(expected_b1_names.into_iter())
            .collect();

        assert_eq!(sorted_names, final_expected_order);

        assert_eq!(sorted_atoms[0].chain_char, 'A');
        assert_eq!(sorted_atoms[0].residue_number, 1);
        assert_eq!(sorted_atoms[0].source.name, "N");

        assert_eq!(sorted_atoms[1].chain_char, 'A');
        assert_eq!(sorted_atoms[1].residue_number, 1);
        assert_eq!(sorted_atoms[1].source.name, "CA");

        assert_eq!(sorted_atoms[2].chain_char, 'A');
        assert_eq!(sorted_atoms[2].residue_number, 2);
        assert_eq!(sorted_atoms[2].source.name, "CA");

        assert_eq!(sorted_atoms[5].chain_char, 'B');
        assert_eq!(sorted_atoms[5].residue_number, 1);
    }

    #[test]
    fn test_compare_atom_names() {
        assert_eq!(compare_atom_names("N", "CA"), Ordering::Less);
        assert_eq!(compare_atom_names("C", "CA"), Ordering::Greater);

        assert_eq!(compare_atom_names("CA", "CB"), Ordering::Less);
        assert_eq!(compare_atom_names("CB", "CG"), Ordering::Less);

        assert_eq!(
            compare_atom_names("HCA", "CB"),
            Ordering::Less,
            "HCA should be treated as HA, which comes before CB"
        );
        assert_eq!(
            compare_atom_names("H", "CA"),
            Ordering::Less,
            "H should be treated as HN, which comes before CA"
        );

        assert_eq!(
            compare_atom_names("XXX", "CA"),
            Ordering::Greater,
            "Unknown names should come after known ones"
        );
        assert_eq!(
            compare_atom_names("XXX", "YYY"),
            Ordering::Less,
            "Two unknowns should be sorted alphabetically"
        );
        assert_eq!(compare_atom_names("YYY", "XXX"), Ordering::Greater);
    }
}
