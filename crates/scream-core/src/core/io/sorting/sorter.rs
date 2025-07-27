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
