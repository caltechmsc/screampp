use phf::{Map, Set, phf_map, phf_set};
use std::cmp::Ordering;

static BACKBONE_ATOM_NAMES: Set<&'static str> = phf_set! {
    "N", "H", "HN", "CA", "HA", "C", "O", "OXT", "H1", "H2", "H3", "NT",
    "HT1", "HT2", "HT3", "OT1", "OT2", "HC", "HOXT", "HA1", "HA2", "1HA", "2HA",
};

static ATOM_ORDER_WEIGHTS: Map<&'static str, i32> = phf_map! {
    "N" => 10, "NT" => 10,
    "HN" => 20, "HN1" => 21, "HN2" => 22, "HN3" => 23, "H1" => 20, "H2" => 21, "H3" => 22,
    "CA" => 30,
    "HCA" => 40, "HA" => 40, "HA1" => 41, "HA2" => 42, "1HA" => 41, "2HA" => 42,
    "CB" => 110,
    "HB" => 120, "HCB" => 120, "HB1" => 121, "HB2" => 122, "HB3" => 123,
    "CG" => 210, "OG1" => 220, "CG1" => 230, "CG2" => 250, "SG" => 270, "OG" => 290,
    "HCG" => 215, "HOG1" => 225, "HCG1" => 240, "HCG2" => 260, "HSG" => 280, "HOG" => 300,
    "CD" => 310, "ND1" => 330, "CD1" => 350, "CD2" => 370, "OD1" => 390, "OD2" => 400, "ND2" => 410, "SD" => 430,
    "HCD" => 320, "HND1" => 340, "HCD1" => 360, "HCD2" => 380, "HOD1" => 391, "HOD2" => 401, "HND2" => 420, "HSD" => 431,
    "CE" => 510, "NE" => 530, "NE1" => 550, "CE1" => 570, "OE1" => 590, "OE2" => 600, "NE2" => 610, "CE2" => 630, "CE3" => 650,
    "HCE" => 520, "HNE" => 540, "HNE1" => 560, "HCE1" => 580, "HOE1" => 591, "HNE2" => 620, "HCE2" => 640, "HCE3" => 660,
    "CZ" => 810, "NZ" => 830, "CZ2" => 850, "CZ3" => 870, "OH" => 880,
    "HCZ" => 820, "HNZ" => 840, "HCZ2" => 860, "HCZ3" => 880, "HOH" => 890,
    "CH2" => 910, "NH1" => 950, "NH2" => 970,
    "HCH2" => 920, "HNH1" => 960, "HH11" => 961, "HH12" => 962, "HNH2" => 980, "HH21" => 981, "HH22" => 982,
    "C" => 5000, "HC" => 5001,
    "O" => 6000, "OX" => 6001, "OXT" => 6002, "HOXT" => 6003,
};

pub fn is_backbone_atom(atom_name: &str) -> bool {
    BACKBONE_ATOM_NAMES.contains(atom_name.trim())
}

pub fn is_heavy_atom(atom_name: &str) -> bool {
    let first_char = atom_name
        .trim()
        .chars()
        .next()
        .map(|c| c.to_ascii_uppercase());
    !matches!(first_char, Some('H') | Some('D'))
}

pub fn residue_atom_order(atom1_name: &str, atom2_name: &str) -> Ordering {
    let weight1 = ATOM_ORDER_WEIGHTS
        .get(atom1_name.trim())
        .unwrap_or(&i32::MAX);
    let weight2 = ATOM_ORDER_WEIGHTS
        .get(atom2_name.trim())
        .unwrap_or(&i32::MAX);
    weight1.cmp(weight2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_backbone_atom_recognizes_standard_backbone_atoms() {
        assert!(is_backbone_atom("N"));
        assert!(is_backbone_atom("CA"));
        assert!(is_backbone_atom("C"));
        assert!(is_backbone_atom("O"));
        assert!(is_backbone_atom("OXT"));
    }

    #[test]
    fn is_backbone_atom_is_case_sensitive_and_trims_whitespace() {
        assert!(!is_backbone_atom("ca"));
        assert!(is_backbone_atom(" CA "));
        assert!(!is_backbone_atom("cb"));
    }

    #[test]
    fn is_backbone_atom_returns_false_for_non_backbone_atoms() {
        assert!(!is_backbone_atom("CB"));
        assert!(!is_backbone_atom("SG"));
        assert!(!is_backbone_atom(""));
    }

    #[test]
    fn is_heavy_atom_returns_false_for_hydrogen_and_deuterium() {
        assert!(!is_heavy_atom("H"));
        assert!(!is_heavy_atom("HA"));
        assert!(!is_heavy_atom("H1"));
        assert!(!is_heavy_atom("D"));
        assert!(!is_heavy_atom("D2"));
    }

    #[test]
    fn is_heavy_atom_returns_true_for_non_hydrogen_atoms() {
        assert!(is_heavy_atom("C"));
        assert!(is_heavy_atom("CA"));
        assert!(is_heavy_atom("O"));
        assert!(is_heavy_atom("N"));
        assert!(is_heavy_atom("SG"));
    }

    #[test]
    fn is_heavy_atom_trims_whitespace_and_is_case_sensitive() {
        assert!(is_heavy_atom(" CA "));
        assert!(!is_heavy_atom(" H "));
        assert!(is_heavy_atom("c"));
        assert!(!is_heavy_atom("h"));
    }

    #[test]
    fn residue_atom_order_returns_ordering_based_on_weights() {
        use std::cmp::Ordering::*;
        assert_eq!(residue_atom_order("N", "CA"), Less);
        assert_eq!(residue_atom_order("O", "C"), Greater);
        assert_eq!(residue_atom_order("CA", "CA"), Equal);
    }

    #[test]
    fn residue_atom_order_handles_unknown_atoms_as_last() {
        use std::cmp::Ordering::*;
        assert_eq!(residue_atom_order("N", "UNKNOWN"), Less);
        assert_eq!(residue_atom_order("UNKNOWN", "N"), Greater);
        assert_eq!(residue_atom_order("UNKNOWN1", "UNKNOWN2"), Equal);
    }

    #[test]
    fn residue_atom_order_trims_whitespace() {
        use std::cmp::Ordering::*;
        assert_eq!(residue_atom_order(" N ", " CA "), Less);
        assert_eq!(residue_atom_order(" O ", " C "), Greater);
    }
}
