use phf::{Map, Set, phf_map, phf_set};
use std::cmp::Ordering;
use thiserror::Error;

static BACKBONE_ATOM_NAMES: Set<&'static str> = phf_set! {
    // Standard backbone
    "N", "H", "HN", "CA", "HA", "C", "O",
    // N-Terminus variants
    "H1", "H2", "H3", "NT", "HT1", "HT2", "HT3",
    // C-Terminus variants
    "OXT", "OT1", "OT2", "HC", "HOXT",
    // Glycine alpha-hydrogens
    "HA1", "HA2", "1HA", "2HA",
};

static STANDARD_AMINO_ACID_NAMES: Set<&'static str> = phf_set! {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "HSE", "HSP", "HSD", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL", "CYX", "LYN", "ARN", "APP", "GLP",
};

static THREE_TO_ONE_CODE: Map<&'static str, char> = phf_map! {
    "ALA" => 'A', "ARG" => 'R', "ASN" => 'N', "ASP" => 'D',
    "CYS" => 'C', "GLN" => 'Q', "GLU" => 'E', "GLY" => 'G',
    "HIS" => 'H', "ILE" => 'I', "LEU" => 'L', "LYS" => 'K',
    "MET" => 'M', "PHE" => 'F', "PRO" => 'P', "SER" => 'S',
    "THR" => 'T', "TRP" => 'W', "TYR" => 'Y', "VAL" => 'V',
    "CYX" => 'C', "HSE" => 'H', "HSP" => 'B', "HSD" => 'J',
    "LYN" => 'K', "ARN" => 'R', "APP" => 'A', "GLP" => 'Q',
};

static ONE_TO_THREE_CODE: Map<char, &'static str> = phf_map! {
    'A' => "ALA", 'R' => "ARG", 'N' => "ASN", 'D' => "ASP",
    'C' => "CYS", 'Q' => "GLN", 'E' => "GLU", 'G' => "GLY",
    'H' => "HIS", 'I' => "ILE", 'L' => "LEU", 'K' => "LYS",
    'M' => "MET", 'F' => "PHE", 'P' => "PRO", 'S' => "SER",
    'T' => "THR", 'W' => "TRP", 'Y' => "TYR", 'V' => "VAL",
};

static ATOM_ORDER_WEIGHTS: Map<&'static str, i32> = phf_map! {
    "N"    => 10, "NT"   => 10,
    "HN"   => 20, "H1"   => 20, "H2"   => 21, "H3"   => 22,
    "HN1"  => 20, "HN2"  => 21, "HN3"  => 22,
    "HT1"  => 20, "HT2"  => 21, "HT3"  => 22,
    "CA"   => 30,
    "HA"   => 40, "HCA"  => 40, "HA1"  => 41, "HA2"  => 42, "1HA" => 41, "2HA" => 42,
    "CB"   => 110,
    "HB"   => 120, "HCB"  => 120, "HB1"  => 121, "HB2"  => 122, "HB3"  => 123,
    "CG"   => 210, "CG1"  => 211, "CG2"  => 212, "SG"   => 213, "OG"   => 214, "OG1"  => 215,
    "HG"   => 220, "HCG"  => 220, "HG1"  => 221, "HG2"  => 222, "HSG"  => 223, "HOG"  => 224, "HOG1" => 225,
    "CD"   => 310, "CD1"  => 311, "CD2"  => 312, "SD"   => 313, "OD1"  => 314, "OD2"  => 315, "ND1"  => 316, "ND2"  => 317,
    "HD1"  => 321, "HD2"  => 322, "HCD"  => 320, "HCD1" => 321, "HCD2" => 322, "HOD1" => 324, "HND1" => 326, "HND2" => 327,
    "CE"   => 410, "CE1"  => 411, "CE2"  => 412, "CE3"  => 413, "NE"   => 414, "NE1"  => 415, "NE2"  => 416, "OE1"  => 417, "OE2"  => 418,
    "HE1"  => 421, "HE2"  => 422, "HE3"  => 423, "HCE"  => 420, "HCE1" => 421, "HCE2" => 422, "HCE3" => 423, "HNE"  => 424, "HNE1" => 425, "HNE2" => 426, "HOE1" => 427,
    "CZ"   => 510, "CZ2"  => 512, "CZ3"  => 513, "NZ"   => 514, "OH"   => 515,
    "HZ"   => 520, "HCZ"  => 520, "HH"   => 525, "HOH"  => 525,
    "CH2"  => 610,
    "HH2"  => 620, "HCH2" => 620,
    "NH1"  => 710, "NH2"  => 711,
    "HH11" => 721, "HH12" => 722, "HH21" => 723, "HH22" => 724,
    "C"    => 8000,
    "O"    => 9000,
    "OXT"  => 9001, "OT1"  => 9001, "OT2"  => 9002,
};

pub fn is_backbone_atom(atom_name: &str) -> bool {
    BACKBONE_ATOM_NAMES.contains(atom_name.trim())
}

pub fn is_sidechain_atom(atom_name: &str) -> bool {
    !is_backbone_atom(atom_name)
}

pub fn is_heavy_atom(atom_name: &str) -> bool {
    !atom_name.trim().starts_with('H')
}

pub fn is_standard_amino_acid(residue_name: &str) -> bool {
    STANDARD_AMINO_ACID_NAMES.contains(residue_name.trim())
}

pub fn three_to_one_letter_code(three_letter: &str) -> Option<char> {
    THREE_TO_ONE_CODE.get(three_letter.trim()).copied()
}

pub fn one_to_three_letter_code(one_letter: char) -> Option<&'static str> {
    ONE_TO_THREE_CODE.get(&one_letter).copied()
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

#[derive(Debug, Error, PartialEq, Eq)]
pub enum ParseMutationInfoError {
    #[error("Invalid format for mutation info string: '{0}'")]
    InvalidFormat(String),
    #[error("Residue ID part is not a valid number in string: '{0}'")]
    InvalidResidueId(String),
}

pub fn parse_mutation_info(s: &str) -> Result<(char, isize, char), ParseMutationInfoError> {
    let trimmed = s.trim();
    if trimmed.len() < 4 || !trimmed.contains('_') {
        return Err(ParseMutationInfoError::InvalidFormat(s.to_string()));
    }

    let aa = trimmed.chars().next().unwrap();
    let chain_id = trimmed.chars().last().unwrap();

    let parts: Vec<_> = trimmed[1..trimmed.len() - 1].split('_').collect();
    if parts.len() != 2 || !parts[1].is_empty() {
        return Err(ParseMutationInfoError::InvalidFormat(s.to_string()));
    }

    let res_id: isize = parts[0]
        .parse()
        .map_err(|_| ParseMutationInfoError::InvalidResidueId(parts[0].to_string()))?;

    Ok((aa, res_id, chain_id))
}
