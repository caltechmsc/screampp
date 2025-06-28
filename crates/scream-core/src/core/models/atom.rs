use nalgebra::Point3;
use std::fmt;
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Element {
    // --- Core Bio-organic ---
    H, // Hydrogen
    C, // Carbon
    N, // Nitrogen
    O, // Oxygen
    P, // Phosphorus
    S, // Sulfur

    // --- Common Halogens ---
    F,  // Fluorine
    Cl, // Chlorine
    Br, // Bromine
    I,  // Iodine

    // --- Common Metal Ions & Metalloids ---
    // Alkali & Alkaline Earth Metals
    Na, // Sodium
    K,  // Potassium
    Mg, // Magnesium
    Ca, // Calcium

    // Transition Metals
    Mn, // Manganese
    Fe, // Iron
    Co, // Cobalt
    Ni, // Nickel
    Cu, // Copper
    Zn, // Zinc
    Mo, // Molybdenum
    W,  // Tungsten
    V,  // Vanadium
    Cd, // Cadmium
    Hg, // Mercury
    Pt, // Platinum
    Au, // Gold
    Pd, // Palladium
    Ru, // Ruthenium
    Rh, // Rhodium
    Ir, // Iridium

    // Other Metals & Metalloids
    Al, // Aluminum
    Si, // Silicon
    Se, // Selenium
    As, // Arsenic
    B,  // Boron
    Li, // Lithium

    // --- Special ---
    Lp,      // Lone Pair
    Du,      // Dummy Atom
    Unknown, // Represents a failure to parse or an unrecognized element (Placeholder for unknown elements)
}

#[derive(Debug, PartialEq, Eq)]
pub struct ParseElementError;

impl FromStr for Element {
    type Err = ParseElementError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_uppercase().as_str() {
            "H" | "1H" | "D" | "2H" | "T" | "3H" => Ok(Self::H),
            "C" => Ok(Self::C),
            "N" => Ok(Self::N),
            "O" => Ok(Self::O),
            "P" => Ok(Self::P),
            "S" => Ok(Self::S),
            "F" => Ok(Self::F),
            "CL" => Ok(Self::Cl),
            "BR" => Ok(Self::Br),
            "I" => Ok(Self::I),
            "NA" => Ok(Self::Na),
            "K" => Ok(Self::K),
            "MG" => Ok(Self::Mg),
            "CA" => Ok(Self::Ca),
            "MN" => Ok(Self::Mn),
            "FE" => Ok(Self::Fe),
            "CO" => Ok(Self::Co),
            "NI" => Ok(Self::Ni),
            "CU" => Ok(Self::Cu),
            "ZN" => Ok(Self::Zn),
            "MO" => Ok(Self::Mo),
            "W" => Ok(Self::W),
            "V" => Ok(Self::V),
            "CD" => Ok(Self::Cd),
            "HG" => Ok(Self::Hg),
            "PT" => Ok(Self::Pt),
            "AU" => Ok(Self::Au),
            "PD" => Ok(Self::Pd),
            "RU" => Ok(Self::Ru),
            "RH" => Ok(Self::Rh),
            "IR" => Ok(Self::Ir),
            "AL" => Ok(Self::Al),
            "SI" => Ok(Self::Si),
            "SE" => Ok(Self::Se),
            "AS" => Ok(Self::As),
            "B" => Ok(Self::B),
            "LI" => Ok(Self::Li),
            "LP" => Ok(Self::Lp),
            "DU" | "DUM" | "DUMMY" => Ok(Self::Du),
            "X" => Ok(Self::Unknown),
            _ => Err(ParseElementError),
        }
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let symbol = match self {
            Self::H => "H",
            Self::C => "C",
            Self::N => "N",
            Self::O => "O",
            Self::P => "P",
            Self::S => "S",
            Self::F => "F",
            Self::Cl => "Cl",
            Self::Br => "Br",
            Self::I => "I",
            Self::Na => "Na",
            Self::K => "K",
            Self::Mg => "Mg",
            Self::Ca => "Ca",
            Self::Mn => "Mn",
            Self::Fe => "Fe",
            Self::Co => "Co",
            Self::Ni => "Ni",
            Self::Cu => "Cu",
            Self::Zn => "Zn",
            Self::Mo => "Mo",
            Self::W => "W",
            Self::V => "V",
            Self::Cd => "Cd",
            Self::Hg => "Hg",
            Self::Pt => "Pt",
            Self::Au => "Au",
            Self::Pd => "Pd",
            Self::Ru => "Ru",
            Self::Rh => "Rh",
            Self::Ir => "Ir",
            Self::Al => "Al",
            Self::Si => "Si",
            Self::Se => "Se",
            Self::As => "As",
            Self::B => "B",
            Self::Li => "Li",
            Self::Lp => "Lp",
            Self::Du => "Du",
            Self::Unknown => "X",
        };
        write!(f, "{}", symbol)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub index: usize,             // Index in the global atom vector
    pub serial: usize,            // Original serial from source file (e.g., PDB)
    pub element: Element,         // Chemical element of the atom
    pub name: String,             // Atom name (e.g., "CA" for alpha carbon)
    pub force_field_type: String, // Force field type (e.g., "C.3" for sp3 carbon)
    pub partial_charge: f64,      // Partial charge of the atom
    pub position: Point3<f64>,    // 3D coordinates of the atom
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    #[test]
    fn parses_valid_element_symbols() {
        assert_eq!(Element::from_str("H").unwrap(), Element::H);
        assert_eq!(Element::from_str("C").unwrap(), Element::C);
        assert_eq!(Element::from_str("O").unwrap(), Element::O);
        assert_eq!(Element::from_str("Cl").unwrap(), Element::Cl);
        assert_eq!(Element::from_str("Fe").unwrap(), Element::Fe);
    }

    #[test]
    fn parses_isotopes_and_aliases() {
        assert_eq!(Element::from_str("1H").unwrap(), Element::H);
        assert_eq!(Element::from_str("D").unwrap(), Element::H);
        assert_eq!(Element::from_str("2H").unwrap(), Element::H);
        assert_eq!(Element::from_str("T").unwrap(), Element::H);
        assert_eq!(Element::from_str("3H").unwrap(), Element::H);
    }

    #[test]
    fn parses_unknown_element() {
        assert_eq!(Element::from_str("X").unwrap(), Element::Unknown);
    }

    #[test]
    fn fails_to_parse_invalid_element() {
        assert!(Element::from_str("Invalid").is_err());
        assert!(Element::from_str("").is_err());
        assert!(Element::from_str("123").is_err());
    }

    #[test]
    fn displays_element_correctly() {
        assert_eq!(Element::H.to_string(), "H");
        assert_eq!(Element::C.to_string(), "C");
        assert_eq!(Element::Cl.to_string(), "Cl");
        assert_eq!(Element::Unknown.to_string(), "X");
    }

    #[test]
    fn creates_atom_with_valid_data() {
        let atom = Atom {
            index: 1,
            serial: 100,
            element: Element::C,
            name: "CA".to_string(),
            force_field_type: "C.3".to_string(),
            partial_charge: -0.12,
            position: Point3::new(1.0, 2.0, 3.0),
        };

        assert_eq!(atom.index, 1);
        assert_eq!(atom.serial, 100);
        assert_eq!(atom.element, Element::C);
        assert_eq!(atom.name, "CA");
        assert_eq!(atom.force_field_type, "C.3");
        assert_eq!(atom.partial_charge, -0.12);
        assert_eq!(atom.position, Point3::new(1.0, 2.0, 3.0));
    }
}
