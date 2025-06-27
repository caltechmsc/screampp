use nalgebra::Point3;

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

impl Element {
    pub fn from_str(symbol: &str) -> Self {
        match symbol.trim().to_uppercase().as_str() {
            "H" | "1H" | "D" | "2H" | "T" | "3H" => Self::H,
            "C" => Self::C,
            "N" => Self::N,
            "O" => Self::O,
            "P" => Self::P,
            "S" => Self::S,

            "F" => Self::F,
            "CL" => Self::Cl,
            "BR" => Self::Br,
            "I" => Self::I,

            "NA" => Self::Na,
            "K" => Self::K,
            "MG" => Self::Mg,
            "CA" => Self::Ca,

            "MN" => Self::Mn,
            "FE" => Self::Fe,
            "CO" => Self::Co,
            "NI" => Self::Ni,
            "CU" => Self::Cu,
            "ZN" => Self::Zn,
            "MO" => Self::Mo,
            "W" => Self::W,
            "V" => Self::V,
            "CD" => Self::Cd,
            "HG" => Self::Hg,
            "PT" => Self::Pt,
            "AU" => Self::Au,
            "PD" => Self::Pd,
            "RU" => Self::Ru,
            "RH" => Self::Rh,
            "IR" => Self::Ir,

            "AL" => Self::Al,
            "SI" => Self::Si,
            "SE" => Self::Se,
            "AS" => Self::As,
            "B" => Self::B,
            "LI" => Self::Li,

            "LP" => Self::Lp,
            "DU" | "DUM" | "DUMMY" => Self::Du,
            "X" => Self::Unknown,

            _ => Self::Unknown,
        }
    }

    pub fn to_str(&self) -> &'static str {
        match self {
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
        }
    }
}

pub struct Atom {
    pub serial: usize,            // Unique identifier for the atom
    pub element: Element,         // Chemical element of the atom
    pub name: String,             // Atom name (e.g., "CA" for alpha carbon)
    pub force_field_type: String, // Force field type (e.g., "C.3" for sp3 carbon)
    pub partial_charge: f64,      // Partial charge of the atom
    pub position: Point3<f64>,    // 3D coordinates of the atom
}

impl Atom {
    pub fn new(
        serial: usize,
        element: Element,
        name: String,
        force_field_type: String,
        partial_charge: f64,
        position: Point3<f64>,
    ) -> Self {
        Self {
            serial,
            element,
            name,
            force_field_type,
            partial_charge,
            position,
        }
    }
}
