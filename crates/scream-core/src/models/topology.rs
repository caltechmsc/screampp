#[repr(u8)]
pub enum BondOrder {
    Single,    // Single bond
    Double,    // Double bond
    Triple,    // Triple bond
    Quadruple, // Quadruple bond
    Quintuple, // Quintuple bond
    Sextuple,  // Sextuple bond
    Aromatic,  // Aromatic bond
    Dative,    // Dative bond
    Hydrogen,  // Hydrogen bond
    Ionic,     // Ionic bond
    Unknown,   // Unspecified or unknown bond order
}

impl BondOrder {
    pub fn default() -> Self {
        BondOrder::Single
    }

    pub fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "single" => BondOrder::Single,
            "double" => BondOrder::Double,
            "triple" => BondOrder::Triple,
            "quadruple" => BondOrder::Quadruple,
            "quintuple" => BondOrder::Quintuple,
            "sextuple" => BondOrder::Sextuple,
            "aromatic" => BondOrder::Aromatic,
            "dative" => BondOrder::Dative,
            "hydrogen" => BondOrder::Hydrogen,
            "ionic" => BondOrder::Ionic,
            _ => BondOrder::Unknown,
        }
    }

    pub fn to_str(&self) -> &str {
        match self {
            BondOrder::Single => "Single",
            BondOrder::Double => "Double",
            BondOrder::Triple => "Triple",
            BondOrder::Quadruple => "Quadruple",
            BondOrder::Quintuple => "Quintuple",
            BondOrder::Sextuple => "Sextuple",
            BondOrder::Aromatic => "Aromatic",
            BondOrder::Dative => "Dative",
            BondOrder::Hydrogen => "Hydrogen",
            BondOrder::Ionic => "Ionic",
            BondOrder::Unknown => "Unknown",
        }
    }
}

pub struct Bond {
    pub atom1_serial: usize, // Serial number of the first atom in the bond
    pub atom2_serial: usize, // Serial number of the second atom in the bond
    pub order: BondOrder,    // Bond order (e.g., single, double, etc.)
}

impl Bond {
    pub fn new(atom1_serial: usize, atom2_serial: usize, order: BondOrder) -> Self {
        if atom1_serial < atom2_serial {
            Self {
                atom1_serial,
                atom2_serial,
                order,
            }
        } else {
            Self {
                atom1_serial: atom2_serial,
                atom2_serial: atom1_serial,
                order,
            }
        }
    }
}
