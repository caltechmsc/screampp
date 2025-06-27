use std::fmt;
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl Default for BondOrder {
    fn default() -> Self {
        BondOrder::Single
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct ParseBondOrderError;

impl FromStr for BondOrder {
    type Err = ParseBondOrderError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "1" | "s" | "single" => Ok(Self::Single),
            "2" | "d" | "double" => Ok(Self::Double),
            "3" | "t" | "triple" => Ok(Self::Triple),
            "ar" | "aromatic" => Ok(Self::Aromatic),
            _ => Err(ParseBondOrderError),
        }
    }
}

impl fmt::Display for BondOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Single => "Single",
                Self::Double => "Double",
                Self::Triple => "Triple",
                Self::Aromatic => "Aromatic",
            }
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Bond {
    pub atom1_idx: usize, // Index of the first atom in the global atom vector
    pub atom2_idx: usize, // Index of the second atom in the global atom vector
    pub order: BondOrder, // Bond order (e.g., single, double, etc.)
}

impl Bond {
    pub fn new(atom1_idx: usize, atom2_idx: usize, order: BondOrder) -> Self {
        if atom1_idx < atom2_idx {
            Self {
                atom1_idx,
                atom2_idx,
                order,
            }
        } else {
            Self {
                atom1_idx: atom2_idx,
                atom2_idx: atom1_idx,
                order,
            }
        }
    }
}
