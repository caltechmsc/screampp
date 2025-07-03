use super::ids::AtomId;
use std::fmt;
use std::str::FromStr;
use thiserror::Error;

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

#[derive(Debug, Error)]
#[error("Invalid bond order string")]
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
    pub atom1_id: AtomId, // ID of the first atom
    pub atom2_id: AtomId, // ID of the second atom
    pub order: BondOrder, // Bond order (e.g., single, double, etc.)
}

impl Bond {
    pub fn new(atom1_id: AtomId, atom2_id: AtomId, order: BondOrder) -> Self {
        Self {
            atom1_id,
            atom2_id,
            order,
        }
    }

    pub fn contains(&self, atom_id: AtomId) -> bool {
        self.atom1_id == atom_id || self.atom2_id == atom_id
    }
}
