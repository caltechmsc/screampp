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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_valid_bond_orders_from_str() {
        assert_eq!(BondOrder::from_str("1").unwrap(), BondOrder::Single);
        assert_eq!(BondOrder::from_str("s").unwrap(), BondOrder::Single);
        assert_eq!(BondOrder::from_str("single").unwrap(), BondOrder::Single);
        assert_eq!(BondOrder::from_str("2").unwrap(), BondOrder::Double);
        assert_eq!(BondOrder::from_str("d").unwrap(), BondOrder::Double);
        assert_eq!(BondOrder::from_str("double").unwrap(), BondOrder::Double);
        assert_eq!(BondOrder::from_str("3").unwrap(), BondOrder::Triple);
        assert_eq!(BondOrder::from_str("t").unwrap(), BondOrder::Triple);
        assert_eq!(BondOrder::from_str("triple").unwrap(), BondOrder::Triple);
        assert_eq!(BondOrder::from_str("ar").unwrap(), BondOrder::Aromatic);
        assert_eq!(
            BondOrder::from_str("aromatic").unwrap(),
            BondOrder::Aromatic
        );
    }

    #[test]
    fn parses_bond_order_case_insensitive() {
        assert_eq!(BondOrder::from_str("SINGLE").unwrap(), BondOrder::Single);
        assert_eq!(BondOrder::from_str("DoUbLe").unwrap(), BondOrder::Double);
        assert_eq!(BondOrder::from_str("TrIpLe").unwrap(), BondOrder::Triple);
        assert_eq!(
            BondOrder::from_str("AROMATIC").unwrap(),
            BondOrder::Aromatic
        );
    }

    #[test]
    fn fails_to_parse_invalid_bond_order() {
        assert!(BondOrder::from_str("quadruple").is_err());
        assert!(BondOrder::from_str("").is_err());
        assert!(BondOrder::from_str("x").is_err());
        assert!(BondOrder::from_str("0").is_err());
    }

    #[test]
    fn displays_bond_order_correctly() {
        assert_eq!(BondOrder::Single.to_string(), "Single");
        assert_eq!(BondOrder::Double.to_string(), "Double");
        assert_eq!(BondOrder::Triple.to_string(), "Triple");
        assert_eq!(BondOrder::Aromatic.to_string(), "Aromatic");
    }

    #[test]
    fn bond_order_default_is_single() {
        assert_eq!(BondOrder::default(), BondOrder::Single);
    }

    #[test]
    fn creates_bond_with_lower_index_first() {
        let bond = Bond::new(2, 5, BondOrder::Double);
        assert_eq!(bond.atom1_idx, 2);
        assert_eq!(bond.atom2_idx, 5);
        assert_eq!(bond.order, BondOrder::Double);
    }

    #[test]
    fn creates_bond_with_indices_swapped_if_needed() {
        let bond = Bond::new(7, 3, BondOrder::Aromatic);
        assert_eq!(bond.atom1_idx, 3);
        assert_eq!(bond.atom2_idx, 7);
        assert_eq!(bond.order, BondOrder::Aromatic);
    }

    #[test]
    fn creates_bond_with_equal_indices() {
        let bond = Bond::new(4, 4, BondOrder::Single);
        assert_eq!(bond.atom1_idx, 4);
        assert_eq!(bond.atom2_idx, 4);
        assert_eq!(bond.order, BondOrder::Single);
    }
}
