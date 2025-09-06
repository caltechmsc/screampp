use super::ids::AtomId;
use std::fmt;
use std::str::FromStr;
use thiserror::Error;

/// Represents the order of a chemical bond between atoms.
///
/// This enum defines the possible bond orders used in molecular topology,
/// supporting common bond types in organic and biological molecules.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum BondOrder {
    /// Single bond.
    Single,
    /// Double bond.
    Double,
    /// Triple bond.
    Triple,
    /// Aromatic bond.
    Aromatic,
}

impl Default for BondOrder {
    /// Returns the default bond order.
    ///
    /// The default is `BondOrder::Single`, representing the most common bond type.
    fn default() -> Self {
        BondOrder::Single
    }
}

/// Error type for failed parsing of bond order strings.
///
/// This error is returned when attempting to parse an invalid
/// string into a `BondOrder`.
#[derive(Debug, Error)]
#[error("Invalid bond order string")]
pub struct ParseBondOrderError;

impl FromStr for BondOrder {
    type Err = ParseBondOrderError;

    /// Parses a string into a `BondOrder`.
    ///
    /// This implementation supports multiple string representations
    /// for each bond order, including numeric and textual forms.
    /// It is case-insensitive.
    ///
    /// # Arguments
    ///
    /// * `s` - The string to parse.
    ///
    /// # Return
    ///
    /// Returns the parsed `BondOrder` if successful.
    ///
    /// # Errors
    ///
    /// Returns `ParseBondOrderError` if the string is invalid.
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
    /// Formats the `BondOrder` as a human-readable string.
    ///
    /// This implementation allows `BondOrder` to be displayed as a string
    /// using capitalized names for each variant.
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

/// Represents a chemical bond between two atoms.
///
/// This struct defines a bond with its constituent atoms and order,
/// providing the basic topology information for molecular structures.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Bond {
    /// ID of the first atom in the bond.
    pub atom1_id: AtomId,
    /// ID of the second atom in the bond.
    pub atom2_id: AtomId,
    /// The order of the bond.
    pub order: BondOrder,
}

impl Bond {
    /// Creates a new `Bond` between two atoms.
    ///
    /// This constructor normalizes the atom IDs to ensure `atom1_id <= atom2_id`,
    /// providing a canonical representation for bond equality and hashing.
    ///
    /// # Arguments
    ///
    /// * `atom1_id` - ID of one atom in the bond.
    /// * `atom2_id` - ID of the other atom in the bond.
    /// * `order` - The order of the bond.
    pub fn new(atom1_id: AtomId, atom2_id: AtomId, order: BondOrder) -> Self {
        // Normalize atom IDs to ensure consistent ordering for equality and hashing.
        let (atom1_id, atom2_id) = if atom1_id <= atom2_id {
            (atom1_id, atom2_id)
        } else {
            (atom2_id, atom1_id)
        };
        Self {
            atom1_id,
            atom2_id,
            order,
        }
    }

    /// Checks if the bond contains a specific atom.
    ///
    /// This method determines whether the given atom ID is one of the
    /// two atoms participating in the bond.
    ///
    /// # Arguments
    ///
    /// * `atom_id` - The atom ID to check.
    ///
    /// # Return
    ///
    /// Returns `true` if the atom is part of the bond, otherwise `false`.
    pub fn contains(&self, atom_id: AtomId) -> bool {
        self.atom1_id == atom_id || self.atom2_id == atom_id
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use slotmap::KeyData;

    fn dummy_atom_id(n: u64) -> AtomId {
        AtomId::from(KeyData::from_ffi(n))
    }

    #[test]
    fn bond_order_from_str_parses_valid_strings() {
        assert_eq!("1".parse::<BondOrder>().unwrap(), BondOrder::Single);
        assert_eq!("single".parse::<BondOrder>().unwrap(), BondOrder::Single);
        assert_eq!("S".parse::<BondOrder>().unwrap(), BondOrder::Single);
        assert_eq!("2".parse::<BondOrder>().unwrap(), BondOrder::Double);
        assert_eq!("double".parse::<BondOrder>().unwrap(), BondOrder::Double);
        assert_eq!("D".parse::<BondOrder>().unwrap(), BondOrder::Double);
        assert_eq!("3".parse::<BondOrder>().unwrap(), BondOrder::Triple);
        assert_eq!("triple".parse::<BondOrder>().unwrap(), BondOrder::Triple);
        assert_eq!("T".parse::<BondOrder>().unwrap(), BondOrder::Triple);
        assert_eq!("ar".parse::<BondOrder>().unwrap(), BondOrder::Aromatic);
        assert_eq!(
            "aromatic".parse::<BondOrder>().unwrap(),
            BondOrder::Aromatic
        );
    }

    #[test]
    fn bond_order_from_str_rejects_invalid_strings() {
        assert!("".parse::<BondOrder>().is_err());
        assert!("quadruple".parse::<BondOrder>().is_err());
        assert!("unknown".parse::<BondOrder>().is_err());
        assert!("0".parse::<BondOrder>().is_err());
    }

    #[test]
    fn bond_order_display_outputs_expected_strings() {
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
    fn bond_new_initializes_fields_correctly() {
        let a1 = dummy_atom_id(1);
        let a2 = dummy_atom_id(2);
        let bond = Bond::new(a1, a2, BondOrder::Double);
        assert_eq!(bond.atom1_id, a1);
        assert_eq!(bond.atom2_id, a2);
        assert_eq!(bond.order, BondOrder::Double);
    }

    #[test]
    fn bond_contains_returns_true_for_both_atoms() {
        let a1 = dummy_atom_id(10);
        let a2 = dummy_atom_id(20);
        let bond = Bond::new(a1, a2, BondOrder::Single);
        assert!(bond.contains(a1));
        assert!(bond.contains(a2));
    }

    #[test]
    fn bond_contains_returns_false_for_unrelated_atom() {
        let a1 = dummy_atom_id(100);
        let a2 = dummy_atom_id(200);
        let unrelated = dummy_atom_id(300);
        let bond = Bond::new(a1, a2, BondOrder::Aromatic);
        assert!(!bond.contains(unrelated));
    }
}
