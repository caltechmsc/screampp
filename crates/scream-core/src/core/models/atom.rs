use super::ids::ResidueId;
use nalgebra::Point3;
use std::str::FromStr;

/// Represents the role or classification of an atom within a molecular structure.
///
/// This enum categorizes atoms based on their functional role in the molecule,
/// which is useful for algorithms that need to distinguish between different
/// types of atoms (e.g., backbone vs. sidechain) for computational efficiency
/// or specific force field applications.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub enum AtomRole {
    /// Backbone atom, typically part of the main chain in proteins (e.g., C, N, O).
    Backbone,
    /// Sidechain atom, part of the side groups attached to the backbone.
    Sidechain,
    /// Ligand atom, associated with small molecules or ligands bound to the structure.
    Ligand,
    /// Water molecule atom, for solvent molecules in the system.
    Water,
    /// Unknown or unclassified atom role.
    #[default]
    Other,
}

/// Caches van der Waals parameters for efficient force field computations.
///
/// This enum stores pre-computed parameters for different van der Waals potentials,
/// allowing for faster energy calculations by avoiding repeated lookups.
/// It supports common potential forms used in molecular simulations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CachedVdwParam {
    /// Lennard-Jones potential parameters.
    LennardJones {
        /// The van der Waals radius in Angstroms.
        radius: f64,
        /// The well depth parameter (epsilon) in kcal/mol.
        well_depth: f64,
    },
    /// Buckingham potential parameters.
    Buckingham {
        /// The van der Waals radius in Angstroms.
        radius: f64,
        /// The well depth parameter in kcal/mol.
        well_depth: f64,
        /// Scaling factor for the exponential term.
        scale: f64,
    },
    /// No cached parameters available.
    None,
}

/// Represents an atom in a molecular structure with its properties and parameters.
///
/// This struct encapsulates all the necessary information about an atom,
/// including its identity, physicochemical properties, and SCREAM-specific
/// parameters used in side-chain placement algorithms. It is designed for
/// high-performance computations in protein modeling.
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    /// The name of the atom (e.g., "CA", "N", "O").
    pub name: String,
    /// The ID of the parent residue this atom belongs to.
    pub residue_id: ResidueId,
    /// The role or classification of the atom in the molecular structure.
    pub role: AtomRole,
    /// The force field atom type (e.g., "C.3", "N.2").
    pub force_field_type: String,
    /// The partial atomic charge in elementary charge units.
    pub partial_charge: f64,
    /// The 3D coordinates of the atom in Angstroms.
    pub position: Point3<f64>,
    /// The "Delta" value for the flat-bottom potential in SCREAM algorithm.
    pub delta: f64,
    /// Cached van der Waals parameters for performance optimization.
    pub vdw_param: CachedVdwParam,
    /// Hydrogen bond type identifier (-1: None, 0: Donor Hydrogen, >0: Acceptor).
    pub hbond_type_id: i8,
}

impl Atom {
    /// Creates a new `Atom` with default values for most fields.
    ///
    /// This constructor initializes an atom with the provided name, residue ID,
    /// and position. Other fields are set to their default values and can be
    /// modified afterward as needed.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the atom.
    /// * `residue_id` - The ID of the residue this atom belongs to.
    /// * `position` - The 3D coordinates of the atom.
    pub fn new(name: &str, residue_id: ResidueId, position: Point3<f64>) -> Self {
        Self {
            name: name.to_string(),
            residue_id,
            position,
            role: AtomRole::default(),
            force_field_type: String::new(),
            partial_charge: 0.0,
            delta: 0.0,
            vdw_param: CachedVdwParam::None,
            hbond_type_id: -1,
        }
    }
}

impl FromStr for AtomRole {
    type Err = ();

    /// Parses a string into an `AtomRole`.
    ///
    /// This implementation allows converting string representations of atom roles
    /// into the corresponding enum variants. It is case-insensitive and supports
    /// common variations (e.g., "side-chain" or "side_chain").
    ///
    /// # Arguments
    ///
    /// * `s` - The string to parse into an `AtomRole`.
    ///
    /// # Return
    ///
    /// Returns the parsed `AtomRole` if successful.
    ///
    /// # Errors
    ///
    /// Returns `()` if the input string does not match any known atom role.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "backbone" => Ok(AtomRole::Backbone),
            "sidechain" | "side-chain" | "side_chain" => Ok(AtomRole::Sidechain),
            "ligand" => Ok(AtomRole::Ligand),
            "water" => Ok(AtomRole::Water),
            "other" | "unknown" => Ok(AtomRole::Other),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::ids::ResidueId;
    use nalgebra::Point3;

    #[test]
    fn new_atom_has_expected_default_fields() {
        let residue_id = ResidueId::default();
        let atom = Atom::new("CA", residue_id, Point3::new(1.0, 2.0, 3.0));

        assert_eq!(atom.name, "CA");
        assert_eq!(atom.residue_id, residue_id);
        assert_eq!(atom.position, Point3::new(1.0, 2.0, 3.0));
        assert_eq!(atom.force_field_type, "");
        assert_eq!(atom.partial_charge, 0.0);
        assert_eq!(atom.delta, 0.0);
        assert!(matches!(atom.vdw_param, CachedVdwParam::None));
        assert_eq!(atom.hbond_type_id, -1);
    }

    #[test]
    fn new_atom_has_correct_default_role() {
        let residue_id = ResidueId::default();
        let atom = Atom::new("CA", residue_id, Point3::new(1.0, 2.0, 3.0));
        assert_eq!(atom.role, AtomRole::Other);
        assert_eq!(atom.role, AtomRole::default());
    }

    #[test]
    fn atom_equality_and_clone_works() {
        let residue_id = ResidueId::default();
        let mut atom1 = Atom::new("N", residue_id, Point3::new(0.0, 0.0, 0.0));
        atom1.role = AtomRole::Backbone; // Also test non-default fields
        let atom2 = atom1.clone();
        assert_eq!(atom1, atom2);
    }

    #[test]
    fn atom_role_derives_expected_traits() {
        let role = AtomRole::Sidechain;
        let role_clone = role.clone();
        assert_eq!(role, role_clone);
        assert_eq!(format!("{:?}", role), "Sidechain");
    }

    #[test]
    fn from_str_parses_valid_roles() {
        assert_eq!(AtomRole::from_str("backbone"), Ok(AtomRole::Backbone));
        assert_eq!(AtomRole::from_str("sidechain"), Ok(AtomRole::Sidechain));
        assert_eq!(AtomRole::from_str("side-chain"), Ok(AtomRole::Sidechain));
        assert_eq!(AtomRole::from_str("side_chain"), Ok(AtomRole::Sidechain));
        assert_eq!(AtomRole::from_str("ligand"), Ok(AtomRole::Ligand));
        assert_eq!(AtomRole::from_str("water"), Ok(AtomRole::Water));
        assert_eq!(AtomRole::from_str("other"), Ok(AtomRole::Other));
        assert_eq!(AtomRole::from_str("unknown"), Ok(AtomRole::Other));
    }

    #[test]
    fn from_str_is_case_insensitive() {
        assert_eq!(AtomRole::from_str("BACKBONE"), Ok(AtomRole::Backbone));
        assert_eq!(AtomRole::from_str("SideChain"), Ok(AtomRole::Sidechain));
        assert_eq!(AtomRole::from_str("LiGaNd"), Ok(AtomRole::Ligand));
        assert_eq!(AtomRole::from_str("wAtEr"), Ok(AtomRole::Water));
        assert_eq!(AtomRole::from_str("OtHeR"), Ok(AtomRole::Other));
    }

    #[test]
    fn from_str_returns_err_for_invalid_role() {
        assert_eq!(AtomRole::from_str("foo"), Err(()));
        assert_eq!(AtomRole::from_str("").unwrap_err(), ());
        assert_eq!(AtomRole::from_str("123"), Err(()));
        assert_eq!(AtomRole::from_str("side chainz"), Err(()));
    }
}
