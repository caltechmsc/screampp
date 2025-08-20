use super::ids::ResidueId;
use nalgebra::Point3;
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub enum AtomRole {
    Backbone,  // Backbone atom (e.g., C, N, O)
    Sidechain, // Sidechain atom (e.g., CH3, OH)
    Ligand,    // Ligand atom (e.g., in a small molecule)
    Water,     // Water molecule atom (e.g., H2O)
    #[default]
    Other, // Unknown or unclassified atom
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CachedVdwParam {
    LennardJones {
        radius: f64,
        well_depth: f64,
    },
    Buckingham {
        radius: f64,
        well_depth: f64,
        scale: f64,
    },
    None,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    // --- Identity & Topology ---
    pub name: String,          // Atom name (e.g., "CA", "N", "O")
    pub residue_id: ResidueId, // ID of the parent residue
    pub role: AtomRole,

    // --- Physicochemical Properties ---
    pub force_field_type: String, // Force field atom type (e.g., "C.3", "N.2")
    pub partial_charge: f64,      // Partial atomic charge
    pub position: Point3<f64>,    // 3D coordinates

    // --- SCREAM Algorithm Specific Parameters ---
    pub delta: f64, // "Delta" value for the flat-bottom potential

    // --- Cached Force Field Parameters for Performance ---
    pub vdw_param: CachedVdwParam, // Cached van der Waals parameters
    pub hbond_type_id: i8, // Hydrogen bond type identifier (-1: None, 0: Donor Hydrogen, >0: Acceptor)
}

impl Atom {
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
}
