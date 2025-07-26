use super::ids::ResidueId;
use nalgebra::Point3;

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
    pub serial: usize,         // Atom serial number from source file
    pub name: String,          // Atom name (e.g., "CA", "N", "O")
    pub residue_id: ResidueId, // ID of the parent residue

    // --- Physicochemical Properties ---
    pub force_field_type: String, // Force field atom type (e.g., "C.3", "N.2")
    pub partial_charge: f64,      // Partial atomic charge
    pub position: Point3<f64>,    // 3D coordinates

    // --- SCREAM Algorithm Specific Parameters ---
    pub delta: f64, // "Delta" value for the flat-bottom potential

    // --- Cached Force Field Parameters for Performance ---
    pub vdw_param: CachedVdwParam, // Cached van der Waals parameters
    pub hbond_type_id: i8,         // Hydrogen bond type identifier
}

impl Atom {
    pub fn new(serial: usize, name: &str, residue_id: ResidueId, position: Point3<f64>) -> Self {
        Self {
            serial,
            name: name.to_string(),
            residue_id,
            position,
            force_field_type: "".to_string(),
            partial_charge: 0.0,
            delta: 0.0,
            vdw_param: CachedVdwParam::None,
            hbond_type_id: -1,
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
        let atom = Atom::new(42, "CA", residue_id, Point3::new(1.0, 2.0, 3.0));
        assert_eq!(atom.serial, 42);
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
    fn atom_equality_and_clone_works() {
        let residue_id = ResidueId::default();
        let atom1 = Atom::new(1, "N", residue_id, Point3::new(0.0, 0.0, 0.0));
        let atom2 = atom1.clone();
        assert_eq!(atom1, atom2);
    }
}
