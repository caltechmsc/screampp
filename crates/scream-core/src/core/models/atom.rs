use nalgebra::Point3;

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub index: usize,             // Index in the global atom vector
    pub serial: usize,            // Original serial from source file (e.g., PDB)
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
    fn test_atom_creation() {
        let atom = Atom {
            index: 0,
            serial: 1,
            name: "CA".to_string(),
            force_field_type: "C.3".to_string(),
            partial_charge: 0.1,
            position: Point3::new(1.0, 2.0, 3.0),
        };

        assert_eq!(atom.index, 0);
        assert_eq!(atom.serial, 1);
        assert_eq!(atom.name, "CA");
        assert_eq!(atom.force_field_type, "C.3");
        assert_eq!(atom.partial_charge, 0.1);
        assert_eq!(atom.position, Point3::new(1.0, 2.0, 3.0));
    }
}
