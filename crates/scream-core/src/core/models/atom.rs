use bitflags::bitflags;
use nalgebra::Point3;

bitflags! {
    #[derive(Default, Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub struct AtomFlags: u8 {
        const IS_FIXED_ROLE      = 0b0000_0001; // The atom's fundamental role is fixed (e.g., a backbone atom)
        const IS_TREATED_AS_FIXED = 0b0000_0010; // The atom is currently being treated as fixed in the energy expression
        const IS_VISIBLE_INTERACTION = 0b0000_0100; // The atom is visible in sidechain-sidechain interaction calculations
        const IS_VISIBLE_LATTICE     = 0b0000_1000; // The atom is visible in empty-lattice (sidechain-environment) calculations
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    // --- Identity & Topology ---
    pub index: usize,     // Global index in the system's atom list
    pub serial: usize,    // Atom serial number from source file
    pub name: String,     // Atom name (e.g., "CA", "N", "O")
    pub res_name: String, // Residue name (e.g., "ALA", "GLY")
    pub res_id: isize,    // Residue sequence number from source file
    pub chain_id: char,   // Chain identifier (e.g., 'A', 'B')

    // --- Physicochemical Properties ---
    pub force_field_type: String, // Force field atom type (e.g., "C.3", "N.2")
    pub partial_charge: f64,      // Partial atomic charge
    pub position: Point3<f64>,    // 3D coordinates

    // --- SCREAM Algorithm Specific Parameters ---
    pub flags: AtomFlags, // Bitflags for various atom properties
    pub delta: f64,       // "Delta" value for the flat-bottom potential

    // --- Cached Force Field Parameters for Performance ---
    pub vdw_radius: f64,     // van der Waals radius
    pub vdw_well_depth: f64, // van der Waals well depth (epsilon)
    pub hbond_type_id: i8,   // Hydrogen bond type identifier
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    #[test]
    fn atom_flags_set_and_check_individual_flags() {
        let mut flags = AtomFlags::empty();
        flags.insert(AtomFlags::IS_FIXED_ROLE);
        assert!(flags.contains(AtomFlags::IS_FIXED_ROLE));
        assert!(!flags.contains(AtomFlags::IS_TREATED_AS_FIXED));
        flags.insert(AtomFlags::IS_TREATED_AS_FIXED);
        assert!(flags.contains(AtomFlags::IS_TREATED_AS_FIXED));
    }

    #[test]
    fn atom_flags_combined_and_removed() {
        let mut flags = AtomFlags::IS_VISIBLE_INTERACTION | AtomFlags::IS_VISIBLE_LATTICE;
        assert!(flags.contains(AtomFlags::IS_VISIBLE_INTERACTION));
        assert!(flags.contains(AtomFlags::IS_VISIBLE_LATTICE));
        flags.remove(AtomFlags::IS_VISIBLE_INTERACTION);
        assert!(!flags.contains(AtomFlags::IS_VISIBLE_INTERACTION));
        assert!(flags.contains(AtomFlags::IS_VISIBLE_LATTICE));
    }

    #[test]
    fn atom_flags_default_is_empty() {
        let flags = AtomFlags::default();
        assert!(!flags.contains(AtomFlags::IS_FIXED_ROLE));
        assert_eq!(flags.bits(), 0);
    }

    #[test]
    fn atom_struct_fields_are_set_correctly() {
        let atom = Atom {
            index: 0,
            serial: 42,
            name: "CA".to_string(),
            res_name: "GLY".to_string(),
            res_id: 5,
            chain_id: 'A',
            force_field_type: "C.3".to_string(),
            partial_charge: -0.123,
            position: Point3::new(1.0, 2.0, 3.0),
            flags: AtomFlags::IS_FIXED_ROLE | AtomFlags::IS_VISIBLE_LATTICE,
            delta: 0.5,
            vdw_radius: 1.7,
            vdw_well_depth: 0.2,
            hbond_type_id: 3,
        };
        assert_eq!(atom.index, 0);
        assert_eq!(atom.serial, 42);
        assert_eq!(atom.name, "CA");
        assert_eq!(atom.res_name, "GLY");
        assert_eq!(atom.res_id, 5);
        assert_eq!(atom.chain_id, 'A');
        assert_eq!(atom.force_field_type, "C.3");
        assert_eq!(atom.partial_charge, -0.123);
        assert_eq!(atom.position, Point3::new(1.0, 2.0, 3.0));
        assert!(atom.flags.contains(AtomFlags::IS_FIXED_ROLE));
        assert!(atom.flags.contains(AtomFlags::IS_VISIBLE_LATTICE));
        assert!(!atom.flags.contains(AtomFlags::IS_TREATED_AS_FIXED));
        assert_eq!(atom.delta, 0.5);
        assert_eq!(atom.vdw_radius, 1.7);
        assert_eq!(atom.vdw_well_depth, 0.2);
        assert_eq!(atom.hbond_type_id, 3);
    }

    #[test]
    fn atom_struct_equality_and_clone() {
        let atom1 = Atom {
            index: 1,
            serial: 2,
            name: "N".to_string(),
            res_name: "ALA".to_string(),
            res_id: 10,
            chain_id: 'B',
            force_field_type: "N.3".to_string(),
            partial_charge: 0.0,
            position: Point3::new(0.0, 0.0, 0.0),
            flags: AtomFlags::empty(),
            delta: 0.0,
            vdw_radius: 1.0,
            vdw_well_depth: 0.1,
            hbond_type_id: 0,
        };
        let atom2 = atom1.clone();
        assert_eq!(atom1, atom2);
    }
}
