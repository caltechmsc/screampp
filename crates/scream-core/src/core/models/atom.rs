use super::ids::ResidueId;
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
    pub serial: usize,         // Atom serial number from source file
    pub name: String,          // Atom name (e.g., "CA", "N", "O")
    pub residue_id: ResidueId, // Stable ID of the parent residue

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

impl Atom {
    pub fn new(serial: usize, name: &str, residue_id: ResidueId, position: Point3<f64>) -> Self {
        Self {
            serial,
            name: name.to_string(),
            residue_id,
            position,
            force_field_type: "".to_string(),
            partial_charge: 0.0,
            flags: AtomFlags::default(),
            delta: 0.0,
            vdw_radius: 0.0,
            vdw_well_depth: 0.0,
            hbond_type_id: -1,
        }
    }
}
