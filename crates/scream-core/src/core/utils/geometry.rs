use nalgebra::{Point3, Rotation3, Unit, Vector3};
use std::collections::HashMap;
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy)]
pub struct CbCreationParams {
    pub off_bisector_angle: f64, // Angle offset from the C-N-CA bisector (in degrees)
    pub off_plane_angle: f64,    // Angle offset from the C-N-CA plane (in degrees)
    pub bond_length: f64,        // The desired CA-CB bond length
}
