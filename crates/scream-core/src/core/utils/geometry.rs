use nalgebra::{Point3, Rotation3, Unit, Vector3};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
pub struct CbCreationParams {
    pub off_bisector_angle: f64,
    pub off_plane_angle: f64,
    pub bond_length: f64,
}

pub fn rotation_to_align(from: &Vector3<f64>, to: &Vector3<f64>) -> Option<Rotation3<f64>> {
    Rotation3::rotation_between(from, to)
}
