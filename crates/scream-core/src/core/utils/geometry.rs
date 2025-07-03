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

pub fn rotation_from_axis_angle(axis: &Vector3<f64>, angle_degrees: f64) -> Rotation3<f64> {
    Rotation3::from_axis_angle(&Unit::new_normalize(*axis), angle_degrees.to_radians())
}

pub fn calculate_cb_position(
    n_pos: &Point3<f64>,
    ca_pos: &Point3<f64>,
    c_pos: &Point3<f64>,
    params: &CbCreationParams,
) -> Point3<f64> {
    let ca_n = (n_pos - ca_pos).normalize();
    let ca_c = (c_pos - ca_pos).normalize();

    let bisector = -(ca_n + ca_c).normalize();
    let plane_normal = Unit::new_normalize(ca_n.cross(&ca_c));

    let rot_off_bisector =
        Rotation3::from_axis_angle(&plane_normal, params.off_bisector_angle.to_radians());
    let cb_vec_in_plane = rot_off_bisector * bisector;

    let in_plane_axis = Unit::new_normalize(cb_vec_in_plane);
    let rot_off_plane =
        Rotation3::from_axis_angle(&in_plane_axis, params.off_plane_angle.to_radians());
    let final_cb_vec = rot_off_plane * cb_vec_in_plane;

    ca_pos + final_cb_vec * params.bond_length
}
