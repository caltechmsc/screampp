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

pub fn calculate_hn_position(
    n_pos: &Point3<f64>,
    ca_pos: &Point3<f64>,
    prev_c_pos: &Point3<f64>,
    bond_length: f64,
) -> Point3<f64> {
    let n_ca = (ca_pos - n_pos).normalize();
    let n_c_prev = (prev_c_pos - n_pos).normalize();

    let hn_dir = -(n_ca + n_c_prev).normalize();

    n_pos + hn_dir * bond_length
}

pub fn generate_sp3_hydrogens(
    base_pos: &Point3<f64>,
    neighbors: &[Point3<f64>],
    bond_length: f64,
) -> Vec<Point3<f64>> {
    let neighbor_vecs: Vec<Vector3<f64>> = neighbors
        .iter()
        .map(|p| (p - base_pos).normalize())
        .collect();

    match neighbor_vecs.len() {
        1 => {
            let n1 = neighbor_vecs[0];
            let mut temp_vec = if n1.x.abs() < 0.9 {
                Vector3::x()
            } else {
                Vector3::y()
            };
            temp_vec = (temp_vec - n1 * n1.dot(&temp_vec)).normalize();

            let rot_axis = Unit::new_normalize(n1);
            let rot = Rotation3::from_axis_angle(&rot_axis, 120.0f64.to_radians());

            let h1_dir = Rotation3::from_axis_angle(
                &Unit::new_normalize(n1.cross(&temp_vec)),
                109.5f64.to_radians(),
            ) * n1;
            let h1 = base_pos + h1_dir.normalize() * bond_length;
            let h2 = base_pos + (rot * (h1 - base_pos));
            let h3 = base_pos + (rot * (h2 - base_pos));
            vec![h1, h2, h3]
        }
        2 => {
            let n1 = neighbor_vecs[0];
            let n2 = neighbor_vecs[1];
            let bisector = (n1 + n2).normalize();
            let cross_prod = n1.cross(&n2).normalize();

            let h_dir_base = cross_prod.cross(&bisector).normalize();
            let rot = Rotation3::from_axis_angle(
                &Unit::new_normalize(bisector),
                109.5f64.to_radians() / 2.0,
            );

            let h1 = base_pos + (rot * h_dir_base) * bond_length;
            let h2 = base_pos
                + (Rotation3::from_axis_angle(
                    &Unit::new_normalize(bisector),
                    -109.5f64.to_radians(),
                ) * (h1 - base_pos));
            vec![h1, h2]
        }
        3 => {
            let n1 = neighbor_vecs[0];
            let n2 = neighbor_vecs[1];
            let n3 = neighbor_vecs[2];
            let h_dir = -(n1 + n2 + n3).normalize();
            vec![base_pos + h_dir * bond_length]
        }
        _ => panic!("generate_sp3_hydrogens expects 1, 2, or 3 neighbors."),
    }
}

pub fn calculate_rmsd(coords1: &[Point3<f64>], coords2: &[Point3<f64>]) -> Option<f64> {
    if coords1.len() != coords2.len() || coords1.is_empty() {
        return None;
    }
    let n = coords1.len() as f64;
    let squared_dist_sum: f64 = coords1
        .iter()
        .zip(coords2.iter())
        .map(|(p1, p2)| (p1 - p2).norm_squared())
        .sum();
    Some((squared_dist_sum / n).sqrt())
}

pub fn calculate_named_rmsd(
    coords1: &HashMap<String, Point3<f64>>,
    coords2: &HashMap<String, Point3<f64>>,
) -> Option<f64> {
    let mut squared_dist_sum = 0.0;
    let mut count = 0;

    for (name, p1) in coords1 {
        if let Some(p2) = coords2.get(name) {
            squared_dist_sum += (p1 - p2).norm_squared();
            count += 1;
        }
    }
    if count == 0 {
        None
    } else {
        Some((squared_dist_sum / count as f64).sqrt())
    }
}

pub fn find_max_atom_deviation<'a>(
    coords1: &'a HashMap<String, Point3<f64>>,
    coords2: &'a HashMap<String, Point3<f64>>,
) -> Option<(f64, &'a str)> {
    coords1
        .iter()
        .filter_map(|(name, p1)| {
            coords2.get(name).map(|p2| {
                let dist = (p1 - p2).norm();
                (dist, name.as_str())
            })
        })
        .max_by(|(dist1, _), (dist2, _)| {
            dist1
                .partial_cmp(dist2)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
}
