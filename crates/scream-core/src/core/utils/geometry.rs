use nalgebra::{Point3, Rotation3, Unit, Vector3};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
pub struct CbCreationParams {
    pub off_bisector_angle: f64,
    pub off_plane_angle: f64,
    pub bond_length: f64,
}

pub fn bond_angle(p1: &Point3<f64>, p2_vertex: &Point3<f64>, p3: &Point3<f64>) -> f64 {
    let v1 = p1 - p2_vertex;
    let v2 = p3 - p2_vertex;
    v1.angle(&v2).to_degrees()
}

pub fn dihedral_angle(
    p1: &Point3<f64>,
    p2: &Point3<f64>,
    p3: &Point3<f64>,
    p4: &Point3<f64>,
) -> f64 {
    let b1 = p1 - p2;
    let b2 = p3 - p2;
    let b3 = p4 - p3;

    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);

    let angle_rad = n1.angle(&n2);

    if n1.cross(&n2).dot(&b2) < 0.0 {
        -angle_rad.to_degrees()
    } else {
        angle_rad.to_degrees()
    }
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
) -> Result<Vec<Point3<f64>>, &'static str> {
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
            Ok(vec![h1, h2, h3])
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
            Ok(vec![h1, h2])
        }
        3 => {
            let n1 = neighbor_vecs[0];
            let n2 = neighbor_vecs[1];
            let n3 = neighbor_vecs[2];
            let h_dir = -(n1 + n2 + n3).normalize();
            Ok(vec![base_pos + h_dir * bond_length])
        }
        _ => Err("generate_sp3_hydrogens expects 1, 2, or 3 neighbors."),
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPSILON: f64 = 1e-9;

    fn points_approx_equal(p1: &Point3<f64>, p2: &Point3<f64>) -> bool {
        (p1 - p2).norm() < EPSILON
    }

    fn f64_approx_equal(a: f64, b: f64) -> bool {
        (a - b).abs() < EPSILON
    }

    #[test]
    fn test_bond_angle() {
        let p1 = Point3::new(1.0, 0.0, 0.0);
        let p2 = Point3::new(0.0, 0.0, 0.0);
        let p3 = Point3::new(0.0, 1.0, 0.0);
        assert!(f64_approx_equal(bond_angle(&p1, &p2, &p3), 90.0));

        let p3_straight = Point3::new(-1.0, 0.0, 0.0);
        assert!(f64_approx_equal(bond_angle(&p1, &p2, &p3_straight), 180.0));

        let p3_cis = Point3::new(1.0, 1.0, 0.0);
        assert!(f64_approx_equal(bond_angle(&p1, &p2, &p3_cis), 45.0));
    }

    #[test]
    fn test_dihedral_angle() {
        let p1 = Point3::new(0.0, -1.0, 0.0);
        let p2 = Point3::new(0.0, 0.0, 0.0);
        let p3 = Point3::new(1.0, 0.0, 0.0);
        let p4 = Point3::new(1.0, 0.0, 1.0);
        assert!(f64_approx_equal(dihedral_angle(&p1, &p2, &p3, &p4), 90.0));

        let p4_neg = Point3::new(1.0, 0.0, -1.0);
        assert!(f64_approx_equal(
            dihedral_angle(&p1, &p2, &p3, &p4_neg),
            -90.0
        ));

        let p4_trans = Point3::new(1.0, -1.0, 0.0);
        assert!(f64_approx_equal(
            dihedral_angle(&p1, &p2, &p3, &p4_trans).abs(),
            180.0
        ));

        let p4_cis = Point3::new(1.0, 1.0, 0.0);
        assert!(f64_approx_equal(
            dihedral_angle(&p1, &p2, &p3, &p4_cis),
            0.0
        ));
    }

    #[test]
    fn test_rotation_to_align() {
        let from = Vector3::x();
        let to = Vector3::y();
        let rot = rotation_to_align(&from, &to).unwrap();

        let rotated = rot * from;
        assert!((rotated - to).norm() < EPSILON);

        let expected_rot = Rotation3::from_axis_angle(&Vector3::z_axis(), PI / 2.0);
        assert!(
            (rot.axis_angle().unwrap().1 - expected_rot.axis_angle().unwrap().1).abs() < EPSILON
        );

        let rot_identity = rotation_to_align(&from, &from).unwrap();
        assert!((rot_identity * from - from).norm() < EPSILON);

        assert!(rotation_to_align(&from, &-from).is_none());
    }

    #[test]
    fn test_rotation_from_axis_angle() {
        let axis = Vector3::new(0.0, 0.0, 2.0);
        let angle_degrees = 90.0;
        let rot = rotation_from_axis_angle(&axis, angle_degrees);

        let vec_to_rotate = Vector3::new(1.0, 0.0, 0.0);
        let rotated_vec = rot * vec_to_rotate;

        let expected_vec = Vector3::new(0.0, 1.0, 0.0);
        assert!((rotated_vec - expected_vec).norm() < EPSILON);
    }

    #[test]
    fn test_calculate_cb_position() {
        let ca_pos = Point3::new(0.0, 0.0, 0.0);
        let n_pos = Point3::new(-1.0, 0.0, 0.0);
        let c_pos = Point3::new(0.0, -1.0, 0.0);

        let params = CbCreationParams {
            off_bisector_angle: 0.0,
            off_plane_angle: 0.0,
            bond_length: 1.5,
        };

        let cb_pos = calculate_cb_position(&n_pos, &ca_pos, &c_pos, &params);
        let expected_dir = (Vector3::new(1.0, 1.0, 0.0)).normalize();
        let expected_pos = ca_pos + expected_dir * 1.5;

        assert!(points_approx_equal(&cb_pos, &expected_pos));
        assert!(f64_approx_equal((cb_pos - ca_pos).norm(), 1.5));
    }

    #[test]
    fn test_calculate_hn_position() {
        let n_pos = Point3::new(0.0, 0.0, 0.0);
        let ca_pos = Point3::new(1.0, 0.0, 0.0);
        let prev_c_pos = Point3::new(0.0, 1.0, 0.0);
        let bond_length = 1.0;

        let hn_pos = calculate_hn_position(&n_pos, &ca_pos, &prev_c_pos, bond_length);

        let expected_dir = (Vector3::new(-1.0, -1.0, 0.0)).normalize();
        let expected_pos = n_pos + expected_dir * bond_length;

        assert!(points_approx_equal(&hn_pos, &expected_pos));
        assert!(f64_approx_equal((hn_pos - n_pos).norm(), bond_length));
    }

    #[test]
    fn generate_sp3_hydrogens_one_neighbor_returns_three_hydrogens() {
        let base_pos = Point3::new(0.0, 0.0, 0.0);
        let neighbors = [Point3::new(1.0, 0.0, 0.0)];
        let bond_length = 1.0;

        let hydrogens = generate_sp3_hydrogens(&base_pos, &neighbors, bond_length).unwrap();
        assert_eq!(hydrogens.len(), 3);

        let neighbor_vec = (neighbors[0] - base_pos).normalize();
        let expected_angle_rad = 109.5f64.to_radians();

        for h_pos in &hydrogens {
            assert!(f64_approx_equal((h_pos - base_pos).norm(), bond_length));
            let h_vec = (h_pos - base_pos).normalize();
            assert!(f64_approx_equal(
                h_vec.angle(&neighbor_vec),
                expected_angle_rad
            ));
        }
    }

    #[test]
    fn generate_sp3_hydrogens_two_neighbors_returns_two_hydrogens() {
        let base_pos = Point3::new(0.0, 0.0, 0.0);
        let neighbors = [Point3::new(1.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)];
        let bond_length = 1.0;

        let hydrogens = generate_sp3_hydrogens(&base_pos, &neighbors, bond_length).unwrap();
        assert_eq!(hydrogens.len(), 2);

        assert!(f64_approx_equal(hydrogens[0].z, -hydrogens[1].z));
        assert!(f64_approx_equal(hydrogens[0].x, hydrogens[1].x));
        assert!(f64_approx_equal(hydrogens[0].y, hydrogens[1].y));
    }

    #[test]
    fn generate_sp3_hydrogens_three_neighbors_returns_one_hydrogen() {
        let base_pos = Point3::new(0.0, 0.0, 0.0);
        let neighbors = [
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ];
        let bond_length = 1.0;

        let hydrogens = generate_sp3_hydrogens(&base_pos, &neighbors, bond_length).unwrap();
        assert_eq!(hydrogens.len(), 1);

        let h_pos = hydrogens[0];
        let expected_dir = -(Vector3::x() + Vector3::y() + Vector3::z()).normalize();
        let expected_pos = base_pos + expected_dir * bond_length;

        assert!(points_approx_equal(&h_pos, &expected_pos));
    }

    #[test]
    #[should_panic]
    fn generate_sp3_hydrogens_panics_on_invalid_count() {
        let base_pos = Point3::new(0.0, 0.0, 0.0);
        let neighbors = [];
        generate_sp3_hydrogens(&base_pos, &neighbors, 1.0).unwrap();
    }

    #[test]
    fn test_calculate_rmsd() {
        let coords1 = [Point3::new(0.0, 0.0, 0.0), Point3::new(3.0, 4.0, 0.0)];
        let coords2 = [Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0)];

        let rmsd = calculate_rmsd(&coords1, &coords2).unwrap();
        assert!(f64_approx_equal(rmsd, 12.5_f64.sqrt()));

        assert!(f64_approx_equal(
            calculate_rmsd(&coords1, &coords1).unwrap(),
            0.0
        ));

        let short_coords = [Point3::origin()];
        assert!(calculate_rmsd(&coords1, &short_coords).is_none());

        assert!(calculate_rmsd(&[], &[]).is_none());
    }

    #[test]
    fn test_calculate_named_rmsd() {
        let mut coords1 = HashMap::new();
        coords1.insert("A".to_string(), Point3::new(0.0, 0.0, 0.0));
        coords1.insert("B".to_string(), Point3::new(3.0, 0.0, 0.0));
        coords1.insert("C".to_string(), Point3::new(0.0, 0.0, 0.0));

        let mut coords2 = HashMap::new();
        coords2.insert("A".to_string(), Point3::new(0.0, 0.0, 0.0));
        coords2.insert("B".to_string(), Point3::new(0.0, 4.0, 0.0));
        coords2.insert("D".to_string(), Point3::new(0.0, 0.0, 0.0));

        let rmsd = calculate_named_rmsd(&coords1, &coords2).unwrap();
        assert!(f64_approx_equal(rmsd, 12.5_f64.sqrt()));

        let mut coords3 = HashMap::new();
        coords3.insert("X".to_string(), Point3::origin());
        assert!(calculate_named_rmsd(&coords1, &coords3).is_none());
    }

    #[test]
    fn test_find_max_atom_deviation() {
        let mut coords1 = HashMap::new();
        coords1.insert("A".to_string(), Point3::new(0.0, 0.0, 0.0));
        coords1.insert("B".to_string(), Point3::new(0.0, 0.0, 0.0));
        coords1.insert("C".to_string(), Point3::new(0.0, 0.0, 0.0));

        let mut coords2 = HashMap::new();
        coords2.insert("A".to_string(), Point3::new(1.0, 0.0, 0.0));
        coords2.insert("B".to_string(), Point3::new(3.0, 4.0, 0.0));
        coords2.insert("C".to_string(), Point3::new(0.0, 3.0, 0.0));
        coords2.insert("D".to_string(), Point3::new(9.0, 9.0, 9.0));

        let result = find_max_atom_deviation(&coords1, &coords2);
        assert!(result.is_some());
        let (max_dev, atom_name) = result.unwrap();

        assert_eq!(atom_name, "B");
        assert!(f64_approx_equal(max_dev, 5.0));

        let coords3 = HashMap::new();
        assert!(find_max_atom_deviation(&coords1, &coords3).is_none());
    }
}
