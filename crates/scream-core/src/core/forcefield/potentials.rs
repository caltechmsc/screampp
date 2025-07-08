use std::f64::consts::PI;

const COULOMB_CONSTANT: f64 = 332.0637; // In kcal·Å/(mol·e²)

#[inline]
pub fn lennard_jones_12_6(dist: f64, r_min: f64, well_depth: f64) -> f64 {
    if dist < 1e-6 {
        return 1e10;
    }
    let rho = r_min / dist;
    let rho6 = rho.powi(6);
    let rho12 = rho6 * rho6;
    well_depth * (rho12 - 2.0 * rho6)
}

#[inline]
pub fn buckingham_exp_6(dist: f64, r_min: f64, well_depth: f64, gamma: f64) -> f64 {
    if dist < 1e-6 {
        return 1e10;
    }
    let rho = dist / r_min;
    if rho < 0.1 {
        return 1e10;
    }

    let factor = gamma / (gamma - 6.0);
    well_depth * (6.0 / (gamma - 6.0) * (gamma * (1.0 - rho)).exp() - factor * rho.powi(-6))
}

#[inline]
pub fn coulomb(dist: f64, q1: f64, q2: f64, dielectric: f64) -> f64 {
    if dist < 1e-6 {
        return q1.signum() * q2.signum() * 1e10;
    }
    COULOMB_CONSTANT * q1 * q2 / (dielectric * dist)
}

#[inline]
pub fn dreiding_hbond_12_10(dist_ad: f64, r_hb: f64, d_hb: f64) -> f64 {
    if dist_ad < 1e-6 {
        return 1e10;
    }
    let rho = r_hb / dist_ad;
    let rho10 = rho.powi(10);
    let rho12 = rho.powi(12);
    d_hb * (5.0 * rho12 - 6.0 * rho10)
}

#[inline]
pub fn dreiding_hbond(dist_ad: f64, angle_ahd_deg: f64, r_hb: f64, d_hb: f64) -> f64 {
    if angle_ahd_deg < 90.0 {
        return 0.0;
    }
    let cos_theta = (angle_ahd_deg * PI / 180.0).cos();
    let angular_term = cos_theta.powi(4);
    let distance_term = dreiding_hbond_12_10(dist_ad, r_hb, d_hb);
    distance_term * angular_term
}

#[inline]
pub fn apply_flat_bottom_vdw<F>(dist: f64, ideal_dist: f64, delta: f64, potential_fn: F) -> f64
where
    F: Fn(f64) -> f64,
{
    if delta <= 1e-9 {
        return potential_fn(dist);
    }

    if dist >= ideal_dist {
        potential_fn(dist)
    } else if dist > ideal_dist - delta {
        potential_fn(ideal_dist)
    } else {
        let effective_dist = dist + delta;
        potential_fn(effective_dist)
    }
}

#[inline]
pub fn apply_flat_bottom_hbond<F>(dist: f64, ideal_dist: f64, delta: f64, potential_fn: F) -> f64
where
    F: Fn(f64) -> f64,
{
    if delta <= 1e-9 {
        return potential_fn(dist);
    }

    if dist >= ideal_dist + delta {
        let effective_dist = dist - delta;
        potential_fn(effective_dist)
    } else if dist <= ideal_dist - delta {
        let effective_dist = dist + delta;
        potential_fn(effective_dist)
    } else {
        potential_fn(ideal_dist)
    }
}
