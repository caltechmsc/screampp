use std::f64::consts::PI;

const COULOMB_CONSTANT: f64 = 332.0637; // In kcal·Å/(mol·e²)

#[inline]
pub fn lennard_jones_12_6(dist: f64, r_min: f64, well_depth: f64) -> f64 {
    if dist < 1e-6 {
        return 1e10; // Return a large positive energy for clashes
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
        return 1e10; // Prevent numerical issues for very small distances
    }

    let factor = gamma / (gamma - 6.0);
    well_depth * (6.0 / (gamma - 6.0) * (gamma * (1.0 - rho)).exp() - factor * rho.powi(-6))
}

#[inline]
pub fn coulomb(dist: f64, q1: f64, q2: f64, dielectric: f64) -> f64 {
    if dist < 1e-6 {
        return q1.signum() * q2.signum() * 1e10; // Return a large energy of the correct sign
    }
    COULOMB_CONSTANT * q1 * q2 / (dielectric * dist)
}

#[inline]
pub fn dreiding_hbond(dist_ad: f64, angle_ahd_deg: f64, r_hb: f64, d_hb: f64) -> f64 {
    if angle_ahd_deg < 90.0 {
        return 0.0;
    }
    let cos_theta = (angle_ahd_deg * PI / 180.0).cos();

    let angular_factor = cos_theta.powi(4);

    if dist_ad < 1e-6 {
        return 1e10;
    }
    let rho = r_hb / dist_ad;
    let rho10 = rho.powi(10);
    let rho12 = rho.powi(12);
    let distance_term = d_hb * (5.0 * rho12 - 6.0 * rho10);

    distance_term * angular_factor
}

#[inline]
pub fn apply_flat_bottom<F>(dist: f64, ideal_dist: f64, delta: f64, base_potential_fn: F) -> f64
where
    F: Fn(f64) -> f64,
{
    if dist >= ideal_dist {
        return base_potential_fn(dist);
    }

    if dist > ideal_dist - delta {
        return base_potential_fn(ideal_dist);
    }

    let effective_dist = dist + delta;
    base_potential_fn(effective_dist)
}
