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

#[cfg(test)]
mod tests {
    use super::*;

    const TOLERANCE: f64 = 1e-9;

    fn f64_approx_equal(a: f64, b: f64) -> bool {
        (a - b).abs() < TOLERANCE
    }

    #[test]
    fn lennard_jones_at_minimum_distance_returns_negative_well_depth() {
        let energy = lennard_jones_12_6(2.0, 2.0, 10.0);
        assert!(f64_approx_equal(energy, -10.0));
    }

    #[test]
    fn lennard_jones_at_very_small_distance_returns_large_positive_energy() {
        let energy = lennard_jones_12_6(1e-7, 2.0, 10.0);
        assert!(f64_approx_equal(energy, 1e10));
    }

    #[test]
    fn buckingham_at_minimum_distance_returns_negative_well_depth() {
        let energy = buckingham_exp_6(2.0, 2.0, 10.0, 12.0);
        assert!(f64_approx_equal(energy, -10.0));
    }

    #[test]
    fn buckingham_at_very_small_distance_returns_large_positive_energy() {
        let energy = buckingham_exp_6(1e-7, 2.0, 10.0, 12.0);
        assert!(f64_approx_equal(energy, 1e10));
    }

    #[test]
    fn coulomb_calculates_repulsive_force_correctly() {
        let energy = coulomb(1.0, 1.0, 1.0, 1.0);
        assert!(f64_approx_equal(energy, COULOMB_CONSTANT));
    }

    #[test]
    fn coulomb_calculates_attractive_force_correctly() {
        let energy = coulomb(2.0, 1.0, -1.0, 1.0);
        assert!(f64_approx_equal(energy, -COULOMB_CONSTANT / 2.0));
    }

    #[test]
    fn coulomb_at_very_small_distance_returns_large_energy_with_correct_sign() {
        assert!(f64_approx_equal(coulomb(1e-7, 1.0, 1.0, 1.0), 1e10));
        assert!(f64_approx_equal(coulomb(1e-7, -1.0, 1.0, 1.0), -1e10));
    }

    #[test]
    fn dreiding_hbond_12_10_at_equilibrium_distance_returns_negative_well_depth() {
        let energy = dreiding_hbond_12_10(2.7, 2.7, 5.0);
        assert!(f64_approx_equal(energy, -5.0));
    }

    #[test]
    fn dreiding_hbond_12_10_at_very_small_distance_returns_large_positive_energy() {
        let energy = dreiding_hbond_12_10(1e-7, 2.7, 5.0);
        assert!(f64_approx_equal(energy, 1e10));
    }

    #[test]
    fn apply_flat_bottom_vdw_uses_base_potential_beyond_ideal_distance() {
        let potential = apply_flat_bottom_vdw(10.0, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential, 100.0));
    }

    #[test]
    fn apply_flat_bottom_vdw_is_flat_within_delta_of_ideal_distance() {
        let potential = apply_flat_bottom_vdw(7.5, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential, 64.0));
    }

    #[test]
    fn apply_flat_bottom_vdw_uses_shifted_distance_below_flat_region() {
        let potential = apply_flat_bottom_vdw(6.0, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential, 49.0));
    }

    #[test]
    fn apply_flat_bottom_vdw_with_zero_delta_is_identity() {
        let potential = apply_flat_bottom_vdw(5.0, 8.0, 0.0, |d| d * d);
        assert!(f64_approx_equal(potential, 25.0));
    }

    #[test]
    fn apply_flat_bottom_hbond_is_flat_around_ideal_distance() {
        let potential = apply_flat_bottom_hbond(8.5, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential, 64.0));
    }

    #[test]
    fn apply_flat_bottom_hbond_uses_shifted_distance_above_flat_region() {
        let potential = apply_flat_bottom_hbond(10.0, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential, 81.0));
    }

    #[test]
    fn apply_flat_bottom_hbond_uses_shifted_distance_below_flat_region() {
        let potential = apply_flat_bottom_hbond(6.0, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential, 49.0));
    }

    #[test]
    fn apply_flat_bottom_hbond_with_zero_delta_is_identity() {
        let potential = apply_flat_bottom_hbond(5.0, 8.0, 0.0, |d| d * d);
        assert!(f64_approx_equal(potential, 25.0));
    }
}
