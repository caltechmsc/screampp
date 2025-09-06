/// The Coulomb constant used in electrostatic potential calculations.
///
/// This constant represents the value of 1/(4πε₀) in units of kcal·Å/(mol·e²),
/// which is the standard unit system used in molecular mechanics force fields.
const COULOMB_CONSTANT: f64 = 332.0637; // In kcal·Å/(mol·e²)

/// Calculates the Lennard-Jones 12-6 potential energy between two atoms.
///
/// This function implements the classic Lennard-Jones potential, which models
/// van der Waals interactions with a repulsive r⁻¹² term and an attractive r⁻⁶ term.
/// The potential reaches its minimum at the specified `r_min` distance.
///
/// # Arguments
///
/// * `dist` - The distance between the two atoms.
/// * `r_min` - The distance at which the potential reaches its minimum.
/// * `well_depth` - The depth of the potential well (negative value).
///
/// # Return
///
/// Returns the potential energy. Positive values indicate repulsion, negative values indicate attraction.
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

/// Calculates the Buckingham exponential-6 potential energy between two atoms.
///
/// This function implements a modified Buckingham potential that switches to
/// pure Lennard-Jones repulsion at short distances to prevent numerical instability.
/// The potential combines exponential repulsion with r⁻⁶ dispersion attraction.
///
/// # Arguments
///
/// * `dist` - The distance between the two atoms.
/// * `r_min` - The distance parameter for the potential.
/// * `well_depth` - The depth of the potential well.
/// * `gamma` - The exponential decay parameter.
///
/// # Return
///
/// Returns the potential energy. The function switches to r⁻¹² repulsion below 60% of `r_min`.
#[inline]
pub fn buckingham_exp_6(dist: f64, r_min: f64, well_depth: f64, gamma: f64) -> f64 {
    const POTENTIAL_SWITCHING_FACTOR: f64 = 0.6;
    let switching_distance = POTENTIAL_SWITCHING_FACTOR * r_min;

    if dist < switching_distance {
        let rho = r_min / dist;
        return well_depth * rho.powi(12);
    }

    if dist < 1e-6 {
        return 1e10;
    }

    let rho = dist / r_min;

    let factor = gamma / (gamma - 6.0);
    well_depth * (6.0 / (gamma - 6.0) * (gamma * (1.0 - rho)).exp() - factor * rho.powi(-6))
}

/// Calculates the Coulomb electrostatic potential energy between two charged atoms.
///
/// This function computes the electrostatic interaction energy using Coulomb's law
/// with the appropriate constant for molecular mechanics units.
///
/// # Arguments
///
/// * `dist` - The distance between the two atoms.
/// * `q1` - The charge of the first atom.
/// * `q2` - The charge of the second atom.
/// * `dielectric` - The dielectric constant of the medium.
///
/// # Return
///
/// Returns the electrostatic potential energy. The sign depends on the charges.
#[inline]
pub fn coulomb(dist: f64, q1: f64, q2: f64, dielectric: f64) -> f64 {
    if dist < 1e-6 {
        return q1.signum() * q2.signum() * 1e10;
    }
    COULOMB_CONSTANT * q1 * q2 / (dielectric * dist)
}

/// Calculates the Dreiding hydrogen bond 12-10 potential energy.
///
/// This function implements the specialized hydrogen bond potential used in the
/// Dreiding force field, which combines r⁻¹² and r⁻¹⁰ terms for hydrogen bonding
/// interactions between donor and acceptor atoms.
///
/// # Arguments
///
/// * `dist_ad` - The distance between the donor and acceptor atoms.
/// * `r_hb` - The equilibrium hydrogen bond distance.
/// * `d_hb` - The hydrogen bond well depth.
///
/// # Return
///
/// Returns the hydrogen bond potential energy.
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

/// Applies a flat-bottom modification to a van der Waals potential function.
///
/// This function modifies the behavior of a potential in the repulsive region by
/// creating a flat energy well around the ideal distance, which can improve
/// numerical stability in molecular dynamics simulations.
///
/// # Arguments
///
/// * `dist` - The actual distance between atoms.
/// * `ideal_dist` - The ideal equilibrium distance.
/// * `delta` - The width of the flat-bottom region.
/// * `potential_fn` - The base potential function to modify.
///
/// # Return
///
/// Returns the modified potential energy.
#[inline]
pub fn apply_flat_bottom_vdw<F>(dist: f64, ideal_dist: f64, delta: f64, potential_fn: F) -> f64
where
    F: Fn(f64) -> f64,
{
    const REPULSIVE_CORE_SCALE_FACTOR: f64 = 0.8;
    let repulsive_core_boundary = REPULSIVE_CORE_SCALE_FACTOR * ideal_dist;

    if dist < repulsive_core_boundary {
        return potential_fn(dist);
    }

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

/// Applies a flat-bottom modification to a hydrogen bond potential function.
///
/// This function creates a flat energy region around the ideal hydrogen bond
/// distance to stabilize the interaction while maintaining the correct asymptotic
/// behavior at long and short ranges.
///
/// # Arguments
///
/// * `dist` - The actual distance between atoms.
/// * `ideal_dist` - The ideal hydrogen bond distance.
/// * `delta` - The width of the flat-bottom region.
/// * `potential_fn` - The base potential function to modify.
///
/// # Return
///
/// Returns the modified potential energy.
#[inline]
pub fn apply_flat_bottom_hbond<F>(dist: f64, ideal_dist: f64, delta: f64, potential_fn: F) -> f64
where
    F: Fn(f64) -> f64,
{
    const REPULSIVE_CORE_SCALE_FACTOR: f64 = 0.75;
    let repulsive_core_boundary = REPULSIVE_CORE_SCALE_FACTOR * ideal_dist;

    if dist < repulsive_core_boundary {
        return potential_fn(dist);
    }

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
    fn buckingham_at_minimum_distance_returns_negative_well_depth() {
        let energy = buckingham_exp_6(2.0, 2.0, 10.0, 12.0);
        assert!(f64_approx_equal(energy, -10.0));
    }

    #[test]
    fn buckingham_is_repulsive_in_the_safe_zone_but_below_minimum() {
        let energy = buckingham_exp_6(1.3, 2.0, 10.0, 12.0);
        assert!(energy > 0.0);

        let lj_repulsion_energy = 10.0 * (2.0_f64 / 1.3).powi(12);
        assert!(
            (energy - lj_repulsion_energy).abs() > 1e-3,
            "Should not be using LJ repulsion here"
        );
    }

    #[test]
    fn buckingham_switches_to_lj_repulsion_to_prevent_catastrophe() {
        let r_min = 2.0;
        let well_depth = 10.0;
        let dist = 0.1;

        let energy = buckingham_exp_6(dist, r_min, well_depth, 12.0);

        let expected_lj_repulsion = well_depth * (r_min / dist).powi(12);

        assert!(
            energy > 1e10,
            "Energy at very short distance should be a huge positive number"
        );
        assert!(
            f64_approx_equal(energy, expected_lj_repulsion),
            "Must switch to LJ r^-12 potential to prevent collapse"
        );
    }

    #[test]
    fn apply_flat_bottom_vdw_is_unchanged_beyond_ideal_distance() {
        let potential = apply_flat_bottom_vdw(10.0, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential, 100.0));
    }

    #[test]
    fn apply_flat_bottom_vdw_is_flat_in_well_region() {
        let potential = apply_flat_bottom_vdw(7.5, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential, 64.0));
    }

    #[test]
    fn apply_flat_bottom_vdw_uses_shifted_distance_in_soft_repulsion_zone() {
        let potential = apply_flat_bottom_vdw(6.8, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential, 7.8 * 7.8));
    }

    #[test]
    fn apply_flat_bottom_vdw_with_zero_delta_is_identity() {
        let potential = apply_flat_bottom_vdw(5.0, 8.0, 0.0, |d| d * d);
        assert!(f64_approx_equal(potential, 25.0));
    }

    #[test]
    fn apply_flat_bottom_vdw_bypasses_softening_in_hardcore_repulsion_zone() {
        let potential = apply_flat_bottom_vdw(5.0, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential, 5.0 * 5.0));
        assert!(!f64_approx_equal(potential, 6.0 * 6.0));
    }

    #[test]
    fn apply_flat_bottom_hbond_is_unchanged_in_attractive_tail() {
        let potential = apply_flat_bottom_hbond(10.0, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential, 9.0 * 9.0));
    }

    #[test]
    fn apply_flat_bottom_hbond_is_flat_around_ideal_distance() {
        let potential1 = apply_flat_bottom_hbond(7.5, 8.0, 1.0, |d| d * d);
        let potential2 = apply_flat_bottom_hbond(8.5, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential1, 64.0));
        assert!(f64_approx_equal(potential2, 64.0));
    }

    #[test]
    fn apply_flat_bottom_hbond_uses_shifted_distance_in_soft_repulsion_zone() {
        let potential = apply_flat_bottom_hbond(6.8, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential, 7.8 * 7.8));
    }

    #[test]
    fn apply_flat_bottom_hbond_with_zero_delta_is_identity() {
        let potential = apply_flat_bottom_hbond(5.0, 8.0, 0.0, |d| d * d);
        assert!(f64_approx_equal(potential, 25.0));
    }

    #[test]
    fn apply_flat_bottom_hbond_bypasses_softening_in_hardcore_repulsion_zone() {
        let potential = apply_flat_bottom_hbond(5.0, 8.0, 1.0, |d| d * d);
        assert!(f64_approx_equal(potential, 5.0 * 5.0));
        assert!(!f64_approx_equal(potential, 6.0 * 6.0));
    }
}
