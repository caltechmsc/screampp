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
