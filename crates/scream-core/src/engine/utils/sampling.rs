use rand::{distributions::WeightedIndex, prelude::*};
use thiserror::Error;
use tracing::instrument;

/// Defines errors that can occur during sampling operations.
///
/// This enum represents various failure modes in the sampling utilities,
/// particularly for Boltzmann-weighted sampling of energy landscapes.
/// Each variant provides specific information about the nature of the error.
#[derive(Debug, Error)]
pub enum SamplingError {
    /// Indicates that the input energies list is empty, preventing any sampling operation.
    #[error("Input energies list is empty, cannot perform sampling")]
    EmptyEnergies,
    /// Occurs when all Boltzmann weights sum to effectively zero, typically due to
    /// very low temperatures (high beta) or extremely large energy differences.
    #[error(
        "All energies are too high or beta is zero, resulting in zero total weight for sampling"
    )]
    ZeroTotalWeight,
    /// Signals that the provided beta value is invalid for Boltzmann sampling.
    ///
    /// Beta must be positive as it represents the inverse temperature parameter.
    #[error("Invalid beta value: {0}. Beta must be positive for Boltzmann sampling")]
    InvalidBeta(f64),
    /// Wraps errors from the underlying random sampling distribution creation.
    ///
    /// This typically occurs when weights contain invalid values (NaN, negative, etc.).
    #[error("Failed to create weighted distribution: {source}")]
    DistributionError {
        #[from]
        source: rand::distributions::WeightedError,
    },
}

/// Performs Boltzmann-weighted sampling from a list of energies.
///
/// This function implements the Boltzmann distribution sampling algorithm commonly used
/// in molecular simulations and optimization. It converts energy values to probabilities
/// using the formula P(i) ∝ exp(-β(E_i - E_min)), where β is the inverse temperature
/// parameter and E_min is the minimum energy in the set.
///
/// The algorithm handles numerical underflow by detecting when total weights become
/// negligibly small and falling back to selecting the minimum energy state.
///
/// # Arguments
///
/// * `energies` - A slice of energy values to sample from. Must not be empty.
/// * `beta` - The inverse temperature parameter (β = 1/kT). Must be positive.
/// * `rng` - A mutable random number generator implementing the `Rng` trait.
///
/// # Return
///
/// Returns the index of the selected energy state, or a `SamplingError` if sampling fails.
///
/// # Errors
///
/// Returns `SamplingError::EmptyEnergies` if the energies slice is empty.
/// Returns `SamplingError::InvalidBeta` if beta is not positive.
/// Returns `SamplingError::ZeroTotalWeight` if all weights underflow to zero.
/// Returns `SamplingError::DistributionError` if weight distribution creation fails.
#[instrument(level = "trace", skip_all, fields(beta))]
pub fn boltzmann_sample(
    energies: &[f64],
    beta: f64,
    rng: &mut impl Rng,
) -> Result<usize, SamplingError> {
    if energies.is_empty() {
        return Err(SamplingError::EmptyEnergies);
    }
    if beta <= 0.0 {
        return Err(SamplingError::InvalidBeta(beta));
    }

    // Find the minimum energy to improve numerical stability
    let min_energy = *energies
        .iter()
        .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap();

    // Compute Boltzmann weights: exp(-β(E_i - E_min))
    let weights: Vec<f64> = energies
        .iter()
        .map(|&e| (-(e - min_energy) * beta).exp())
        .collect();

    let total_weight: f64 = weights.iter().sum();
    // Handle numerical underflow by checking if total weight is effectively zero
    if total_weight <= f64::EPSILON {
        tracing::warn!(
            "Total Boltzmann weight is near zero ({}). This might indicate a very low temperature or large energy differences, leading to numerical underflow. Returning first index as fallback.",
            total_weight
        );
        // Fallback: return the index of the minimum energy state
        if let Some(idx) = energies
            .iter()
            .position(|&e| (e - min_energy).abs() < f64::EPSILON)
        {
            return Ok(idx);
        } else {
            return Err(SamplingError::ZeroTotalWeight);
        }
    }

    let dist = WeightedIndex::new(&weights)?;
    let selected_index = dist.sample(rng);

    Ok(selected_index)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn boltzmann_sample_with_empty_energies_returns_error() {
        let energies = [];
        let beta = 1.0;
        let mut rng = StdRng::seed_from_u64(42);
        let result = boltzmann_sample(&energies, beta, &mut rng);
        assert!(matches!(result, Err(SamplingError::EmptyEnergies)));
    }

    #[test]
    fn boltzmann_sample_with_zero_beta_returns_error() {
        let energies = [1.0, 2.0];
        let beta = 0.0;
        let mut rng = StdRng::seed_from_u64(42);
        let result = boltzmann_sample(&energies, beta, &mut rng);
        assert!(matches!(result, Err(SamplingError::InvalidBeta(b)) if b == 0.0));
    }

    #[test]
    fn boltzmann_sample_with_negative_beta_returns_error() {
        let energies = [1.0, 2.0];
        let beta = -1.0;
        let mut rng = StdRng::seed_from_u64(42);
        let result = boltzmann_sample(&energies, beta, &mut rng);
        assert!(matches!(result, Err(SamplingError::InvalidBeta(b)) if b == -1.0));
    }

    #[test]
    fn boltzmann_sample_selects_lowest_energy_at_very_low_temperature() {
        let energies = [10.0, 1.0, 20.0, 100.0];
        let beta = 1000.0;
        let mut rng = StdRng::seed_from_u64(42);
        let result = boltzmann_sample(&energies, beta, &mut rng).unwrap();
        assert_eq!(result, 1);
    }

    #[test]
    fn boltzmann_sample_is_within_bounds_for_valid_input() {
        let energies = [10.0, 1.0, 20.0, 100.0];
        let beta = 1.0;
        let mut rng = StdRng::seed_from_u64(42);
        let result = boltzmann_sample(&energies, beta, &mut rng).unwrap();
        assert!(result < energies.len());
    }

    #[test]
    fn boltzmann_sample_handles_numerical_underflow_by_returning_min_energy_index() {
        let energies = [0.0, 1000.0, 2000.0];
        let beta = 1000.0;
        let mut rng = StdRng::seed_from_u64(42);
        let result = boltzmann_sample(&energies, beta, &mut rng).unwrap();
        assert_eq!(result, 0);
    }

    #[test]
    fn boltzmann_sample_with_all_equal_energies_is_uniform() {
        let energies = [5.0, 5.0, 5.0, 5.0];
        let beta = 1.0;
        let mut rng = StdRng::seed_from_u64(42);
        let mut counts = vec![0; energies.len()];
        let n_samples = 10000;

        for _ in 0..n_samples {
            let index = boltzmann_sample(&energies, beta, &mut rng).unwrap();
            counts[index] += 1;
        }

        let expected_count = (n_samples / energies.len()) as f64;
        for count in counts {
            let deviation = (count as f64 - expected_count).abs() / expected_count;
            assert!(deviation < 0.1, "Distribution is not uniform enough");
        }
    }

    #[test]
    fn boltzmann_sample_with_single_energy_returns_zero() {
        let energies = [123.45];
        let beta = 1.0;
        let mut rng = StdRng::seed_from_u64(42);
        let result = boltzmann_sample(&energies, beta, &mut rng).unwrap();
        assert_eq!(result, 0);
    }
}
