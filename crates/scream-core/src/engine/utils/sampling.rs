use rand::{distributions::WeightedIndex, prelude::*};
use thiserror::Error;
use tracing::instrument;

#[derive(Debug, Error)]
pub enum SamplingError {
    #[error("Input energies list is empty, cannot perform sampling")]
    EmptyEnergies,
    #[error(
        "All energies are too high or beta is zero, resulting in zero total weight for sampling"
    )]
    ZeroTotalWeight,
    #[error("Invalid beta value: {0}. Beta must be positive for Boltzmann sampling")]
    InvalidBeta(f64),
    #[error("Failed to create weighted distribution: {source}")]
    DistributionError {
        #[from]
        source: rand::distributions::WeightedError,
    },
}

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

    let min_energy = *energies
        .iter()
        .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap();

    let weights: Vec<f64> = energies
        .iter()
        .map(|&e| (-(e - min_energy) * beta).exp())
        .collect();

    let total_weight: f64 = weights.iter().sum();
    if total_weight <= f64::EPSILON {
        tracing::warn!(
            "Total Boltzmann weight is near zero ({}). This might indicate a very low temperature or large energy differences, leading to numerical underflow. Returning first index as fallback.",
            total_weight
        );
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
