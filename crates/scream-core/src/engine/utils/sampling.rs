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
