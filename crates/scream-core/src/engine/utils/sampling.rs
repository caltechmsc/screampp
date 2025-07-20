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
