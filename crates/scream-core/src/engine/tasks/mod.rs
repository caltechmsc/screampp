//! Tasks for side-chain placement optimization and energy calculations.
//!
//! Tasks are the core computational units that perform specific calculations or optimizations
//! during protein side-chain placement. Each submodule implements a different type of task,
//! such as energy calculations, optimization algorithms, and clash detection. Tasks are
//! designed to be modular and composable, allowing for flexible optimization workflows.

pub mod clash_detection;
pub mod doublet_optimization;
pub mod el_energy;
pub mod fixed_energy;
pub mod interaction_energy;
