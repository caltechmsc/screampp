//! # Engine Module
//!
//! This module implements the optimization engine for automated protein side-chain placement
//! in SCREAM++, providing the computational framework for conformational optimization and
//! energy minimization workflows.
//!
//! ## Overview
//!
//! The engine module orchestrates the complete optimization process for protein side-chain
//! conformation prediction. It manages optimization state, coordinates computational tasks,
//! handles energy calculations, and provides a flexible framework for different optimization
//! strategies and algorithms.
//!
//! ## Architecture
//!
//! The module is organized into specialized submodules that handle different aspects
//! of the optimization process:
//!
//! - **Configuration** ([`config`]) - Optimization parameters, convergence criteria, and settings
//! - **State Tracking** ([`state`]) - Solution states, optimization progress, and result management
//! - **Progress Monitoring** ([`progress`]) - Progress reporting and user feedback mechanisms
//! - **Error Handling** ([`error`]) - Engine-specific error types and error propagation
//!
//! ## Key Capabilities
//!
//! - **Multi-strategy optimization** supporting various algorithms and convergence criteria
//! - **Parallel computation** for energy calculations and optimization tasks
//! - **Energy caching** to avoid redundant calculations during optimization
//! - **Flexible residue selection** for targeted optimization of specific protein regions
//! - **Progress monitoring** with detailed reporting and convergence tracking
//! - **Transactional modifications** ensuring system consistency during optimization
//! - **Extensible task system** allowing custom optimization algorithms
//! - **Comprehensive error handling** with detailed diagnostic information

pub(crate) mod cache;
pub mod config;
pub(crate) mod context;
pub(crate) mod energy_grid;
pub mod error;
pub(crate) mod placement;
pub mod progress;
pub mod state;
pub(crate) mod tasks;
pub(crate) mod transaction;
pub(crate) mod utils;
