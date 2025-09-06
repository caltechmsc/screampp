//! # Force Field Module
//!
//! This module provides the core functionality for molecular mechanics force field calculations
//! in the SCREAM++ protein side-chain placement library. It implements energy evaluation,
//! parameter management, and scoring algorithms for molecular interactions.
//!
//! ## Overview
//!
//! The force field module is responsible for computing interaction energies between atoms
//! and molecular groups using classical molecular mechanics potentials. It supports:
//!
//! - **Van der Waals interactions** using Lennard-Jones and Buckingham potentials
//! - **Electrostatic interactions** with Coulomb's law and distance-dependent dielectric
//! - **Hydrogen bonding** using specialized 12-10 potentials
//! - **Energy weighting** based on atom roles (backbone vs sidechain)
//! - **Parameter management** for different force field types and atom types
//!
//! ## Key Components
//!
//! - [`params`] - Force field parameter structures and configuration
//! - [`scoring`] - High-level energy scoring interface for molecular systems
//! - [`term`] - Energy term aggregation and reporting
//! - [`parameterization`] - Automatic assignment of force field parameters to atoms
//!
//! ## Usage
//!
//! The main entry point for energy calculations is the [`scoring::Scorer`] struct,
//! which provides methods to compute interaction energies between atom groups while
//! properly handling bonded exclusions and energy weighting.
//!
//! ```ignore
//! use screampp::core::forcefield::scoring::Scorer;
//!
//! let scorer = Scorer::new(&system, &forcefield);
//! let energy = scorer.score_interaction(query_atoms, environment_atoms)?;
//! ```

pub(crate) mod energy;
pub mod parameterization;
pub mod params;
pub(crate) mod potentials;
pub mod scoring;
pub mod term;
