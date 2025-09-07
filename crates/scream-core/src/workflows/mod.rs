//! # Workflows Module
//!
//! This module provides high-level workflow implementations that orchestrate complete
//! optimization processes for protein side-chain placement in SCREAM++.
//!
//! ## Overview
//!
//! Workflows are the top-level entry points for users of SCREAM++. They encapsulate
//! the entire optimization pipeline, from initial setup through final result generation.
//! Each workflow handles resource loading, parameter validation, progress reporting,
//! and result organization, providing a clean and simple API for complex optimization tasks.
//!
//! ## Architecture
//!
//! The module is organized around specific optimization workflows:
//!
//! - **Placement Workflow** ([`place`]) - Complete side-chain conformation optimization
//!   including clash resolution, simulated annealing, and refinement phases.
//!
//! ## Key Capabilities
//!
//! - **End-to-end optimization** from molecular input to optimized conformations
//! - **Resource management** including forcefield, topology, and rotamer library loading
//! - **Progress monitoring** with detailed phase and task reporting
//! - **Result organization** with sorted solutions and energy analysis
//! - **Error handling** with comprehensive diagnostic information
//! - **Flexible configuration** supporting various optimization strategies

pub mod place;
