//! # Core Module
//!
//! This module provides the fundamental building blocks and algorithms for protein
//! side-chain placement and molecular modeling in SCREAM++, serving as the computational
//! core of the library.
//!
//! ## Overview
//!
//! The core module implements the essential data structures, algorithms, and utilities
//! required for automated protein side-chain conformation prediction. It provides a
//! complete framework for representing molecular systems, computing interaction energies,
//! and managing conformational libraries.
//!
//! ## Architecture
//!
//! The module is organized into specialized submodules that handle different aspects
//! of molecular modeling:
//!
//! - **Molecular Representation** ([`models`]) - Data structures for atoms, residues, chains, and systems
//! - **Energy Calculations** ([`forcefield`]) - Force field parameters and energy computation
//! - **File I/O** ([`io`]) - Reading/writing molecular file formats with canonical ordering
//! - **Structural Knowledge** ([`topology`]) - Residue topology definitions and atom classification
//! - **Conformational Libraries** ([`rotamers`]) - Pre-computed side-chain rotamer collections
//!
//! ## Key Capabilities
//!
//! - **Complete molecular system representation** with efficient data structures
//! - **Molecular mechanics energy calculations** using classical force fields
//! - **Multi-format file I/O** with consistent atom/residue ordering
//! - **Rotamer library management** for side-chain conformation sampling
//! - **Topology-aware atom classification** for backbone/sidechain distinction
//! - **Extensible force field support** for different parameter sets
//!
//! ## Scientific Foundation
//!
//! The core module implements algorithms based on established computational chemistry
//! principles:
//!
//! - **Molecular mechanics** for energy minimization and conformational analysis
//! - **Rotamer libraries** derived from statistical analysis of protein structures
//! - **Force field methods** including Lennard-Jones, Coulomb, and hydrogen bonding potentials
//! - **Topology-based modeling** respecting molecular connectivity and stereochemistry

pub mod forcefield;
pub mod io;
pub mod models;
pub mod rotamers;
pub mod topology;
