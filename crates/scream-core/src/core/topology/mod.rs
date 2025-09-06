//! # Topology Module
//!
//! This module provides functionality for managing molecular topology information,
//! including residue structures and atom classifications for protein systems.
//!
//! ## Overview
//!
//! The topology module defines the structural organization of molecular systems,
//! particularly focusing on amino acid residues and their constituent atoms.
//! It provides data structures and utilities for:
//!
//! - **Residue topology definitions** - Specifying anchor and sidechain atoms for each residue type
//! - **Topology registry** - Loading and accessing residue topology information from configuration files
//! - **Atom classification** - Distinguishing between backbone and sidechain atoms
//!
//! ## Key Components
//!
//! - [`registry`] - Topology registry for loading and accessing residue definitions
//!
//! ## Usage
//!
//! Residue topologies are typically loaded from TOML configuration files and used
//! to classify atoms during molecular system construction and analysis.
//!
//! ```ignore
//! use screampp::core::topology::registry::TopologyRegistry;
//!
//! let registry = TopologyRegistry::load("topology.toml")?;
//! let ala_topology = registry.get("ALA").unwrap();
//! ```

pub mod registry;
