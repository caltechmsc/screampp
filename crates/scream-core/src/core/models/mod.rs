//! # Core Models Module
//!
//! This module contains the fundamental data structures and models used to represent
//! molecular systems in SCREAM++, providing the foundation for all molecular modeling operations.
//!
//! ## Overview
//!
//! The models module defines the core abstractions for representing molecular structures,
//! including atoms, residues, chains, and their topological relationships. These models
//! are designed to:
//!
//! - **Represent molecular structure** - Complete description of atomic coordinates and connectivity
//! - **Support efficient operations** - Optimized data structures for computational algorithms
//! - **Enable extensibility** - Flexible design for different molecular types and properties
//! - **Maintain type safety** - Strong typing for molecular data integrity
//!
//! ## Key Components
//!
//! - [`atom`] - Individual atom representation with coordinates, types, and properties
//! - [`residue`] - Amino acid residue structure and classification
//! - [`chain`] - Polypeptide chain organization and metadata
//! - [`system`] - Complete molecular system with all components and relationships
//! - [`topology`] - Bond connectivity and molecular topology information
//! - [`ids`] - Unique identifier types for atoms, residues, and chains
//!
//! ## Usage
//!
//! The models form the backbone of molecular data representation in SCREAM++.
//! Most operations start with constructing or manipulating these core structures.
//!
//! ```ignore
//! use screampp::core::models::{system::MolecularSystem, atom::Atom};
//!
//! let mut system = MolecularSystem::new();
//! let chain_id = system.add_chain('A', ChainType::Protein);
//! let residue_id = system.add_residue(chain_id, 1, "ALA", None)?;
//!
//! let atom = Atom::new("CA", residue_id, Point3::new(0.0, 0.0, 0.0));
//! system.add_atom_to_residue(residue_id, atom)?;
//! ```

pub mod atom;
pub mod chain;
pub mod ids;
pub mod residue;
pub mod system;
pub mod topology;
