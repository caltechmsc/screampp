//! # Rotamers Module
//!
//! This module provides comprehensive functionality for managing protein side-chain rotamers
//! in SCREAM++, enabling efficient and accurate protein structure prediction and refinement.
//!
//! ## Overview
//!
//! The rotamers module implements support for discrete conformational states of amino acid
//! side chains, which are crucial for protein structure modeling. Rotamers represent the
//! most probable orientations of side chains around their rotatable bonds, based on
//! statistical analysis of known protein structures.
//!
//! ## Key Components
//!
//! - [`library`] - Rotamer library management and loading from external files
//! - [`rotamer`] - Core data structures for representing individual rotamers
//!
//! ## Scientific Background
//!
//! Protein side chains can adopt multiple discrete conformations (rotamers) due to
//! rotations around single bonds. The most common rotamers are pre-computed from
//! structural databases and used in:
//!
//! - **Protein structure prediction** - Selecting optimal side-chain conformations
//! - **Protein design** - Exploring sequence-conformation relationships
//! - **Molecular docking** - Avoiding steric clashes during ligand binding
//! - **Structure refinement** - Optimizing side-chain packing in crystal structures
//!
//! ## Usage
//!
//! The rotamers module is typically used through the [`library::RotamerLibrary`] which
//! provides access to pre-computed rotamer conformations for different amino acid types.
//!
//! ```ignore
//! use screampp::core::rotamers::library::RotamerLibrary;
//!
//! let library = RotamerLibrary::load("rotamers.toml", &topology, &forcefield, 1.0)?;
//! let ala_rotamers = library.get_rotamers_for(ResidueType::Alanine);
//! ```

pub mod library;
pub mod rotamer;
