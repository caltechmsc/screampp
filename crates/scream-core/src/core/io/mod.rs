//! # Input/Output Module
//!
//! This module provides comprehensive input/output functionality for molecular file formats
//! used in computational chemistry and structural biology applications.
//!
//! ## Overview
//!
//! The I/O module enables SCREAM++ to read from and write to various molecular structure
//! file formats, providing a unified interface for molecular data exchange. It supports:
//!
//! - **File format parsing** - Reading molecular structures from standard formats
//! - **File format writing** - Exporting molecular systems to various output formats
//! - **Canonical ordering** - Consistent sorting of atoms and residues for reproducible output
//! - **Trait-based design** - Extensible interface for adding new file format support
//!
//! ## Key Components
//!
//! - [`bgf`] - Implementation for BGF (BioGraf) file format I/O
//! - [`traits`] - Common traits defining the molecular file I/O interface
//!
//! ## Usage
//!
//! The module provides a trait-based approach to file I/O, allowing different formats
//! to implement a common interface for reading and writing molecular structures.
//!
//! ```ignore
//! use screampp::core::io::{bgf::BgfFile, traits::MolecularFile};
//! use std::fs::File;
//! use std::io::BufReader;
//!
//! // Read a molecular structure from a BGF file
//! let file = File::open("molecule.bgf")?;
//! let mut reader = BufReader::new(file);
//! let (system, metadata) = BgfFile::read_from(&mut reader)?;
//!
//! // Write a molecular system to a BGF file
//! let mut file = File::create("output.bgf")?;
//! BgfFile::write_system_to(&system, &mut file)?;
//! ```

pub mod bgf;
pub(crate) mod sorting;
pub mod traits;
