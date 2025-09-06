//! Provides input/output functionality for molecular file formats.
//!
//! This module contains implementations for reading and writing various molecular
//! structure file formats commonly used in computational chemistry and structural
//! biology. It provides a unified trait-based interface for file I/O operations
//! and includes utilities for canonical ordering of molecular components.

pub mod bgf;
pub(crate) mod sorting;
pub mod traits;
