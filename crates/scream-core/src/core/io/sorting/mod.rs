//! Utilities for sorting molecular system components into canonical order.
//!
//! This module contains functionality for organizing atoms, residues, and other
//! molecular components into a consistent, biologically meaningful order for
//! file output and processing. The sorting ensures reproducible results across
//! different molecular file formats and maintains the expected ordering used
//! in computational chemistry and structural biology applications.

pub mod rules;
pub mod sorter;
