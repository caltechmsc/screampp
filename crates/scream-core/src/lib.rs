//! # SCREAM++ Core Library
//!
//! A modernized, high-performance library for protein side-chain placement and structure redesign.
//!
//! ## Architecture
//!
//! The library is structured into three primary public modules, designed to be used at different
//! levels of abstraction:
//!
//! - [`workflows`]: The highest-level API.
//!
//! - [`engine`]: The configuration and state management layer.
//!
//! - [`core`]: The foundational layer. It contains the fundamental data models, file I/O traits, and lower-level computational tools.

pub mod core;
pub mod engine;
pub mod workflows;
