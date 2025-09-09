//! # SCREAM++ Core Library
//!
//! A modernized, high-performance library for protein side-chain placement and structure redesign,
//! based on the scientific principles of the SCREAM method.
//!
//! ## Architectural Philosophy
//!
//! The library is designed with a strict three-layer architecture to ensure a clear separation of concerns,
//! making it modular, testable, and extensible.
//!
//! - **[`core`]: The Foundation.** Contains stateless data models (`MolecularSystem`),
//!   pure mathematical representations of the forcefield (`potentials`, `scoring`), and I/O utilities.
//!
//! - **[`engine`]: The Logic Core.** This stateful layer orchestrates the optimization process.
//!   It includes high-performance data structures like `EnergyGrid` for incremental updates,
//!   `SystemView` for transactional modifications, and the implementation of optimization
//!   algorithms (e.g., `doublet_optimization`).
//!
//! - **[`workflows`]: The Public API.** This is the highest-level, user-facing layer. It ties the
//!   `engine` and `core` together to execute complete scientific procedures, such as side-chain
//!   placement. It provides a simple and powerful entry point for end-users of the library.

pub mod core;
pub mod engine;
pub mod workflows;
