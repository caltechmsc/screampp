# `scream-core` Developer Documentation

This documentation details the internal architecture, data models, and algorithms of the `scream-core` library. It is intended for developers contributing to or building upon the library's scientific functionalities. For API-level details, refer to the rustdoc documentation.

## Core Design Principles

The design of `scream-core` adheres to three principles:

1. **Separation of Concerns**: A strict three-layer architecture isolates data representation (`core`), stateful logic (`engine`), and high-level procedures (`workflows`).
2. **Performance by Design**: The architecture employs an incremental energy update model, parallelism, and a transactional memory system to ensure CPU and memory efficiency.
3. **Safety and Correctness**: Rust's ownership model is leveraged to eliminate memory-related bugs common in scientific computing.

## Documentation Index

The following documents provide a comprehensive overview of the library's internals.

---

### **1. [Architecture and Data Models](./01_architecture_and_data_models.md)**

- The **Three-Layer Architecture** (Core, Engine, Workflows).
- Core data structures, including `MolecularSystem` and the rationale for `slotmap`-based identifiers.
- The role of `TopologyRegistry` and `RotamerLibrary`.

---

### **2. [Forcefield and Energy Calculation](./02_energy_calculation.md)**

- The decomposition of the physical pairwise energy model into algorithmically efficient components: `Fixed Energy`, `Interaction Energy`, and `Empty Lattice (EL) Energy`.
- The theoretical basis and implementation of the **Flat-Bottom Strategy**.

---

### **3. [Algorithms and Workflows](./03_algorithms_and_workflows.md)**

- Step-by-step analysis of the main `place::run` workflow.
- Detailed breakdown of core algorithms: **Doublet Optimization**, **Simulated Annealing**, and the **Incremental Energy Update Model** (`EnergyGrid`).

---

### **4. [Performance and Memory](./04_performance_and_memory.md)**

- CPU optimization strategies, including **parallelism with Rayon** and heuristic pruning.
- Memory efficiency through the **transactional model (`SystemView`)**, which avoids expensive cloning of the `MolecularSystem`.
