# `scream-core` Developer Documentation

Welcome to the developer documentation for `scream-core`, the scientific heart of the SCREAM++ project. This library contains all the fundamental data structures, forcefield implementations, and optimization algorithms for protein side-chain placement.

The purpose of this documentation is to explain the **"why"** behind the code—our architectural decisions, algorithmic strategies, and performance considerations. It is intended to be read alongside the in-code API documentation (`cargo doc`).

## Core Philosophy

The design of `scream-core` is guided by three core principles:

1. **Layered Architecture**: A strict separation of concerns into three layers (`core`, `engine`, `workflows`) makes the library modular, testable, and extensible. Data representation is cleanly separated from algorithmic logic and high-level procedures.

2. **Data-Centric Design**: Complex algorithms are built around high-performance data structures. The `MolecularSystem` provides a robust and safe representation of molecular data, while the `EnergyGrid` enables efficient, incremental energy calculations, which is the key to the engine's speed.

3. **Performance by Design**: This is achieved through aggressive parallelism of independent tasks with `rayon`, efficient memory management via Rust's ownership model, and transactional, clone-free system modifications.

## Recommended Reading Path

For developers new to `scream-core`, we recommend reading the detailed documentation in the following order to build a comprehensive understanding of the library:

1. **[Architecture and Data Models](./01_architecture_and_data_models.md)**: Start here to understand the foundational structure of the library and how molecular systems are represented in memory.

2. **[Forcefield and Energy Calculation](./02_energy_calculation.md)**: Once you understand the data structures, this document explains the scientific and mathematical basis for how we calculate energies, including the core "Flat-Bottom Strategy."

3. **[Algorithms and Workflows](./03_algorithms_and_workflows.md)**: This document connects the static data models and energy functions into a dynamic process. It details the step-by-step logic of the main side-chain placement algorithm.

4. **[Performance and Memory](./04_performance_and_memory.md)**: Finally, read this to understand the key optimizations, parallelism strategies, and memory management techniques that make `scream-core` fast and efficient.

## Documentation Index

Here is a quick reference to the detailed documentation for each major component of the library.

---

### **[1. Architecture and Data Models](./01_architecture_and_data_models.md)**

- Explains the three-layer (`core`, `engine`, `workflows`) design philosophy.
- Details the structure of `MolecularSystem` and the rationale behind using `slotmap` for stable identifiers.
- Describes how structural knowledge is encoded in `TopologyRegistry` and `RotamerLibrary`.

---

### **[2. Forcefield and Energy Calculation](./02_energy_calculation.md)**

- Starts from the fundamental physics of pairwise non-bonded energies.
- Breaks down the total energy into algorithmically efficient components: `E_fixed`, `E_interaction`, and `E_EL`.
- Provides a detailed explanation of the core scientific concept: the **Flat-Bottom Strategy** and the role of the `Δ` parameter.

---

### **[3. Algorithms and Workflows](./03_algorithms_and_workflows.md)**

- Provides a step-by-step walkthrough of the main `place::run` workflow, from setup to final results.
- Details the core optimization loop, including clash detection, doublet optimization, and optional simulated annealing.
- Explains the key algorithms that enable performance, such as the incremental energy updates via `EnergyGrid` and the transactional `SystemView` model.

---

### **[4. Performance and Memory](./04_performance_and_memory.md)**

- Discusses CPU performance optimizations, including the $O(N)$ incremental energy model and parallelism with `rayon`.
- Explains memory efficiency strategies, focusing on the clone-free transactional model (`SystemView`) that avoids costly memory allocations.
