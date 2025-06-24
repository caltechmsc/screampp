# SCREAM++

**SCREAM++** is an enhanced, high-performance software package for automated protein side-chain placement. It builds upon the scientific foundation of the original SCREAM (Side-Chain Rotamer Excitation Analysis Method), which predicts accurate side-chain conformations using rotamer libraries and a unique flat-bottom potential strategy. This new generation of SCREAM has been completely re-engineered from the ground up in **Rust** for superior memory safety, performance, and modern development practices, replacing the original C++ 98/Python 2 implementation.

The core mission of SCREAM++ is to provide a robust, reliable, and easy-to-use tool for researchers in computational biology, structural biology, and drug design.

## Features

- **High-Accuracy Side-Chain Prediction**: Implements the core **Flat-Bottom Strategy** to accurately model side-chain conformations, compensating for the discrete nature of rotamer libraries.
- **Memory Safety & Performance**: Written in Rust to eliminate entire classes of memory-related bugs and deliver performance comparable to, or exceeding, traditional C++ implementations.
- **Multi-Scheme Support**: Natively handles various charge schemes (CHARMM, AMBER, QEq) and force field parameters, using pre-optimized rotamer libraries for maximum self-consistency.
- **Flexible Optimization Modes**: Supports both global side-chain optimization and focused refinement within a defined **binding site** or region of interest.
- **Modern Tooling**: Managed by Cargo for simple, reproducible builds and dependency management.
- **Multiple Interfaces**:
  - **Standalone CLI**: An easy-to-use command-line interface for standard prediction tasks.
  - **Python Package**: A user-friendly Python API (via PyO3) for scripting, integration, and advanced workflows.
  - **C-Compatible Library**: A C FFI layer for integration with C, C++, and other programming languages.
  - **Native Rust Crate**: A core Rust library (`screampp`) available on crates.io for direct use in other Rust-based scientific computing projects.

## Tech Stack

- **Language**: Rust
- **Supported Languages**: Rust, Python (via PyO3), C/C++ (via FFI)
- **Build System**: Cargo

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
