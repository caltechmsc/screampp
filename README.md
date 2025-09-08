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
  - **Standalone CLI**: A powerful and easy-to-use command-line interface for all prediction tasks.
  - **Native Rust Crate**: The core library (`screampp`) is available on [crates.io](https://crates.io/crates/screampp) for direct use in other Rust-based scientific computing projects.
  - **Python Package (Future)**: A user-friendly Python API (via PyO3) for scripting and integration is planned for future releases.
  - **C-Compatible Library (Future)**: A C FFI layer for integration with C, C++, and other languages is planned for future releases.

## Getting Started

There are two main ways to use SCREAM++: as a standalone command-line tool or as a library in your own Rust project.

### 1. Using the Command-Line Interface (CLI)

This is the recommended method for most users.

#### Step 1: Download the executable

Go to the [**GitHub Releases page**](https://github.com/caltechmsc/screampp/releases) and download the pre-compiled binary for your operating system (Linux, macOS, Windows). Unzip the archive.

#### Step 2: Download the required data files

The first time you run the CLI, you must download the necessary forcefield and rotamer library files. This is a one-time setup.

```bash
# On Linux/macOS
./scream data download

# On Windows (Command Prompt)
scream.exe data download
```

This command will fetch and unpack the data to a default location on your system.

#### Step 3: Run a prediction

You are now ready to run a side-chain placement job.

```bash
# Optimize all side-chains in input.bgf and save the result
./scream place -i path/to/input.bgf -o path/to/output.bgf
```

For detailed instructions on all commands, options, and advanced configuration, please refer to the [**CLI User Manual**](docs/cli/USAGE.md).

### 2. Using as a Rust Library

If you are a Rust developer, you can add `screampp` as a dependency to your project.

#### Step 1: Add to `Cargo.toml`

```toml
[dependencies]
screampp = "0.5.0"
```

#### Step 2: Use in your code

You can now use the high-level `workflows` API to perform side-chain placement programmatically.

```rust
use screampp::core::io::bgf::BgfFile;
use screampp::engine::config::{PlacementConfigBuilder, ResidueSelection, ConvergenceConfig};
use screampp::engine::progress::ProgressReporter;
use screampp::workflows::place;

fn run_placement() -> Result<(), Box<dyn std::error::Error>> {
    let (system, metadata) = BgfFile::read_from_path("path/to/input.bgf")?;

    let config = PlacementConfigBuilder::new()
        .forcefield_path("path/to/data/forcefield/dreiding-lj-12-6-0.4.toml")
        .delta_params_path("path/to/data/delta/delta-rmsd-1.0.csv")
        .s_factor(1.1)
        .rotamer_library_path("path/to/data/rotamers/charmm@rmsd-1.0.toml")
        .topology_registry_path("path/to/data/topology/registry.toml")
        .residues_to_optimize(ResidueSelection::All)
        .max_iterations(100)
        .num_solutions(1)
        .include_input_conformation(true)
        .final_refinement_iterations(2)
        .convergence_config(ConvergenceConfig {
            energy_threshold: 0.01,
            patience_iterations: 5,
        })
        .build()?;

    let reporter = ProgressReporter::new(); // Or provide a callback for progress updates
    let result = place::run(&system, &config, &reporter)?;

    if let Some(best_solution) = result.solutions.first() {
        println!("Best energy: {:.2} kcal/mol", best_solution.total_energy);
    }

    Ok(())
}
```

## Documentation

Comprehensive documentation is available for both users and developers.

- **For Users**:

  - [**CLI User Manual**](docs/cli/USAGE.md): A complete guide to installing and using the `scream` command-line tool, including detailed explanations of all commands, configuration options, and practical examples.

- **For Developers**:

  - [**Rust Library API Docs (docs.rs)**](https://docs.rs/screampp): The official, versioned API documentation for the `scream-core` (`screampp`) crate, generated by `rustdoc`. This is the best resource for understanding the public API of the library.
  - **Developer Documentation**: In-depth documentation covering the architecture, data models, algorithms, and design philosophy of the entire project. This is essential reading for anyone looking to contribute to or deeply understand the internals of SCREAM++.
    - [**`scream-core` Developer Docs**](docs/dev/core/README.md): Details the internal architecture, data models, and algorithms of the `scream-core` library.
    - [**`scream-cli` Developer Docs**](docs/dev/cli/README.md): Provides a comprehensive technical breakdown of the `scream-cli` crate, including execution flow, configuration handling, and data management.

## Tech Stack

- **Core Language**: Rust
- **Build System**: Cargo
- **Planned Interfaces**: Python (via PyO3), C/C++ (via FFI)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
