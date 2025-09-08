# SCREAM++ CLI: User Manual

Welcome to the SCREAM++ command-line interface! This manual will guide you through configuring and using the `scream` CLI for protein side-chain placement and structural optimization.

## Table of Contents

- [SCREAM++ CLI: User Manual](#scream-cli-user-manual)
  - [Table of Contents](#table-of-contents)
  - [Setting Up the Data Directory](#setting-up-the-data-directory)
    - [Downloading Data Files](#downloading-data-files)
    - [Managing the Data Path](#managing-the-data-path)
  - [Core Functionality: Side-Chain Placement (`scream place`)](#core-functionality-side-chain-placement-scream-place)
    - [Basic Usage](#basic-usage)
    - [Argument Reference](#argument-reference)
      - [Core Arguments](#core-arguments)
      - [Forcefield Overrides](#forcefield-overrides)
      - [Optimization Overrides](#optimization-overrides)
      - [General Arguments](#general-arguments)
  - [The Configuration File (`config.toml`)](#the-configuration-file-configtoml)
    - [Configuration Structure](#configuration-structure)
    - [Example Configuration File](#example-configuration-file)
    - [Detailed Configuration Options](#detailed-configuration-options)
      - [The `[forcefield]` Table](#the-forcefield-table)
      - [The `[sampling]` Table](#the-sampling-table)
      - [The `[optimization]` Table](#the-optimization-table)
      - [The `[residues-to-optimize]` Table](#the-residues-to-optimize-table)
  - [Practical Examples (Use Cases)](#practical-examples-use-cases)
    - [Example 1: Simple Global Optimization](#example-1-simple-global-optimization)
    - [Example 2: Using a Specific Rotamer Library and s-factor](#example-2-using-a-specific-rotamer-library-and-s-factor)
    - [Example 3: Optimizing a Ligand Binding Pocket](#example-3-optimizing-a-ligand-binding-pocket)
    - [Example 4: Generating Multiple Solutions with Templated Naming](#example-4-generating-multiple-solutions-with-templated-naming)
  - [Configuration Reference Table](#configuration-reference-table)

---

## Setting Up the Data Directory

SCREAM++ relies on a set of data files, including forcefield parameters, rotamer libraries, and topology definitions. Before you can use the `place` command, you must download this data.

### Downloading Data Files

> **Important:** This is a **mandatory first step**. SCREAM++ cannot run without these data files. (Unless you plan to provide all necessary files manually, which is not recommended for typical users.)

This is the **mandatory first step** for using SCREAM++. Execute the following command to download and automatically unpack all required data files to their default location:

```sh
scream data download
```

This command fetches the same versioned data bundle from the official GitHub Releases page. If you have downloaded the data before but wish to force an update, use the `--force` flag:

```sh
scream data download --force
```

### Managing the Data Path

By default, data files are stored in your operating system's standard user data directory (e.g., `~/.local/share/screampp` on Linux).

- **View the current data path**:

  ```sh
  scream data path
  ```

- **Set a custom data path**:

  If you wish to store the data files in a different location (e.g., a project directory or a shared server path), use `set-path`:

  ```sh
  scream data set-path /path/to/my/custom/data/directory
  ```

  After setting this, future `download` and `place` commands will use this new path.

- **Reset to the default path**:

  To revert to the default data storage location, run:

  ```sh
  scream data reset-path
  ```

---

## Core Functionality: Side-Chain Placement (`scream place`)

The `place` command is the heart of SCREAM++, used to perform side-chain optimization on an input protein structure.

### Basic Usage

The simplest usage involves providing an input file and an output file path:

```sh
scream place -i input_protein.bgf -o optimized_protein.bgf
```

This command will optimize all side-chains in `input_protein.bgf` using built-in default parameters and save the lowest-energy structure to `optimized_protein.bgf`.

### Argument Reference

You can precisely control the optimization process via command-line arguments or a configuration file. Command-line arguments will always override settings in a configuration file.

#### Core Arguments

- `-i, --input <PATH>` (**Required**):
  Specifies the path to the input molecular structure file. Currently, only BGF format is supported.

- `-o, --output <PATH_TEMPLATE>` (**Required**):
  Specifies the path for the output file(s). This is a powerful template that allows for dynamic file naming based on the optimization results. Supported placeholders include:

  - `{i}` or `{n}`: The rank of the solution (1 is best).
  - `{total}` or `{N}`: The total number of solutions generated.
  - `{energy}` or `{total_energy}`: The total energy of the solution (kcal/mol).
  - `{score}` or `{opt_score}`: The internal optimization score of the solution.
  - _Example_: `-o "solution_{i}_of_{total}_E_{energy}.bgf"`

- `-c, --config <PATH>`:
  Specifies the path to a configuration file in TOML format. This is the recommended method for complex setups (e.g., selecting specific residues).

#### Forcefield Overrides

These arguments allow you to quickly override forcefield settings from the configuration file.

- `-s, --s-factor <FLOAT>`:
  Sets the `s-factor` for the flat-bottom potential. This value influences the energy function's tolerance for minor atomic clashes. The default is typically `1.1`.

- `--forcefield-path <NAME_OR_PATH>`:
  Specifies the forcefield parameter file. Can be a local file path or a logical name (e.g., `'lj-12-6@0.4'`).

- `--delta-params-path <NAME_OR_PATH>`:
  Specifies the `delta` parameter file for the flat-bottom potential correction. Can be a local file path or a logical name (e.g., `'rmsd-1.0'`).

- `-t, --topology-registry <NAME_OR_PATH>`:
  Specifies the residue topology registry file. The default `'default'` is usually sufficient.

#### Optimization Overrides

- `-l, --rotamer-library <NAME_OR_PATH>`:
  Specifies the rotamer library to use. This is a critical parameter as it defines the conformational space for side-chain sampling. Can be a local file path or a logical name (e.g., `'charmm@rmsd-1.0'`). Available charge schemes include `amber`, `charmm`, `qeq`, `amber-n`, `charmm-n`, and `qeq-n`. Available resolutions (diversity) include `rmsd-0.1` to `rmsd-5.0` (0.1 increments) and `all-torsion`.

- `-n, --num-solutions <INT>`:
  Specifies the number of lowest-energy solutions to generate and save. Defaults to `1`.

- `--max-iterations <INT>`:
  Sets the maximum number of iterations for the clash resolution phase.

- `--with-input-conformation` / `--no-input-conformation`:
  Force the inclusion or exclusion of the original side-chain conformation from the input structure as a candidate.

- `--no-refinement`:
  Disables the final singlet optimization (refinement) stage. This can speed up calculations but may slightly reduce accuracy.

- `--no-annealing`:
  Disables the simulated annealing process, even if it is enabled in the configuration file.

- `-S, --set <KEY=VALUE>`:
  Directly set a configuration value from the command line for quick overrides. Example: `-S optimization.max-iterations=200`.

#### General Arguments

- `-v, --verbose`: Increase the log verbosity. `-v` (INFO), `-vv` (DEBUG), `-vvv` (TRACE).
- `-q, --quiet`: Suppress all log output except for errors.
- `--log-file <PATH>`: Write logs to a specified file in addition to console output.
- `-j, --threads <NUM>`: Set the number of threads for parallel computation. Defaults to all available CPU cores. (If you are using a HPC with a multi-threaded job, please ignore this and let the system manage CPU usage - SCREAM++ will automatically detect the number of available cores available to your job.)

---

## The Configuration File (`config.toml`)

For complex or reproducible tasks, using a TOML configuration file is highly recommended.

### Configuration Structure

The configuration file is organized into four main sections (TOML tables):

1. `[forcefield]`: Parameters related to the energy function.
2. `[sampling]`: Parameters related to conformational sampling.
3. `[optimization]`: Parameters to control the optimization algorithm.
4. `[residues-to-optimize]`: Defines the scope of the optimization.

### Example Configuration File

Here is an example configuration file with all common options and detailed comments. You can copy this and modify it to suit your needs.

```toml
# =============================================================================
# SCREAM++ Example Configuration File
# =============================================================================

# -----------------------------------------------------------------------------
# [forcefield] - Energy Function & Parameters
# -----------------------------------------------------------------------------
[forcefield]

# -- Core Parameters --

# Path or logical name for the main forcefield parameter file.
# Logical names: 'exp-6@0.4', 'lj-12-6@0.4', etc.
# Default: "exp-6@0.4"
forcefield-path = "exp-6@0.4"

# Path or logical name for the flat-bottom delta parameters.
# The diversity (e.g., "rmsd-1.0") should match your rotamer library.
# Default: "rmsd-1.0"
delta-params-path = "rmsd-1.0"

# The 's-factor' for the flat-bottom potential. This is a critical parameter
# that tunes the tolerance for atomic clashes.
# Default: 1.1
s-factor = 1.1

# -- [Optional] Advanced Energy Weighting --
#
# This section allows you to scale energy terms for interactions between
# different types of atoms (Atom Roles: Backbone, Sidechain, Ligand, Water, Other).
# By default, all weights are 1.0.
#
# [[forcefield.energy-weights.rules]]
# groups = ["Backbone", "Sidechain"]
# weights = { vdw = 0.8, coulomb = 0.8, hbond = 1.0 }
#
# [[forcefield.energy-weights.rules]]
# groups = ["Sidechain", "Ligand"]
# weights = { vdw = 0.5, coulomb = 1.0, hbond = 1.2 }


# -----------------------------------------------------------------------------
# [sampling] - Conformational Sampling
# -----------------------------------------------------------------------------
[sampling]

# Path or logical name for the rotamer library.
# The diversity (e.g., "rmsd-1.0") should match your delta-params-path.
# Logical names: 'charmm@rmsd-1.0', 'amber@rmsd-1.0', etc.
# Default: "charmm@rmsd-1.0"
rotamer-library = "charmm@rmsd-1.0"


# -----------------------------------------------------------------------------
# [optimization] - Algorithm Control
# -----------------------------------------------------------------------------
[optimization]

# Number of lowest-energy, unique solutions to generate and save.
# Default: 1
num-solutions = 1

# Maximum number of iterations for the primary clash-resolution loop.
# Default: 100
max-iterations = 100

# If true, the original side-chain conformation from the input structure
# will be included as a candidate during the optimization.
# Default: true
include-input-conformation = true

# Number of refinement iterations (singlet optimization) to perform after the
# main clash-resolution loop has converged. Set to 0 to disable.
# Default: 2
final-refinement-iterations = 2

# -- [Optional] Simulated Annealing --
#
# To enable simulated annealing for better global energy landscape exploration,
# uncomment this entire section. This may improve results but will increase runtime.
#
# [optimization.simulated-annealing]
# initial-temperature = 5.0      # Starting temperature (in energy units).
# final-temperature = 0.1        # Temperature at which to stop the annealing.
# cooling-rate = 0.9             # Multiplicative factor to decrease temperature (e.g., T_new = T_old * 0.9).
# steps-per-temperature = 100    # Number of Monte Carlo moves to attempt at each temperature step.

# -- Convergence Criteria --
#
# Defines the conditions for stopping the clash-resolution loop.
#
[optimization.convergence]
# The loop will stop if the best energy found does not improve by at least
# this amount (in kcal/mol) over a 'patience' number of iterations.
# Default: 0.01
energy-threshold = 0.01

# The number of consecutive iterations without sufficient energy improvement
# before the algorithm is considered to have converged.
# Default: 5
patience-iterations = 5


# -----------------------------------------------------------------------------
# [residues-to-optimize] - Defines the Scope of Optimization
# -----------------------------------------------------------------------------
[residues-to-optimize]

# TYPE 1: Optimize all residues in the protein.
type = "all"

# TYPE 2: Optimize a specific list of residues.
# type = "list"
# # 'include' defines a whitelist. If 'include' is empty, all residues are selected.
# include = [
#   { chain-id = 'A', residue-number = 25 },
#   { chain-id = 'A', residue-number = 101 },
# ]
# # 'exclude' defines a blacklist that overrides the selection.
# exclude = [
#   { chain-id = 'A', residue-number = 50 },
# ]

# TYPE 3: Optimize residues within a radius of a ligand.
# type = "ligand-binding-site"
# # Specify the ligand's location.
# [residues-to-optimize.ligand-residue]
# chain-id = 'X'
# residue-number = 999
# # Define the radius in Angstroms from any heavy atom of the ligand.
# radius-angstroms = 5.0
```

### Detailed Configuration Options

#### The `[forcefield]` Table

- `forcefield-path`: Specifies forcefield parameters. Logical name format: `<potential_type>@<version>`.
- `delta-params-path`: Specifies `delta` parameters. Logical name format: `rmsd-<value>` or `all-torsion`. The `<value>` should match the rotamer library's `diversity`.
- `s-factor`: A key parameter. It is strongly recommended to use a value consistent with `delta-params-path` and `rotamer-library`.
- `energy-weights`: An advanced option. Allows you to assign different weights to interactions between different atom roles (e.g., Backbone-Sidechain).

#### The `[sampling]` Table

- `rotamer-library`: Specifies the rotamer library. Logical name format: `<scheme>@<diversity>`, e.g., `charmm@rmsd-1.0`. `scheme` can be `charmm`, `amber`, etc. The `diversity` should match `delta-params-path`.

#### The `[optimization]` Table

- `simulated-annealing`: If this section is present, simulated annealing will be enabled. This helps escape local energy minima but increases computation time.
- `convergence`: Controls when the iterative algorithm stops.

#### The `[residues-to-optimize]` Table

This is the core section for defining the scope of your optimization.

- `type = "all"`: The simplest case; all residues are optimized.
- `type = "list"`:
  - `include`: Defines a whitelist of residues to optimize. If `include` is empty, all residues are selected by default.
  - `exclude`: Defines a blacklist. These residues will **not** be optimized, even if they are in the `include` list.
- `type = "ligand-binding-site"`:
  - `ligand-residue`: Specifies the ligand residue.
  - `radius-angstroms`: Defines a spherical region. Any protein residue with a heavy atom (non-hydrogen) inside this sphere will be selected for optimization.

---

## Practical Examples (Use Cases)

Assume we have an input file `protein.bgf`.

### Example 1: Simple Global Optimization

**Goal**: Quickly optimize all side-chains in the protein and save the best result.

```sh
scream place -i protein.bgf -o protein_optimized.bgf
```

> This uses all default parameters and is ideal for a quick initial assessment.

### Example 2: Using a Specific Rotamer Library and s-factor

**Goal**: Perform optimization using the CHARMM charge scheme and a finer rotamer library (`rmsd-0.6`).

```sh
scream place \
    -i protein.bgf \
    -o protein_charmm_0.6.bgf \
    -l charmm@rmsd-0.6 \
    --delta-params-path rmsd-0.6 \
    -s 1.2
```

> **Note**: The `diversity` (`rmsd-0.6`) of the `rotamer-library` and `delta-params-path` must match. The `s-factor` (`1.2`) is a recommended value optimized for this diversity.

### Example 3: Optimizing a Ligand Binding Pocket

**Goal**: The input file `complex.bgf` contains a protein and a ligand (chain X, residue number 999). We want to optimize only the protein side-chains within 5 Å of the ligand.

First, create a configuration file `pocket_opt.toml`:

```toml
# pocket_opt.toml
[residues-to-optimize]
type = "ligand-binding-site"
radius-angstroms = 5.0

[residues-to-optimize.ligand-residue]
chain-id = 'X'
residue-number = 999
```

Then, run the command:

```sh
scream place -i complex.bgf -o complex_pocket_optimized.bgf -c pocket_opt.toml
```

> **Note**: Using energy weights is recommended for binding site optimization, as it allows advanced control and fine-tuning. For example, you can increase the strength of sidechain-ligand hydrogen bonds or reduce the van der Waals interactions between sidechains and backbone to achieve specific optimization goals.

### Example 4: Generating Multiple Solutions with Templated Naming

**Goal**: We are uncertain which conformation is best and want to generate the top 3 lowest-energy solutions, naming the output files based on their rank and energy.

```sh
scream place \
    -i protein.bgf \
    -o "protein_sol_{i}_E_{energy}.bgf" \
    -n 3
```

This will generate files with names like:

- `protein_sol_1_E_-1234.56.bgf` (The lowest energy solution)
- `protein_sol_2_E_-1232.10.bgf`
- `protein_sol_3_E_-1230.05.bgf`

---

## Configuration Reference Table

This table provides a comprehensive mapping between the `scream place` command-line arguments and their equivalent settings in the `config.toml` file. Use it as a quick reference for all available configuration options. Command-line arguments always take precedence over the configuration file.

| CLI Argument (Short) | CLI Argument (Long)         | `config.toml` Key                                        | Value Type            | Default           | Description                                                                                |
| :------------------- | :-------------------------- | :------------------------------------------------------- | :-------------------- | :---------------- | :----------------------------------------------------------------------------------------- |
| `-i`                 | `--input`                   | _(N/A)_                                                  | File Path             | **Required**      | Path to the input molecular structure file (.bgf).                                         |
| `-o`                 | `--output`                  | _(N/A)_                                                  | File Path Template    | **Required**      | Path for the output file(s), supports templating.                                          |
| `-c`                 | `--config`                  | _(N/A)_                                                  | File Path             | (None)            | Path to the main TOML configuration file.                                                  |
| **---**              | **---**                     | **Forcefield Settings**                                  | **---**               | **---**           | **---**                                                                                    |
| `-s`                 | `--s-factor`                | `forcefield.s-factor`                                    | Float                 | `1.1`             | The scaling factor (`s`) for the flat-bottom potential.                                    |
|                      | `--forcefield-path`         | `forcefield.forcefield-path`                             | String (Path/Logical) | `exp-6@0.4`       | Path or logical name for the forcefield parameters.                                        |
|                      | `--delta-params-path`       | `forcefield.delta-params-path`                           | String (Path/Logical) | `rmsd-1.0`        | Path or logical name for the flat-bottom delta parameters.                                 |
| `-t`                 | `--topology-registry`       | `topology-registry-path`                                 | String (Path/Logical) | `default`         | Path or logical name for the residue topology registry.                                    |
| **---**              | **---**                     | **Sampling & Optimization**                              | **---**               | **---**           | **---**                                                                                    |
| `-l`                 | `--rotamer-library`         | `sampling.rotamer-library`                               | String (Path/Logical) | `charmm@rmsd-1.0` | Path or logical name for the rotamer library.                                              |
| `-n`                 | `--num-solutions`           | `optimization.num-solutions`                             | Integer               | `1`               | The number of top solutions to generate and save.                                          |
|                      | `--max-iterations`          | `optimization.max-iterations`                            | Integer               | `100`             | Maximum number of iterations for the clash resolution phase.                               |
|                      | `--with-input-conformation` | `optimization.include-input-conformation`                | Boolean Flag          | `true`            | Forces **inclusion** of the input conformation as a candidate.                             |
|                      | `--no-input-conformation`   | `optimization.include-input-conformation`                | Boolean Flag          | `true`            | Forces **exclusion** of the input conformation.                                            |
|                      | `--no-refinement`           | `optimization.final-refinement-iterations`               | Boolean Flag          | `2`               | Disables the final refinement stage (sets iterations to 0).                                |
|                      | `--no-annealing`            | `optimization.simulated-annealing`                       | Boolean Flag          | (Disabled)        | Disables simulated annealing, even if set in the config file.                              |
| **---**              | **---**                     | **Convergence Settings**                                 | **---**               | **---**           | **---**                                                                                    |
|                      | _(N/A)_                     | `optimization.convergence.energy-threshold`              | Float                 | `0.01`            | Stop if energy improvement is less than this value (kcal/mol).                             |
|                      | _(N/A)_                     | `optimization.convergence.patience-iterations`           | Integer               | `5`               | Number of iterations without improvement before stopping.                                  |
| **---**              | **---**                     | **Simulated Annealing**                                  | **---**               | **---**           | **---**                                                                                    |
|                      | _(N/A)_                     | `optimization.simulated-annealing.initial-temperature`   | Float                 | (N/A)             | Starting temperature for the annealing schedule.                                           |
|                      | _(N/A)_                     | `optimization.simulated-annealing.final-temperature`     | Float                 | (N/A)             | Ending temperature for the annealing schedule.                                             |
|                      | _(N/A)_                     | `optimization.simulated-annealing.cooling-rate`          | Float (0-1)           | (N/A)             | Multiplicative factor for decreasing temperature.                                          |
|                      | _(N/A)_                     | `optimization.simulated-annealing.steps-per-temperature` | Integer               | (N/A)             | Number of moves to attempt at each temperature step.                                       |
| **---**              | **---**                     | **Residue Selection**                                    | **---**               | **---**           | **---**                                                                                    |
|                      | _(N/A)_                     | `residues-to-optimize`                                   | Table                 | `type = "all"`    | Defines which residues to optimize. See manual for `list` and `ligand-binding-site` types. |
| **---**              | **---**                     | **Advanced Overrides**                                   | **---**               | **---**           | **---**                                                                                    |
| `-S`                 | `--set`                     | _(Various)_                                              | String (`KEY=VALUE`)  | (None)            | Overrides a specific config value directly. E.g., `-S optimization.num-solutions=5`.       |
| **---**              | **---**                     | **General Settings**                                     | **---**               | **---**           | **---**                                                                                    |
| `-j`                 | `--threads`                 | _(N/A)_                                                  | Integer               | (CPU Cores)       | Number of threads for parallel computation.                                                |
| `-v`                 | `--verbose`                 | _(N/A)_                                                  | Count Flag            | (Off)             | Increase verbosity (`-v`, `-vv`, `-vvv`).                                                  |
| `-q`                 | `--quiet`                   | _(N/A)_                                                  | Boolean Flag          | (Off)             | Suppress all output except errors.                                                         |
|                      | `--log-file`                | _(N/A)_                                                  | File Path             | (None)            | Path to write a detailed log file.                                                         |
