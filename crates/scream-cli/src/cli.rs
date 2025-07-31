use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

const HELP_TEMPLATE: &str = "\
{before-help}{name} {version}
{author-with-newline}{about-with-newline}
{usage-heading} {usage}

{all-args}{after-help}
";

#[derive(Parser, Debug)]
#[command(
    author = "Tony Kan, Ted Yu, William A. Goddard III, Victor Wai Tak Kam",
    version,
    about = "SCREAM++ CLI - A command-line interface for SCREAM++, an enhanced software package for automated protein side-chain placement and redesign.",
    help_template = HELP_TEMPLATE,
)]
#[command(propagate_version = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,

    /// Increase verbosity level (-v for INFO, -vv for DEBUG, -vvv for TRACE)
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    pub verbose: u8,

    /// Suppress all log output except for errors
    #[arg(short, long, global = true, conflicts_with = "verbose")]
    pub quiet: bool,

    /// Write logs to a specified file in addition to the console output
    #[arg(long, global = true, value_name = "PATH")]
    pub log_file: Option<PathBuf>,

    /// Set the number of threads for parallel computation.
    /// Defaults to the number of available logical cores.
    #[arg(short = 'j', long, global = true, value_name = "NUM")]
    pub threads: Option<usize>,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Place or redesign side-chains on a protein structure using the SCREAM algorithm.
    Place(PlaceArgs),
    /// Manage local data files (forcefields, rotamer libraries) required by SCREAM++.
    Data(DataArgs),
}

/// Arguments for the `place` subcommand.
#[derive(Args, Debug)]
pub struct PlaceArgs {
    // --- Core Arguments ---
    /// Path to the input molecular structure file (e.g., protein.bgf).
    #[arg(short, long, required = true, value_name = "PATH")]
    pub input: PathBuf,

    /// Path for the output molecular structure file.
    #[arg(short, long, required = true, value_name = "PATH")]
    pub output: PathBuf,

    /// Path to the main configuration file in TOML format.
    #[arg(short, long, required = true, value_name = "PATH")]
    pub config: PathBuf,

    // --- Forcefield Overrides ---
    /// Override the s-factor for the flat-bottom potential from the config file.
    #[arg(short = 's', long, value_name = "FLOAT")]
    pub s_factor: Option<f64>,

    /// Override the forcefield parameter file.
    /// Can be a path or a logical name (e.g., 'lj-12-6@0.4').
    #[arg(long, value_name = "NAME_OR_PATH")]
    pub forcefield_path: Option<String>,

    /// Override the delta parameter file for the flat-bottom potential.
    /// Can be a path or a logical name (e.g., 'rmsd-1.0').
    #[arg(long, value_name = "NAME_OR_PATH")]
    pub delta_params_path: Option<String>,

    // --- Sampling Overrides ---
    /// Override the rotamer library.
    /// Can be a path or a logical name (e.g., 'charmm@rmsd-1.0').
    #[arg(short = 'l', long = "rotamer-library", value_name = "NAME_OR_PATH")]
    pub rotamer_library: Option<String>,

    /// Override the rotamer placement registry file.
    /// Can be a path or a logical name (e.g., 'default').
    #[arg(short = 'p', long = "placement-registry", value_name = "NAME_OR_PATH")]
    pub placement_registry: Option<String>,

    // --- Optimization Overrides ---
    /// Override the number of solutions to generate.
    #[arg(short, long, value_name = "INT")]
    pub num_solutions: Option<usize>,

    /// Override the maximum number of clash resolution iterations.
    #[arg(long, value_name = "INT")]
    pub max_iterations: Option<usize>,

    /// Override `optimization.include-input-conformation` from the config file.
    #[command(flatten)]
    pub include_input_conformation: IncludeInputConformation,

    /// Disable the final refinement stage, overriding the config file.
    #[arg(long)]
    pub no_refinement: bool,

    /// Disable simulated annealing, even if it is defined in the config file.
    #[arg(long)]
    pub no_annealing: bool,

    /// Set a specific configuration value, overriding the config file.
    /// Can be used multiple times. Example: -S optimization.num_solutions=5
    #[arg(short = 'S', long = "set", value_name = "KEY=VALUE", num_args(0..))]
    pub set_values: Vec<String>,
}

/// A group to handle mutually exclusive boolean flags for including the input conformation.
#[derive(Args, Debug, Clone, Copy)]
#[group(required = false, multiple = false)]
pub struct IncludeInputConformation {
    /// Force inclusion of the input side-chain conformation in sampling.
    #[arg(long)]
    pub with_input_conformation: bool,
    /// Force exclusion of the input side-chain conformation from sampling.
    #[arg(long)]
    pub no_input_conformation: bool,
}

/// Arguments for the `data` subcommand.
#[derive(Args, Debug)]
pub struct DataArgs {
    #[command(subcommand)]
    pub command: DataCommands,
}

/// Available commands for data management.
#[derive(Subcommand, Debug)]
pub enum DataCommands {
    /// Download and unpack the default data files from the official repository.
    Download {
        /// Force re-download and overwrite existing data.
        #[arg(long)]
        force: bool,
    },
    /// Show the absolute path to the local data directory.
    Path,
    /// Set a custom absolute path for the local data directory.
    SetPath {
        /// The new path to use for storing data files.
        #[arg(required = true)]
        path: PathBuf,
    },
    /// Reset the data path to its default, OS-specific location.
    ResetPath,
}
