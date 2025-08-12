use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

const AUTHORS: &str = "Tony Kan, Ted Yu, Soo-Kyung Kim, William A. Goddard III, Victor Wai Tak Kam";
const ABOUT: &str = "SCREAM++ CLI - A command-line interface for SCREAM++, an advanced computational framework for automated protein side-chain placement and structural redesign.";
const COPYRIGHT: &str = "Copyright (c) 2025 California Institute of Technology, Materials and Process Simulation Center (MSC)";
const HELP_TEMPLATE: &str = "\
{before-help}{name} {version}
{author-with-newline}{about-with-newline}
{usage-heading} {usage}

{all-args}{after-help}
";

#[derive(Parser, Debug)]
#[command(
    author = AUTHORS,
    version,
    about = ABOUT,
    after_help = COPYRIGHT,
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

    /// Path for the output molecular structure file(s).
    ///
    /// This path acts as a template for naming output files, especially when generating
    /// multiple solutions (e.g., using --num-solutions). The following placeholders
    /// can be used in the filename to create dynamic names:
    ///
    ///   {n}, {i}        - The rank/index of the solution (1-based, 1 is the best).
    ///   {N}, {total}    - The total number of solutions being written.
    ///   {total_energy}  - The final total system energy of the solution.
    ///   {energy}        - An alias for {total_energy}.
    ///   {opt_score}     - The internal optimization score of the solution.
    ///   {score}         - An alias for {opt_score}.
    ///
    /// All energy/score values are formatted to two decimal places.
    ///
    /// Example Usage:
    ///   -o "solution_{i}.bgf"          -> solution_1.bgf, solution_2.bgf, ...
    ///   -o "run_A_sol_{i}_of_{N}.bgf"  -> run_A_sol_1_of_3.bgf, ...
    ///   -o "result_E_{energy}.bgf"     -> result_E_-15343.92.bgf, ...
    ///
    /// If multiple solutions are generated and no placeholder is found in the path, a
    /// "-best-{n}" suffix will be automatically added before the file extension to
    /// prevent overwriting files. For a single solution, the path is used as-is.
    #[arg(short, long, required = true, value_name = "PATH_TEMPLATE")]
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

    /// Override the residue topology registry file.
    /// Can be a path or a logical name (e.g., 'default').
    #[arg(short = 't', long = "topology-registry", value_name = "NAME_OR_PATH")]
    pub topology_registry: Option<String>,

    // --- Sampling Overrides ---
    /// Override the rotamer library.
    /// Can be a path or a logical name (e.g., 'charmm@rmsd-1.0').
    #[arg(short = 'l', long = "rotamer-library", value_name = "NAME_OR_PATH")]
    pub rotamer_library: Option<String>,

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verify_cli_parsing() {
        use clap::CommandFactory;
        Cli::command().debug_assert();
    }

    #[test]
    fn test_place_args_parsing() {
        let args = [
            "scream",
            "place",
            "-i",
            "in.bgf",
            "-o",
            "out.bgf",
            "-c",
            "config.toml",
            "-s",
            "0.8",
            "-l",
            "amber@rmsd-1.2",
            "-n",
            "5",
            "--with-input-conformation",
            "-S",
            "optimization.max-iterations=200",
        ];
        let cli = Cli::parse_from(args);
        match cli.command {
            Commands::Place(place_args) => {
                assert_eq!(place_args.input, PathBuf::from("in.bgf"));
                assert_eq!(place_args.output, PathBuf::from("out.bgf"));
                assert_eq!(place_args.config, PathBuf::from("config.toml"));
                assert_eq!(place_args.s_factor, Some(0.8));
                assert_eq!(
                    place_args.rotamer_library,
                    Some("amber@rmsd-1.2".to_string())
                );
                assert_eq!(place_args.num_solutions, Some(5));
                assert!(
                    place_args
                        .include_input_conformation
                        .with_input_conformation
                );
                assert!(!place_args.include_input_conformation.no_input_conformation);
                assert_eq!(
                    place_args.set_values,
                    vec!["optimization.max-iterations=200".to_string()]
                );
            }
            _ => panic!("Expected Place subcommand"),
        }
    }

    #[test]
    fn test_mutually_exclusive_flags() {
        let args = [
            "scream",
            "place",
            "-i",
            "in.bgf",
            "-o",
            "out.bgf",
            "-c",
            "c.toml",
            "--with-input-conformation",
            "--no-input-conformation",
        ];
        let result = Cli::try_parse_from(args);
        assert!(result.is_err(), "Clap should reject conflicting flags");
    }

    #[test]
    fn test_data_download_parsing() {
        let args = ["scream", "data", "download", "--force"];
        let cli = Cli::parse_from(args);
        match cli.command {
            Commands::Data(data_args) => match data_args.command {
                DataCommands::Download { force } => {
                    assert!(force);
                }
                _ => panic!("Expected data download subcommand"),
            },
            _ => panic!("Expected data subcommand"),
        }
    }
}
