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
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Place(PlaceArgs),
    Data(DataArgs),
}

#[derive(Args, Debug)]
struct PlaceArgs {
    #[arg(short, long)]
    input: PathBuf,

    #[arg(short, long)]
    output: PathBuf,
}

#[derive(Args, Debug)]
struct DataArgs {
    #[command(subcommand)]
    command: DataCommands,
}

#[derive(Subcommand, Debug)]
enum DataCommands {
    Download,
    Path,
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Place(args) => {
            println!("'place' command is not yet implemented.");
            println!("Input path: {:?}", args.input);
            println!("Output path: {:?}", args.output);
        }
        Commands::Data(args) => match args.command {
            DataCommands::Download => {
                println!("'data download' command is not yet implemented.");
            }
            DataCommands::Path => {
                println!("'data path' command is not yet implemented.");
            }
        },
    }
}
