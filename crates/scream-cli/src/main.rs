mod cli;
mod commands;
mod config;
mod data;
mod error;
mod logging;
mod utils;

use crate::cli::{Cli, Commands};
use crate::error::{CliError, Result};
use clap::Parser;
use tracing::{debug, error, info};

#[tokio::main]
async fn main() {
    if let Err(e) = run().await {
        error!("An error occurred: {}", e);
        eprintln!("\nError: {}", e);
        std::process::exit(1);
    }
}

async fn run() -> Result<()> {
    let cli = Cli::parse();

    logging::setup_logging(cli.verbose, cli.quiet, &cli.log_file)?;

    info!("SCREAM++ CLI v{} starting up.", env!("CARGO_PKG_VERSION"));
    debug!("Full CLI arguments parsed: {:?}", &cli);

    if let Some(num_threads) = cli.threads {
        info!(
            "Setting Rayon global thread pool to {} threads.",
            num_threads
        );
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .map_err(|e| {
                CliError::Other(anyhow::anyhow!("Failed to build global thread pool: {}", e))
            })?;
    }

    match cli.command {
        Commands::Place(args) => {
            info!("Dispatching to 'place' command.");
            commands::place::run(args).await?;
        }
        Commands::Data(args) => {
            info!("Dispatching to 'data' command.");
            commands::data::run(args).await?;
        }
    }

    info!("Command executed successfully. Shutting down.");
    Ok(())
}
