mod cli;
mod commands;
mod config;
mod data;
mod error;
mod logging;
mod ui;
mod utils;

use crate::cli::{Cli, Commands};
use crate::error::{CliError, Result};
use crate::ui::UiManager;
use clap::Parser;
use tokio::task;
use tracing::{debug, error, info};

#[tokio::main]
async fn main() {
    if let Err(e) = run().await {
        tokio::time::sleep(std::time::Duration::from_millis(100)).await;
        eprintln!("\nError: {}", e);
        std::process::exit(1);
    }
}

async fn run() -> Result<()> {
    let cli = Cli::parse();

    let (ui_manager, ui_sender) = UiManager::new();

    let ui_handle = task::spawn(ui_manager.run());

    logging::setup_logging(cli.verbose, cli.quiet, &cli.log_file, ui_sender.clone())?;

    let (panic_hook, eyre_hook) = color_eyre::config::HookBuilder::default().into_hooks();
    eyre_hook.install().unwrap();
    std::panic::set_hook(Box::new(move |pi| {
        error!("{}", panic_hook.panic_report(pi));
    }));

    info!(
        "🚀 SCREAM++ CLI v{} starting up.",
        env!("CARGO_PKG_VERSION")
    );
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

    let command_result = match cli.command {
        Commands::Place(args) => {
            info!("Dispatching to 'place' command.");
            commands::place::run(args, ui_sender.clone()).await
        }
        Commands::Data(args) => {
            info!("Dispatching to 'data' command.");
            commands::data::run(args).await
        }
    };

    drop(ui_sender);
    ui_handle
        .await
        .map_err(|e| CliError::Other(anyhow::anyhow!("UI manager task failed: {}", e)))?;

    if command_result.is_ok() {
        info!("🎉 Command executed successfully. Shutting down.");
    }

    command_result
}
