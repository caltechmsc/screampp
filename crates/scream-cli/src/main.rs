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
use tracing::{debug, error, info, warn};

#[tokio::main]
async fn main() {
    if let Err(e) = run_app().await {
        tokio::time::sleep(std::time::Duration::from_millis(50)).await;
        eprintln!("\n❌ Error: {}", e);
        std::process::exit(1);
    }
}

async fn run_app() -> Result<()> {
    let (ui_manager, ui_sender, shutdown_sender) = UiManager::new();
    let ui_handle = task::spawn(ui_manager.run());

    let cli = Cli::parse();
    logging::setup_logging(cli.verbose, cli.quiet, &cli.log_file, ui_sender.clone())?;

    let (panic_hook, eyre_hook) = color_eyre::config::HookBuilder::default().into_hooks();
    eyre_hook.install().map_err(|e| CliError::Other(e.into()))?;
    std::panic::set_hook(Box::new(move |pi| {
        error!("{}", panic_hook.panic_report(pi));
    }));

    let command_result = async {
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

        match cli.command {
            Commands::Place(args) => {
                info!("Dispatching to 'place' command.");
                commands::place::run(args, ui_sender).await
            }
            Commands::Data(args) => {
                info!("Dispatching to 'data' command.");
                commands::data::run(args).await
            }
        }
    }
    .await;

    match &command_result {
        Ok(_) => {
            info!("✅ Command completed successfully.");
            println!("✅ Command completed successfully.");
        }
        Err(e) => {
            error!("❌ Command failed: {}", e);
            eprintln!("❌ Command failed: {}", e);
        }
    }

    if shutdown_sender.send(true).is_err() {
        warn!("UI manager may have already exited before shutdown signal.");
    }

    ui_handle
        .await
        .map_err(|e| CliError::Other(anyhow::anyhow!("UI manager task failed: {}", e)))?;

    command_result
}
