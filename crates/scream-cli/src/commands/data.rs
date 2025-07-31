use crate::cli::{DataArgs, DataCommands};
use crate::data::{DataManager, DataProgress};
use crate::error::Result;
use indicatif::{ProgressBar, ProgressStyle};
use tracing::info;

pub async fn run(args: DataArgs) -> Result<()> {
    match args.command {
        DataCommands::Download { force } => {
            handle_download(force).await?;
        }
        DataCommands::Path => {
            handle_path()?;
        }
        DataCommands::SetPath { path } => {
            handle_set_path(path)?;
        }
        DataCommands::ResetPath => {
            handle_reset_path()?;
        }
    }
    Ok(())
}

async fn handle_download(force: bool) -> Result<()> {
    println!("Initializing data manager...");
    let manager = DataManager::new()?;

    let pb = ProgressBar::new(0);
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({bytes_per_sec}, {eta})",
        )
            .unwrap()
            .progress_chars("#>-"),
    );
    pb.set_draw_target(indicatif::ProgressDrawTarget::stderr_with_hz(2));

    println!("Downloading default data to: {:?}", manager.get_data_path());

    let progress_callback = |progress: DataProgress| match progress {
        DataProgress::DownloadStarted { total_size } => {
            if let Some(size) = total_size {
                pb.set_length(size);
            }
            pb.set_message("Downloading...");
        }
        DataProgress::Downloading { downloaded } => {
            pb.set_position(downloaded);
        }
        DataProgress::Unpacking => {
            pb.set_style(ProgressStyle::with_template("{spinner:.green} {msg}").unwrap());
            pb.set_message("Unpacking archive...");
        }
    };

    match manager.download_data(force, progress_callback).await {
        Ok(_) => {
            pb.finish_with_message("✓ Data download and setup complete.");
            Ok(())
        }
        Err(e) => {
            pb.finish_with_message("✗ Download failed.");
            Err(e)
        }
    }
}

fn handle_path() -> Result<()> {
    let manager = DataManager::new()?;
    println!(
        "Current SCREAM++ data directory: {}",
        manager.get_data_path().display()
    );
    Ok(())
}

fn handle_set_path(path: std::path::PathBuf) -> Result<()> {
    let absolute_path = if path.is_absolute() {
        path
    } else {
        std::env::current_dir()?.join(path)
    };

    info!("Setting custom data path to: {:?}", &absolute_path);
    DataManager::set_custom_path(&absolute_path)?;
    println!("Successfully set data path to: {}", absolute_path.display());
    println!("Future runs will now use this directory.");
    Ok(())
}

fn handle_reset_path() -> Result<()> {
    info!("Resetting data path to default.");
    DataManager::reset_path()?;
    let manager = DataManager::new()?;
    println!(
        "Data path has been reset to the default: {}",
        manager.get_data_path().display()
    );
    Ok(())
}
