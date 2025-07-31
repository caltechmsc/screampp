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
