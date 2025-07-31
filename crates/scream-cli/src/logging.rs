use crate::error::{CliError, Result};
use std::fs::File;
use std::path::PathBuf;
use tracing_subscriber::{
    filter::LevelFilter,
    fmt::{self},
    prelude::*,
};

pub fn setup_logging(verbosity: u8, quiet: bool, log_file: Option<PathBuf>) -> Result<()> {
    let level_filter = if quiet {
        LevelFilter::OFF
    } else {
        match verbosity {
            0 => LevelFilter::WARN,
            1 => LevelFilter::INFO,
            2 => LevelFilter::DEBUG,
            _ => LevelFilter::TRACE,
        }
    };

    let stderr_layer = fmt::layer()
        .with_writer(std::io::stderr)
        .with_ansi(true)
        .with_target(false)
        .compact();

    let subscriber = tracing_subscriber::registry()
        .with(level_filter)
        .with(stderr_layer);

    if let Some(path) = log_file {
        let file = File::create(&path).map_err(|e| CliError::Io(e))?;

        let file_layer = fmt::layer()
            .with_writer(file)
            .with_ansi(false)
            .with_thread_ids(true)
            .with_target(true);

        subscriber.with(file_layer).init();
    } else {
        subscriber.init();
    }

    Ok(())
}
