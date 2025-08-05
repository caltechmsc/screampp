use screampp::engine::error::EngineError;
use std::path::PathBuf;
use thiserror::Error;

pub type Result<T> = std::result::Result<T, CliError>;

#[derive(Debug, Error)]
pub enum CliError {
    #[error(transparent)]
    ScreamCore(#[from] EngineError),

    #[error("Configuration error: {0}")]
    Config(String),

    #[error("Data management error: {0}")]
    Data(String),

    #[error("Failed to parse file '{path}': {source}", path = path.display())]
    FileParsing {
        path: PathBuf,
        #[source]
        source: anyhow::Error,
    },

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Network error: {0}")]
    Network(#[from] reqwest::Error),

    #[error("Invalid argument: {0}")]
    Argument(String),

    #[error(transparent)]
    Other(#[from] anyhow::Error),
}
