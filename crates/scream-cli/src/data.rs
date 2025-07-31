use crate::error::{CliError, Result};
use crate::utils::parser::{ForcefieldName, ParsedLogicalName, RotamerLibraryName};
use directories::ProjectDirs;
use futures_util::StreamExt;
use std::fs::{self};
use std::path::{Path, PathBuf};
use tracing::{debug, info, warn};

const DATA_URL: &str = "https://github.com/caltechmsc/screampp/releases/download/XXX.tar.zst"; // TODO: Replace with actual URL

#[derive(Debug, Clone, Copy)]
pub enum DataProgress {
    DownloadStarted { total_size: Option<u64> },
    Downloading { downloaded: u64 },
    Unpacking,
}

#[derive(Debug)]
pub struct DataManager {
    base_path: PathBuf,
}
