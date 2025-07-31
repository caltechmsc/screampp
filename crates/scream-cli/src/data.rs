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

impl DataManager {
    pub fn new() -> Result<Self> {
        let path = Self::determine_data_path()?;
        debug!("DataManager initialized with path: {:?}", &path);
        Ok(Self { base_path: path })
    }

    pub fn get_data_path(&self) -> &Path {
        &self.base_path
    }

    pub async fn download_data(
        &self,
        force: bool,
        mut progress_callback: impl FnMut(DataProgress),
    ) -> Result<()> {
        info!("Preparing to download data to {:?}", &self.base_path);
        if self.base_path.exists() {
            if force {
                info!("--force specified, removing existing data directory.");
                fs::remove_dir_all(&self.base_path)?;
            } else {
                return Err(CliError::Data(
                    "Data directory already exists. Use --force to overwrite.".to_string(),
                ));
            }
        }
        fs::create_dir_all(&self.base_path)?;

        info!("Sending request to {}", DATA_URL);
        let client = reqwest::Client::new();
        let response = client.get(DATA_URL).send().await?.error_for_status()?;

        let total_size = response.content_length();
        progress_callback(DataProgress::DownloadStarted { total_size });

        let mut downloaded: u64 = 0;
        let mut stream = response.bytes_stream();
        let mut buffer: Vec<u8> = Vec::with_capacity(total_size.unwrap_or(0) as usize);

        while let Some(item) = stream.next().await {
            let chunk = item?;
            buffer.extend_from_slice(&chunk);
            downloaded += chunk.len() as u64;
            progress_callback(DataProgress::Downloading { downloaded });
        }

        progress_callback(DataProgress::Unpacking);
        info!("Download complete. Decompressing and unpacking archive...");

        let cursor = std::io::Cursor::new(buffer);
        let zstd_decoder = zstd::stream::read::Decoder::new(cursor)?;
        let mut archive = tar::Archive::new(zstd_decoder);
        archive
            .unpack(&self.base_path)
            .map_err(|e| CliError::Io(e))?;

        info!("Data successfully unpacked to {:?}", &self.base_path);
        Ok(())
    }

    pub fn set_custom_path(path: &Path) -> Result<()> {
        let config_path = Self::get_path_config_file()?;
        if let Some(parent) = config_path.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::write(config_path, path.to_str().unwrap()).map_err(CliError::from)
    }

    pub fn reset_path() -> Result<()> {
        if let Ok(config_path) = Self::get_path_config_file() {
            if config_path.exists() {
                fs::remove_file(config_path)?;
            }
        }
        Ok(())
    }

    pub fn resolve_logical_name(&self, parsed_name: &ParsedLogicalName) -> Result<PathBuf> {
        let resolved_path = match parsed_name {
            ParsedLogicalName::RotamerLibrary(RotamerLibraryName { scheme, diversity }) => self
                .base_path
                .join("rotamers")
                .join(scheme)
                .join(format!("{}.toml", diversity)),

            ParsedLogicalName::Forcefield(ForcefieldName {
                potential_type,
                version,
            }) => self
                .base_path
                .join("forcefield")
                .join(format!("dreiding-{}-{}.toml", potential_type, version)),

            ParsedLogicalName::DeltaParams { diversity } => self
                .base_path
                .join("delta")
                .join(format!("delta-{}.csv", diversity)),

            ParsedLogicalName::PlacementRegistry => {
                self.base_path.join("rotamers").join("placement.toml")
            }
        };

        Ok(resolved_path)
    }

    fn determine_data_path() -> Result<PathBuf> {
        match Self::get_path_config_file() {
            Ok(config_path) if config_path.exists() => {
                let custom_path_str = fs::read_to_string(&config_path)?.trim().to_string();
                if custom_path_str.is_empty() {
                    warn!("Custom path config file is empty, falling back to default path.");
                    Self::get_default_data_path()
                } else {
                    Ok(PathBuf::from(custom_path_str))
                }
            }
            _ => Self::get_default_data_path(),
        }
    }

    fn get_path_config_file() -> Result<PathBuf> {
        ProjectDirs::from("edu", "caltech", "screampp")
            .map(|dirs| dirs.config_dir().join("path.conf"))
            .ok_or_else(|| CliError::Data("Could not determine config directory path.".to_string()))
    }

    fn get_default_data_path() -> Result<PathBuf> {
        ProjectDirs::from("edu", "caltech", "screampp")
            .map(|dirs| dirs.data_dir().to_path_buf())
            .ok_or_else(|| {
                CliError::Data("Could not determine default data directory path.".to_string())
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::parser;
    use tempfile::tempdir;

    #[test]
    fn resolve_logical_names_constructs_correct_paths() {
        let temp_dir = tempdir().unwrap();
        let base_path = temp_dir.path();

        fs::create_dir_all(base_path.join("rotamers/amber")).unwrap();
        fs::write(base_path.join("rotamers/amber/rmsd-1.0.toml"), "").unwrap();
        fs::create_dir_all(base_path.join("delta")).unwrap();
        fs::write(base_path.join("delta/delta-all-torsion.csv"), "").unwrap();
        fs::create_dir_all(base_path.join("forcefield")).unwrap();
        fs::write(base_path.join("forcefield/dreiding-lj-12-6-0.4.toml"), "").unwrap();
        fs::write(base_path.join("rotamers/placement.toml"), "").unwrap();

        let manager = DataManager {
            base_path: base_path.to_path_buf(),
        };

        let rot_parsed = parser::parse_logical_name("amber@rmsd-1.0", "rotamer-library").unwrap();
        let rot_path = manager.resolve_logical_name(&rot_parsed).unwrap();
        assert_eq!(rot_path, base_path.join("rotamers/amber/rmsd-1.0.toml"));

        let delta_parsed = parser::parse_logical_name("all-torsion", "delta-params").unwrap();
        let delta_path = manager.resolve_logical_name(&delta_parsed).unwrap();
        assert_eq!(delta_path, base_path.join("delta/delta-all-torsion.csv"));

        let ff_parsed = parser::parse_logical_name("lj-12-6@0.4", "forcefield").unwrap();
        let ff_path = manager.resolve_logical_name(&ff_parsed).unwrap();
        assert_eq!(
            ff_path,
            base_path.join("forcefield/dreiding-lj-12-6-0.4.toml")
        );

        let reg_parsed = parser::parse_logical_name("default", "placement-registry").unwrap();
        let reg_path = manager.resolve_logical_name(&reg_parsed).unwrap();
        assert_eq!(reg_path, base_path.join("rotamers/placement.toml"));
    }

    #[test]
    fn resolve_logical_name_fails_if_file_does_not_exist() {
        let temp_dir = tempdir().unwrap();
        let manager = DataManager {
            base_path: temp_dir.path().to_path_buf(),
        };

        let parsed = parser::parse_logical_name("amber@rmsd-1.0", "rotamer-library").unwrap();
        let result = manager.resolve_logical_name(&parsed);

        assert!(result.is_ok());

        let resolved_path = result.unwrap();

        assert!(!resolved_path.exists());
        assert!(
            resolved_path
                .to_string_lossy()
                .contains("rotamers/amber/rmsd-1.0.toml")
        );
    }
}
