use screampp::engine::config as core_config;
use std::path::PathBuf;

pub struct AppConfig {
    pub input_path: PathBuf,
    pub output_template: PathBuf,
    pub core_config: core_config::PlacementConfig,
}
