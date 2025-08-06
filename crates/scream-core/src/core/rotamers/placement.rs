use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct PlacementInfo {
    pub anchor_atoms: Vec<String>,
    pub sidechain_atoms: Vec<String>,
}

pub type PlacementRegistry = HashMap<String, PlacementInfo>;

#[derive(Debug, Error)]
pub enum PlacementLoadError {
    #[error("File I/O error for '{path}': {source}")]
    Io {
        path: String,
        source: std::io::Error,
    },
    #[error("TOML parsing error for '{path}': {source}")]
    Toml {
        path: String,
        source: toml::de::Error,
    },
}

pub fn load_placement_registry(path: &Path) -> Result<PlacementRegistry, PlacementLoadError> {
    let content = std::fs::read_to_string(path).map_err(|e| PlacementLoadError::Io {
        path: path.to_string_lossy().to_string(),
        source: e,
    })?;
    toml::from_str(&content).map_err(|e| PlacementLoadError::Toml {
        path: path.to_string_lossy().to_string(),
        source: e,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn loads_registry_successfully_from_valid_file() {
        let content = r#"
[ALA]
anchor_atoms = ["N", "CA", "C"]
sidechain_atoms = ["CB", "HB1", "HB2", "HB3"]

[GLY]
anchor_atoms = ["N", "CA", "C"]
sidechain_atoms = ["HA1", "HA2"]
"#;
        let mut file = NamedTempFile::new().unwrap();
        write!(file, "{}", content).unwrap();

        let registry = load_placement_registry(file.path()).unwrap();

        assert_eq!(registry.len(), 2);

        let ala_info = registry.get("ALA").expect("ALA info should be present");
        assert_eq!(ala_info.anchor_atoms, vec!["N", "CA", "C"]);
        assert_eq!(ala_info.sidechain_atoms, vec!["CB", "HB1", "HB2", "HB3"]);

        let gly_info = registry.get("GLY").expect("GLY info should be present");
        assert_eq!(gly_info.anchor_atoms, vec!["N", "CA", "C"]);
        assert_eq!(gly_info.sidechain_atoms, vec!["HA1", "HA2"]);
    }

    #[test]
    fn loads_empty_registry_from_empty_file() {
        let mut file = NamedTempFile::new().unwrap();
        write!(file, "").unwrap();

        let registry = load_placement_registry(file.path()).unwrap();

        assert!(registry.is_empty());
    }

    #[test]
    fn returns_io_error_for_nonexistent_file() {
        let path = Path::new("nonexistent_placement_file.toml");

        let result = load_placement_registry(path);

        assert!(matches!(result, Err(PlacementLoadError::Io { .. })));
    }

    #[test]
    fn returns_toml_error_for_malformed_file() {
        let content = "this is not valid toml";
        let mut file = NamedTempFile::new().unwrap();
        write!(file, "{}", content).unwrap();

        let result = load_placement_registry(file.path());

        assert!(matches!(result, Err(PlacementLoadError::Toml { .. })));
    }

    #[test]
    fn returns_toml_error_for_incorrect_structure() {
        let content = r#"
[[ALA]]
anchor_atoms = ["N", "CA", "C"]
"#;
        let mut file = NamedTempFile::new().unwrap();
        write!(file, "{}", content).unwrap();

        let result = load_placement_registry(file.path());

        assert!(matches!(result, Err(PlacementLoadError::Toml { .. })));
    }
}
