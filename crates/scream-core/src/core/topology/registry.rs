use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Deserialize, Clone, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ResidueTopology {
    pub anchor_atoms: Vec<String>,
    pub sidechain_atoms: Vec<String>,
}

#[derive(Debug, Clone, Default)]
pub struct TopologyRegistry {
    registry: HashMap<String, ResidueTopology>,
}

impl TopologyRegistry {
    pub fn load(path: &Path) -> Result<Self, TopologyLoadError> {
        let content = std::fs::read_to_string(path).map_err(|e| TopologyLoadError::Io {
            path: path.to_string_lossy().to_string(),
            source: e,
        })?;
        let registry: HashMap<String, ResidueTopology> =
            toml::from_str(&content).map_err(|e| TopologyLoadError::Toml {
                path: path.to_string_lossy().to_string(),
                source: e,
            })?;
        Ok(Self { registry })
    }

    pub fn get(&self, residue_name: &str) -> Option<&ResidueTopology> {
        self.registry.get(residue_name)
    }
}

#[derive(Debug, Error)]
pub enum TopologyLoadError {
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_registry_file(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        write!(file, "{}", content).unwrap();
        file
    }

    #[test]
    fn loads_registry_successfully_and_get_works() {
        let content = r#"
[ALA]
anchor_atoms = ["N", "CA", "C"]
sidechain_atoms = ["CB"]

[GLY]
anchor_atoms = ["N", "CA", "C"]
sidechain_atoms = ["HA"]
"#;
        let file = create_test_registry_file(content);

        let registry = TopologyRegistry::load(file.path()).unwrap();

        assert_eq!(registry.registry.len(), 2);

        let ala_topo = registry.get("ALA").expect("ALA topology should be present");
        assert_eq!(ala_topo.anchor_atoms, vec!["N", "CA", "C"]);
        assert_eq!(ala_topo.sidechain_atoms, vec!["CB"]);

        assert!(registry.get("LEU").is_none());
    }

    #[test]
    fn loads_empty_registry_from_empty_file() {
        let file = create_test_registry_file("");
        let registry = TopologyRegistry::load(file.path()).unwrap();
        assert!(registry.registry.is_empty());
    }

    #[test]
    fn load_returns_io_error_for_nonexistent_file() {
        let path = Path::new("nonexistent_topology_file.toml");
        let result = TopologyRegistry::load(path);
        assert!(matches!(result, Err(TopologyLoadError::Io { .. })));
    }

    #[test]
    fn load_returns_toml_error_for_malformed_file() {
        let file = create_test_registry_file("this is not valid toml");
        let result = TopologyRegistry::load(file.path());
        assert!(matches!(result, Err(TopologyLoadError::Toml { .. })));
    }
}
