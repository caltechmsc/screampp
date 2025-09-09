use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;

/// Represents the topological structure of an amino acid residue.
///
/// This struct defines the atom composition and classification for a specific
/// residue type, distinguishing between anchor (backbone) atoms that define
/// the peptide chain connectivity and sidechain atoms that contribute to
/// the residue's unique chemical properties.
#[derive(Debug, Deserialize, Clone, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ResidueTopology {
    /// Names of atoms that serve as anchor points for the peptide backbone.
    ///
    /// These atoms (typically N, CA, C) define the connectivity and geometry
    /// of the polypeptide chain and are conserved across most amino acids.
    pub anchor_atoms: Vec<String>,
    /// Names of atoms that belong to the residue's sidechain.
    ///
    /// These atoms vary between different amino acid types and determine
    /// the chemical properties and functionality of each residue.
    pub sidechain_atoms: Vec<String>,
}

/// Manages a collection of residue topology definitions for molecular systems.
///
/// This registry provides centralized access to topological information for different
/// amino acid residue types, enabling consistent atom classification and structural
/// analysis across the molecular modeling pipeline.
#[derive(Debug, Clone, Default)]
pub struct TopologyRegistry {
    /// Internal storage mapping residue names to their topology definitions.
    registry: HashMap<String, ResidueTopology>,
}

impl TopologyRegistry {
    /// Loads residue topology definitions from a TOML configuration file.
    ///
    /// This method reads and parses a TOML file containing topology definitions
    /// for various amino acid residues, populating the registry with the parsed data.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the TOML file containing topology definitions
    ///
    /// # Return
    ///
    /// Returns a new `TopologyRegistry` instance populated with the loaded topologies,
    /// or an error if the file cannot be read or parsed.
    ///
    /// # Errors
    ///
    /// Returns `TopologyLoadError::Io` if the file cannot be read.
    /// Returns `TopologyLoadError::Toml` if the file content is not valid TOML or
    /// contains invalid topology definitions.
    pub fn load(path: &Path) -> Result<Self, TopologyLoadError> {
        // Read the entire file content into memory for TOML parsing
        let content = std::fs::read_to_string(path).map_err(|e| TopologyLoadError::Io {
            path: path.to_string_lossy().to_string(),
            source: e,
        })?;
        // Parse the TOML content into a HashMap of residue topologies
        let registry: HashMap<String, ResidueTopology> =
            toml::from_str(&content).map_err(|e| TopologyLoadError::Toml {
                path: path.to_string_lossy().to_string(),
                source: e,
            })?;
        Ok(Self { registry })
    }

    /// Retrieves the topology definition for a specific residue type.
    ///
    /// This method provides access to the topology information for a given
    /// amino acid residue, allowing classification of its constituent atoms.
    ///
    /// # Arguments
    ///
    /// * `residue_name` - Three-letter code of the amino acid residue (e.g., "ALA", "GLY")
    ///
    /// # Return
    ///
    /// Returns `Some(&ResidueTopology)` if the residue is found in the registry,
    /// or `None` if the residue type is not defined.
    pub fn get(&self, residue_name: &str) -> Option<&ResidueTopology> {
        self.registry.get(residue_name)
    }
}

/// Represents errors that can occur when loading topology definitions from files.
///
/// This enum encapsulates various failure modes during the topology loading process,
/// providing detailed context about what went wrong and where.
#[derive(Debug, Error)]
pub enum TopologyLoadError {
    /// Indicates that the topology file could not be read from disk.
    ///
    /// This error occurs when there are permission issues, the file doesn't exist,
    /// or other I/O-related problems prevent reading the topology configuration.
    #[error("File I/O error for '{path}': {source}")]
    Io {
        /// The path to the file that could not be read.
        path: String,
        /// The underlying I/O error that occurred.
        source: std::io::Error,
    },
    /// Indicates that the topology file content is not valid TOML or contains invalid data.
    ///
    /// This error occurs when the file exists but cannot be parsed as TOML,
    /// or when the parsed data doesn't match the expected `ResidueTopology` structure.
    #[error("TOML parsing error for '{path}': {source}")]
    Toml {
        /// The path to the file that could not be parsed.
        path: String,
        /// The underlying TOML parsing error that occurred.
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
