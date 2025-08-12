use thiserror::Error;

#[derive(Debug, Error, PartialEq, Eq)]
pub enum ParseError {
    #[error(
        "Invalid rotamer library format for '{0}'. Expected 'scheme@diversity' (e.g., 'charmm@rmsd-1.0')."
    )]
    InvalidRotamerLibraryFormat(String),

    #[error("Invalid forcefield format for '{0}'. Expected 'type@version' (e.g., 'lj-12-6@0.3').")]
    InvalidForcefieldFormat(String),

    #[error("Invalid topology registry name '{0}'. Only 'default' is recognized.")]
    InvalidTopologyRegistryName(String),

    #[error("Unknown logical name kind: '{0}'. Expected 'rotamer-library', 'forcefield', etc.")]
    UnknownKind(String),

    #[error("Component '{component}' cannot be empty in logical name '{name}'.")]
    EmptyComponent {
        component: &'static str,
        name: String,
    },
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct RotamerLibraryName {
    pub scheme: String,
    pub diversity: String,
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct ForcefieldName {
    pub potential_type: String,
    pub version: String,
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub enum ParsedLogicalName {
    RotamerLibrary(RotamerLibraryName),
    Forcefield(ForcefieldName),
    DeltaParams { diversity: String },
    TopologyRegistry,
}

pub fn parse_logical_name(name: &str, kind: &str) -> Result<ParsedLogicalName, ParseError> {
    match kind {
        "rotamer-library" => parse_rotamer_library(name).map(ParsedLogicalName::RotamerLibrary),
        "forcefield" => parse_forcefield(name).map(ParsedLogicalName::Forcefield),
        "delta-params" => Ok(ParsedLogicalName::DeltaParams {
            diversity: name.to_string(),
        }),
        "topology-registry" => {
            if name == "default" {
                Ok(ParsedLogicalName::TopologyRegistry)
            } else {
                Err(ParseError::InvalidTopologyRegistryName(name.to_string()))
            }
        }
        _ => Err(ParseError::UnknownKind(kind.to_string())),
    }
}

fn parse_rotamer_library(name: &str) -> Result<RotamerLibraryName, ParseError> {
    let parts: Vec<&str> = name.splitn(2, '@').collect();

    if parts.len() != 2 {
        return Err(ParseError::InvalidRotamerLibraryFormat(name.to_string()));
    }

    let scheme = parts[0].trim().to_string();
    let diversity = parts[1].trim().to_string();

    if scheme.is_empty() {
        return Err(ParseError::EmptyComponent {
            component: "scheme",
            name: name.to_string(),
        });
    }
    if diversity.is_empty() {
        return Err(ParseError::EmptyComponent {
            component: "diversity",
            name: name.to_string(),
        });
    }

    Ok(RotamerLibraryName { scheme, diversity })
}

fn parse_forcefield(name: &str) -> Result<ForcefieldName, ParseError> {
    let parts: Vec<&str> = name.splitn(2, '@').collect();

    if parts.len() != 2 {
        return Err(ParseError::InvalidForcefieldFormat(name.to_string()));
    }

    let potential_type = parts[0].trim().to_string();
    let version = parts[1].trim().to_string();

    if potential_type.is_empty() {
        return Err(ParseError::EmptyComponent {
            component: "type",
            name: name.to_string(),
        });
    }
    if version.is_empty() {
        return Err(ParseError::EmptyComponent {
            component: "version",
            name: name.to_string(),
        });
    }

    Ok(ForcefieldName {
        potential_type,
        version,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_rotamer_library_success() {
        let result = parse_logical_name("amber@rmsd-1.0", "rotamer-library").unwrap();
        assert_eq!(
            result,
            ParsedLogicalName::RotamerLibrary(RotamerLibraryName {
                scheme: "amber".to_string(),
                diversity: "rmsd-1.0".to_string(),
            })
        );
    }

    #[test]
    fn test_parse_rotamer_library_with_spaces() {
        let result = parse_logical_name(" charmm @ all-torsion ", "rotamer-library").unwrap();
        assert_eq!(
            result,
            ParsedLogicalName::RotamerLibrary(RotamerLibraryName {
                scheme: "charmm".to_string(),
                diversity: "all-torsion".to_string(),
            })
        );
    }

    #[test]
    fn test_parse_rotamer_library_fails_no_at() {
        let result = parse_logical_name("charmm-rmsd-1.0", "rotamer-library");
        assert!(matches!(
            result,
            Err(ParseError::InvalidRotamerLibraryFormat(_))
        ));
    }

    #[test]
    fn test_parse_rotamer_library_fails_empty_scheme() {
        let result = parse_logical_name("@rmsd-1.0", "rotamer-library");
        assert_eq!(
            result,
            Err(ParseError::EmptyComponent {
                component: "scheme",
                name: "@rmsd-1.0".to_string()
            })
        );
    }

    #[test]
    fn test_parse_rotamer_library_fails_empty_diversity() {
        let result = parse_logical_name("charmm@", "rotamer-library");
        assert_eq!(
            result,
            Err(ParseError::EmptyComponent {
                component: "diversity",
                name: "charmm@".to_string()
            })
        );
    }

    #[test]
    fn test_parse_forcefield_success() {
        let result = parse_logical_name("lj-12-6@0.3", "forcefield").unwrap();
        assert_eq!(
            result,
            ParsedLogicalName::Forcefield(ForcefieldName {
                potential_type: "lj-12-6".to_string(),
                version: "0.3".to_string(),
            })
        );
    }

    #[test]
    fn test_parse_forcefield_fails_no_at() {
        let result = parse_logical_name("exp-6-0.4", "forcefield");
        assert!(matches!(
            result,
            Err(ParseError::InvalidForcefieldFormat(_))
        ));
    }

    #[test]
    fn test_parse_forcefield_fails_empty_version() {
        let result = parse_logical_name("lj-12-6@", "forcefield");
        assert_eq!(
            result,
            Err(ParseError::EmptyComponent {
                component: "version",
                name: "lj-12-6@".to_string()
            })
        );
    }

    #[test]
    fn test_parse_delta_params() {
        let result = parse_logical_name("rmsd-1.4", "delta-params").unwrap();
        assert_eq!(
            result,
            ParsedLogicalName::DeltaParams {
                diversity: "rmsd-1.4".to_string()
            }
        );
    }

    #[test]
    fn test_parse_topology_registry_success() {
        let result = parse_logical_name("default", "topology-registry").unwrap();
        assert_eq!(result, ParsedLogicalName::TopologyRegistry);
    }

    #[test]
    fn test_parse_topology_registry_fails_invalid_name() {
        let result = parse_logical_name("custom", "topology-registry");
        assert!(matches!(
            result,
            Err(ParseError::InvalidTopologyRegistryName(_))
        ));
    }

    #[test]
    fn test_parse_unknown_kind() {
        let result = parse_logical_name("some-value", "unknown-kind");
        assert!(matches!(result, Err(ParseError::UnknownKind(_))));
    }
}
