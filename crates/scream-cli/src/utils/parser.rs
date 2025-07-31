use thiserror::Error;

#[derive(Debug, Error, PartialEq, Eq)]
pub enum ParseError {
    #[error(
        "Invalid rotamer library format for '{0}'. Expected 'scheme@diversity' (e.g., 'charmm@rmsd-1.0')."
    )]
    InvalidRotamerLibraryFormat(String),

    #[error("Invalid forcefield format for '{0}'. Expected 'type@version' (e.g., 'lj-12-6@0.3').")]
    InvalidForcefieldFormat(String),

    #[error("Invalid placement registry name '{0}'. Only 'default' is recognized.")]
    InvalidPlacementRegistryName(String),

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
    PlacementRegistry,
}

pub fn parse_logical_name(name: &str, kind: &str) -> Result<ParsedLogicalName, ParseError> {
    match kind {
        "rotamer-library" => parse_rotamer_library(name).map(ParsedLogicalName::RotamerLibrary),
        "forcefield" => parse_forcefield(name).map(ParsedLogicalName::Forcefield),
        "delta-params" => Ok(ParsedLogicalName::DeltaParams {
            diversity: name.to_string(),
        }),
        "placement-registry" => {
            if name == "default" {
                Ok(ParsedLogicalName::PlacementRegistry)
            } else {
                Err(ParseError::InvalidPlacementRegistryName(name.to_string()))
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
