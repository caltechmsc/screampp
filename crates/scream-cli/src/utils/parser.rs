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
