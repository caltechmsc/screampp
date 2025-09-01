use crate::error::{CliError, Result};
use screampp::core::models::atom::AtomRole;
use screampp::engine::config as core_config;
use serde::Deserialize;
use std::path::Path;
use tracing::debug;

#[derive(Deserialize, Debug, Default, Clone)]
#[serde(deny_unknown_fields)]
pub struct FileResidueSpecifier {
    #[serde(rename = "chain-id")]
    pub chain_id: char,
    #[serde(rename = "residue-number")]
    pub residue_number: isize,
}

impl From<FileResidueSpecifier> for core_config::ResidueSpecifier {
    fn from(p: FileResidueSpecifier) -> Self {
        Self {
            chain_id: p.chain_id,
            residue_number: p.residue_number,
        }
    }
}
