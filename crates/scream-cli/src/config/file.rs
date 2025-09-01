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

#[derive(Deserialize, Debug, Clone)]
#[serde(rename_all = "kebab-case", tag = "type")]
pub enum FileResidueSelection {
    All,
    List {
        #[serde(default)]
        include: Vec<FileResidueSpecifier>,
        #[serde(default)]
        exclude: Vec<FileResidueSpecifier>,
    },
    LigandBindingSite {
        #[serde(rename = "ligand-residue")]
        ligand_residue: FileResidueSpecifier,
        #[serde(rename = "radius-angstroms")]
        radius_angstroms: f64,
    },
}

impl From<FileResidueSelection> for core_config::ResidueSelection {
    fn from(p: FileResidueSelection) -> Self {
        match p {
            FileResidueSelection::All => core_config::ResidueSelection::All,
            FileResidueSelection::List { include, exclude } => {
                core_config::ResidueSelection::List {
                    include: include.into_iter().map(Into::into).collect(),
                    exclude: exclude.into_iter().map(Into::into).collect(),
                }
            }
            FileResidueSelection::LigandBindingSite {
                ligand_residue,
                radius_angstroms,
            } => core_config::ResidueSelection::LigandBindingSite {
                ligand_residue: ligand_residue.into(),
                radius_angstroms,
            },
        }
    }
}
