use crate::core::models::residue::ResidueType;
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct InputFiles {
    pub structure: PathBuf,
    pub non_bonded_params: PathBuf,
    pub topology_params: PathBuf,
    pub placement_params: PathBuf,
    pub charge_params: PathBuf,
    pub delta_params: PathBuf,
    pub rotamer_library: PathBuf,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TargetResidues {
    All,
    AllExcept(Vec<(char, isize)>),
    Explicit(Vec<(char, isize)>),
}
