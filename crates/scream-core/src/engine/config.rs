use crate::core::models::residue::ResidueType;
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct InputPaths {
    pub structure_file: PathBuf,
    pub forcefield_file: PathBuf,
    pub topology_file: PathBuf,
    pub placement_file: PathBuf,
    pub rotamer_library_file: PathBuf,
    pub charge_file: PathBuf,
    pub delta_file: PathBuf,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PlacementSelection {
    All,
    AllExcept(Vec<(char, isize)>),
    Explicit(Vec<(char, isize)>),
}

pub type DesignMap = HashMap<(char, isize), Vec<ResidueType>>;

#[derive(Debug, Clone)]
pub enum PlacementMode {
    Standard(PlacementSelection),
    Design(DesignMap),
}

#[derive(Debug, Clone)]
pub struct EngineConfig {
    pub mode: PlacementMode,
    pub s_factor: f64,
    pub max_iterations: usize,
    pub convergence_tolerance: f64,
}
