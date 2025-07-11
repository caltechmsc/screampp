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

#[derive(Debug, Clone)]
pub struct DesignTask {
    pub positions: HashMap<(char, isize), Vec<ResidueType>>,
}

#[derive(Debug, Clone)]
pub struct InteractionAnalysisTask {
    pub group1_selector: String,
    pub group2_selector: String,
}

#[derive(Debug, Clone)]
pub enum ScreamTask {
    Place(TargetResidues),
    Design(DesignTask),
    Analyze(InteractionAnalysisTask),
}

#[derive(Debug, Clone, Copy)]
pub struct AlgorithmConfig {
    pub s_factor: f64,
    pub max_iterations: usize,
    pub convergence_tolerance: f64,
}

impl Default for AlgorithmConfig {
    fn default() -> Self {
        Self {
            s_factor: 1.0,
            max_iterations: 10,
            convergence_tolerance: 0.01,
        }
    }
}
