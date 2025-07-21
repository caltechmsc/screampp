use crate::core::models::residue::ResidueType;
use itertools::Itertools;
use std::collections::HashMap;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error, PartialEq, Eq, Clone)]
pub enum ConfigError {
    #[error("Missing required parameter: {0}")]
    MissingParameter(&'static str),
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ResidueSpecifier {
    pub chain_id: char,
    pub residue_number: isize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ResidueSelection {
    All,
    List {
        include: Vec<ResidueSpecifier>,
        exclude: Vec<ResidueSpecifier>,
    },
    LigandBindingSite {
        ligand_residue: ResidueSpecifier,
        radius_angstroms: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct ForcefieldConfig {
    pub forcefield_path: PathBuf,
    pub delta_params_path: PathBuf,
    pub s_factor: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SamplingConfig {
    pub rotamer_library_path: PathBuf,
    pub placement_registry_path: PathBuf,
}

#[derive(Debug, Clone, PartialEq)]
pub struct OptimizationConfig {
    pub max_iterations: usize,
    pub convergence_threshold: f64,
    pub num_solutions: usize,
    pub include_input_conformation: bool,
    pub use_simulated_annealing: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PlacementConfig {
    pub forcefield: ForcefieldConfig,
    pub sampling: SamplingConfig,
    pub optimization: OptimizationConfig,
    pub residues_to_optimize: ResidueSelection,
}

pub type DesignSpec = HashMap<ResidueSpecifier, Vec<ResidueType>>;

pub trait DesignSpecExt {
    fn get_by_specifier(&self, chain_id: char, residue_number: isize) -> Option<&Vec<ResidueType>>;
}

impl DesignSpecExt for DesignSpec {
    fn get_by_specifier(&self, chain_id: char, residue_number: isize) -> Option<&Vec<ResidueType>> {
        self.get(&ResidueSpecifier {
            chain_id,
            residue_number,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct DesignConfig {
    pub forcefield: ForcefieldConfig,
    pub sampling: SamplingConfig,
    pub optimization: OptimizationConfig,
    pub design_spec: DesignSpec,
    pub neighbors_to_repack: ResidueSelection,
}

#[derive(Debug, Clone, PartialEq)]
pub enum AtomSelection {
    Residue(ResidueSpecifier),
    Chain(char),
    All,
}

#[derive(Debug, Clone, PartialEq)]
pub enum AnalysisType {
    Interaction {
        group1: AtomSelection,
        group2: AtomSelection,
    },
    ClashDetection {
        threshold_kcal_mol: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct AnalyzeConfig {
    pub forcefield: ForcefieldConfig,
    pub analysis_type: AnalysisType,
}
