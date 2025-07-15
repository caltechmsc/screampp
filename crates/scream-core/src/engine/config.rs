use crate::core::models::residue::ResidueType;
use std::collections::HashMap;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error, PartialEq, Eq)]
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
        ligand_residue_name: String,
        radius_angstroms: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub enum AtomSelection {
    Residue(ResidueSpecifier),
    Chain(char),
    Ligand(String),
    All,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GlobalConfig {
    pub forcefield_path: PathBuf,
    pub rotamer_library_path: PathBuf,
    pub delta_params_path: PathBuf,
    pub s_factor: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PlacementConfig {
    pub global: GlobalConfig,
    pub input_structure_path: PathBuf,
    pub output_path: PathBuf,
    pub residues_to_optimize: ResidueSelection,
    pub max_iterations: usize,
    pub convergence_threshold: f64,
    pub num_solutions: usize,
}

#[derive(Default)]
pub struct PlacementConfigBuilder {
    forcefield_path: Option<PathBuf>,
    rotamer_library_path: Option<PathBuf>,
    delta_params_path: Option<PathBuf>,
    s_factor: Option<f64>,
    input_structure_path: Option<PathBuf>,
    output_path: Option<PathBuf>,
    residues_to_optimize: Option<ResidueSelection>,
    max_iterations: Option<usize>,
    convergence_threshold: Option<f64>,
    num_solutions: Option<usize>,
}

impl PlacementConfigBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn forcefield_path(mut self, path: PathBuf) -> Self {
        self.forcefield_path = Some(path);
        self
    }
    pub fn rotamer_library_path(mut self, path: PathBuf) -> Self {
        self.rotamer_library_path = Some(path);
        self
    }
    pub fn delta_params_path(mut self, path: PathBuf) -> Self {
        self.delta_params_path = Some(path);
        self
    }
    pub fn s_factor(mut self, factor: f64) -> Self {
        self.s_factor = Some(factor);
        self
    }
    pub fn input_structure_path(mut self, path: PathBuf) -> Self {
        self.input_structure_path = Some(path);
        self
    }
    pub fn output_path(mut self, path: PathBuf) -> Self {
        self.output_path = Some(path);
        self
    }
    pub fn residues_to_optimize(mut self, selection: ResidueSelection) -> Self {
        self.residues_to_optimize = Some(selection);
        self
    }
    pub fn max_iterations(mut self, iterations: usize) -> Self {
        self.max_iterations = Some(iterations);
        self
    }
    pub fn convergence_threshold(mut self, threshold: f64) -> Self {
        self.convergence_threshold = Some(threshold);
        self
    }
    pub fn num_solutions(mut self, n: usize) -> Self {
        self.num_solutions = Some(n);
        self
    }

    pub fn build(self) -> Result<PlacementConfig, ConfigError> {
        let global = GlobalConfig {
            forcefield_path: self
                .forcefield_path
                .ok_or(ConfigError::MissingParameter("forcefield_path"))?,
            rotamer_library_path: self
                .rotamer_library_path
                .ok_or(ConfigError::MissingParameter("rotamer_library_path"))?,
            delta_params_path: self
                .delta_params_path
                .ok_or(ConfigError::MissingParameter("delta_params_path"))?,
            s_factor: self
                .s_factor
                .ok_or(ConfigError::MissingParameter("s_factor"))?,
        };

        Ok(PlacementConfig {
            global,
            input_structure_path: self
                .input_structure_path
                .ok_or(ConfigError::MissingParameter("input_structure_path"))?,
            output_path: self
                .output_path
                .ok_or(ConfigError::MissingParameter("output_path"))?,
            residues_to_optimize: self
                .residues_to_optimize
                .ok_or(ConfigError::MissingParameter("residues_to_optimize"))?,
            max_iterations: self
                .max_iterations
                .ok_or(ConfigError::MissingParameter("max_iterations"))?,
            convergence_threshold: self
                .convergence_threshold
                .ok_or(ConfigError::MissingParameter("convergence_threshold"))?,
            num_solutions: self
                .num_solutions
                .ok_or(ConfigError::MissingParameter("num_solutions"))?,
        })
    }
}

pub type DesignSpec = HashMap<ResidueSpecifier, Vec<ResidueType>>;

#[derive(Debug, Clone, PartialEq)]
pub struct DesignConfig {
    pub global: GlobalConfig,
    pub template_structure_path: PathBuf,
    pub output_prefix: PathBuf,
    pub design_spec: DesignSpec,
    pub neighbors_to_repack: ResidueSelection,
    pub num_solutions: usize,
}

// TODO: Implement a builder for DesignConfig similar to PlacementConfigBuilder.

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
    pub global: GlobalConfig,
    pub input_structure_path: PathBuf,
    pub analysis_type: AnalysisType,
}

// TODO: Implement a builder for AnalyzeConfig similar to PlacementConfigBuilder.
