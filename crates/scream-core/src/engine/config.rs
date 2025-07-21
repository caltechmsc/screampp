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

#[derive(Default)]
pub struct PlacementConfigBuilder {
    forcefield_path: Option<PathBuf>,
    delta_params_path: Option<PathBuf>,
    s_factor: Option<f64>,
    rotamer_library_path: Option<PathBuf>,
    placement_registry_path: Option<PathBuf>,
    max_iterations: Option<usize>,
    convergence_threshold: Option<f64>,
    num_solutions: Option<usize>,
    include_input_conformation: Option<bool>,
    use_simulated_annealing: Option<bool>,
    residues_to_optimize: Option<ResidueSelection>,
}

impl PlacementConfigBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn forcefield_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.forcefield_path = Some(path.into());
        self
    }
    pub fn delta_params_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.delta_params_path = Some(path.into());
        self
    }
    pub fn s_factor(mut self, factor: f64) -> Self {
        self.s_factor = Some(factor);
        self
    }
    pub fn rotamer_library_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.rotamer_library_path = Some(path.into());
        self
    }
    pub fn placement_registry_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.placement_registry_path = Some(path.into());
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
    pub fn include_input_conformation(mut self, include: bool) -> Self {
        self.include_input_conformation = Some(include);
        self
    }
    pub fn use_simulated_annealing(mut self, use_sa: bool) -> Self {
        self.use_simulated_annealing = Some(use_sa);
        self
    }
    pub fn residues_to_optimize(mut self, selection: ResidueSelection) -> Self {
        self.residues_to_optimize = Some(selection);
        self
    }

    pub fn build(self) -> Result<PlacementConfig, ConfigError> {
        let forcefield = ForcefieldConfig {
            forcefield_path: self
                .forcefield_path
                .ok_or(ConfigError::MissingParameter("forcefield_path"))?,
            delta_params_path: self
                .delta_params_path
                .ok_or(ConfigError::MissingParameter("delta_params_path"))?,
            s_factor: self
                .s_factor
                .ok_or(ConfigError::MissingParameter("s_factor"))?,
        };
        let sampling = SamplingConfig {
            rotamer_library_path: self
                .rotamer_library_path
                .ok_or(ConfigError::MissingParameter("rotamer_library_path"))?,
            placement_registry_path: self
                .placement_registry_path
                .ok_or(ConfigError::MissingParameter("placement_registry_path"))?,
        };
        let optimization = OptimizationConfig {
            max_iterations: self
                .max_iterations
                .ok_or(ConfigError::MissingParameter("max_iterations"))?,
            convergence_threshold: self
                .convergence_threshold
                .ok_or(ConfigError::MissingParameter("convergence_threshold"))?,
            num_solutions: self
                .num_solutions
                .ok_or(ConfigError::MissingParameter("num_solutions"))?,
            include_input_conformation: self.include_input_conformation.unwrap_or(false),
            use_simulated_annealing: self.use_simulated_annealing.unwrap_or(false),
        };

        Ok(PlacementConfig {
            forcefield,
            sampling,
            optimization,
            residues_to_optimize: self
                .residues_to_optimize
                .ok_or(ConfigError::MissingParameter("residues_to_optimize"))?,
        })
    }
}

#[derive(Default)]
pub struct DesignConfigBuilder {
    forcefield_path: Option<PathBuf>,
    delta_params_path: Option<PathBuf>,
    s_factor: Option<f64>,
    rotamer_library_path: Option<PathBuf>,
    placement_registry_path: Option<PathBuf>,
    max_iterations: Option<usize>,
    convergence_threshold: Option<f64>,
    num_solutions: Option<usize>,
    include_input_conformation: Option<bool>,
    use_simulated_annealing: Option<bool>,
    design_spec: Option<DesignSpec>,
    neighbors_to_repack: Option<ResidueSelection>,
}

impl DesignConfigBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn forcefield_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.forcefield_path = Some(path.into());
        self
    }
    pub fn delta_params_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.delta_params_path = Some(path.into());
        self
    }
    pub fn s_factor(mut self, factor: f64) -> Self {
        self.s_factor = Some(factor);
        self
    }

    pub fn rotamer_library_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.rotamer_library_path = Some(path.into());
        self
    }
    pub fn placement_registry_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.placement_registry_path = Some(path.into());
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
    pub fn include_input_conformation(mut self, include: bool) -> Self {
        self.include_input_conformation = Some(include);
        self
    }
    pub fn use_simulated_annealing(mut self, use_sa: bool) -> Self {
        self.use_simulated_annealing = Some(use_sa);
        self
    }

    pub fn design_spec(mut self, spec: DesignSpec) -> Self {
        self.design_spec = Some(spec);
        self
    }
    pub fn neighbors_to_repack(mut self, selection: ResidueSelection) -> Self {
        self.neighbors_to_repack = Some(selection);
        self
    }

    pub fn build(self) -> Result<DesignConfig, ConfigError> {
        let forcefield = ForcefieldConfig {
            forcefield_path: self
                .forcefield_path
                .ok_or(ConfigError::MissingParameter("forcefield_path"))?,
            delta_params_path: self
                .delta_params_path
                .ok_or(ConfigError::MissingParameter("delta_params_path"))?,
            s_factor: self
                .s_factor
                .ok_or(ConfigError::MissingParameter("s_factor"))?,
        };
        let sampling = SamplingConfig {
            rotamer_library_path: self
                .rotamer_library_path
                .ok_or(ConfigError::MissingParameter("rotamer_library_path"))?,
            placement_registry_path: self
                .placement_registry_path
                .ok_or(ConfigError::MissingParameter("placement_registry_path"))?,
        };
        let optimization = OptimizationConfig {
            max_iterations: self
                .max_iterations
                .ok_or(ConfigError::MissingParameter("max_iterations"))?,
            convergence_threshold: self
                .convergence_threshold
                .ok_or(ConfigError::MissingParameter("convergence_threshold"))?,
            num_solutions: self
                .num_solutions
                .ok_or(ConfigError::MissingParameter("num_solutions"))?,
            include_input_conformation: self.include_input_conformation.unwrap_or(false),
            use_simulated_annealing: self.use_simulated_annealing.unwrap_or(false),
        };

        Ok(DesignConfig {
            forcefield,
            sampling,
            optimization,
            design_spec: self
                .design_spec
                .ok_or(ConfigError::MissingParameter("design_spec"))?,
            neighbors_to_repack: self
                .neighbors_to_repack
                .ok_or(ConfigError::MissingParameter("neighbors_to_repack"))?,
        })
    }
}

#[derive(Default)]
pub struct AnalyzeConfigBuilder {
    forcefield_path: Option<PathBuf>,
    delta_params_path: Option<PathBuf>,
    s_factor: Option<f64>,
    analysis_type: Option<AnalysisType>,
}

impl AnalyzeConfigBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn forcefield_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.forcefield_path = Some(path.into());
        self
    }
    pub fn delta_params_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.delta_params_path = Some(path.into());
        self
    }
    pub fn s_factor(mut self, factor: f64) -> Self {
        self.s_factor = Some(factor);
        self
    }
    pub fn analysis_type(mut self, analysis: AnalysisType) -> Self {
        self.analysis_type = Some(analysis);
        self
    }

    pub fn build(self) -> Result<AnalyzeConfig, ConfigError> {
        let forcefield = ForcefieldConfig {
            forcefield_path: self
                .forcefield_path
                .ok_or(ConfigError::MissingParameter("forcefield_path"))?,
            delta_params_path: self
                .delta_params_path
                .ok_or(ConfigError::MissingParameter("delta_params_path"))?,
            s_factor: self
                .s_factor
                .ok_or(ConfigError::MissingParameter("s_factor"))?,
        };

        Ok(AnalyzeConfig {
            forcefield,
            analysis_type: self
                .analysis_type
                .ok_or(ConfigError::MissingParameter("analysis_type"))?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn design_spec_ext_get_by_specifier_finds_existing_entry() {
        let specifier = ResidueSpecifier {
            chain_id: 'A',
            residue_number: 10,
        };
        let residue_types = vec![ResidueType::Alanine, ResidueType::Glycine];
        let mut design_spec: DesignSpec = HashMap::new();
        design_spec.insert(specifier.clone(), residue_types.clone());

        let result = design_spec.get_by_specifier('A', 10);
        assert_eq!(result, Some(&residue_types));
    }

    #[test]
    fn design_spec_ext_get_by_specifier_returns_none_for_missing_entry() {
        let design_spec: DesignSpec = HashMap::new();
        let result = design_spec.get_by_specifier('B', 20);
        assert_eq!(result, None);
    }

    #[test]
    fn placement_config_builder_builds_successfully_with_all_parameters() {
        let builder = PlacementConfigBuilder::new()
            .forcefield_path("ff.dat")
            .delta_params_path("delta.dat")
            .s_factor(0.5)
            .rotamer_library_path("rot.lib")
            .placement_registry_path("reg.json")
            .max_iterations(100)
            .convergence_threshold(0.01)
            .num_solutions(10)
            .include_input_conformation(true)
            .use_simulated_annealing(true)
            .residues_to_optimize(ResidueSelection::All);

        assert!(builder.build().is_ok());
    }

    #[test]
    fn placement_config_builder_fails_on_missing_required_parameter() {
        assert_eq!(
            PlacementConfigBuilder::new().build().unwrap_err(),
            ConfigError::MissingParameter("forcefield_path")
        );

        let builder_missing_residues = PlacementConfigBuilder::new()
            .forcefield_path("ff.dat")
            .delta_params_path("delta.dat")
            .s_factor(0.5)
            .rotamer_library_path("rot.lib")
            .placement_registry_path("reg.json")
            .max_iterations(100)
            .convergence_threshold(0.01)
            .num_solutions(10);

        assert_eq!(
            builder_missing_residues.build().unwrap_err(),
            ConfigError::MissingParameter("residues_to_optimize")
        );
    }

    #[test]
    fn placement_config_builder_uses_default_values_for_optional_booleans() {
        let config = PlacementConfigBuilder::new()
            .forcefield_path("ff.dat")
            .delta_params_path("delta.dat")
            .s_factor(0.5)
            .rotamer_library_path("rot.lib")
            .placement_registry_path("reg.json")
            .max_iterations(100)
            .convergence_threshold(0.01)
            .num_solutions(10)
            .residues_to_optimize(ResidueSelection::All)
            .build()
            .unwrap();

        assert!(!config.optimization.include_input_conformation);
        assert!(!config.optimization.use_simulated_annealing);
    }

    #[test]
    fn design_config_builder_builds_successfully_with_all_parameters() {
        let design_spec = DesignSpec::new();
        let builder = DesignConfigBuilder::new()
            .forcefield_path("ff.dat")
            .delta_params_path("delta.dat")
            .s_factor(0.5)
            .rotamer_library_path("rot.lib")
            .placement_registry_path("reg.json")
            .max_iterations(100)
            .convergence_threshold(0.01)
            .num_solutions(10)
            .design_spec(design_spec)
            .neighbors_to_repack(ResidueSelection::All);

        assert!(builder.build().is_ok());
    }

    #[test]
    fn design_config_builder_fails_on_missing_design_spec() {
        let builder = DesignConfigBuilder::new()
            .forcefield_path("ff.dat")
            .delta_params_path("delta.dat")
            .s_factor(0.5)
            .rotamer_library_path("rot.lib")
            .placement_registry_path("reg.json")
            .max_iterations(100)
            .convergence_threshold(0.01)
            .num_solutions(10)
            .neighbors_to_repack(ResidueSelection::All);

        assert_eq!(
            builder.build().unwrap_err(),
            ConfigError::MissingParameter("design_spec")
        );
    }

    #[test]
    fn analyze_config_builder_builds_successfully() {
        let builder = AnalyzeConfigBuilder::new()
            .forcefield_path("ff.dat")
            .delta_params_path("delta.dat")
            .s_factor(0.5)
            .analysis_type(AnalysisType::ClashDetection {
                threshold_kcal_mol: 1.0,
            });

        let config = builder.build().unwrap();

        assert_eq!(config.forcefield.forcefield_path, Path::new("ff.dat"));
        assert_eq!(
            config.analysis_type,
            AnalysisType::ClashDetection {
                threshold_kcal_mol: 1.0
            }
        );
    }

    #[test]
    fn analyze_config_builder_fails_on_missing_analysis_type() {
        let builder = AnalyzeConfigBuilder::new()
            .forcefield_path("ff.dat")
            .delta_params_path("delta.dat")
            .s_factor(0.5);

        assert_eq!(
            builder.build().unwrap_err(),
            ConfigError::MissingParameter("analysis_type")
        );
    }
}
