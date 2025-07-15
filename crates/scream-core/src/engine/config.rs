use crate::core::models::residue::ResidueType;
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
pub struct ScoringConfig {
    pub forcefield_path: PathBuf,
    pub rotamer_library_path: PathBuf,
    pub delta_params_path: PathBuf,
    pub s_factor: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct OptimizationConfig {
    pub max_iterations: usize,
    pub convergence_threshold: f64,
    pub num_solutions: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PlacementConfig {
    pub scoring: ScoringConfig,
    pub optimization: OptimizationConfig,
    pub residues_to_optimize: ResidueSelection,
}

#[derive(Default)]
pub struct PlacementConfigBuilder {
    forcefield_path: Option<PathBuf>,
    rotamer_library_path: Option<PathBuf>,
    delta_params_path: Option<PathBuf>,
    s_factor: Option<f64>,
    max_iterations: Option<usize>,
    convergence_threshold: Option<f64>,
    num_solutions: Option<usize>,
    residues_to_optimize: Option<ResidueSelection>,
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
    pub fn residues_to_optimize(mut self, selection: ResidueSelection) -> Self {
        self.residues_to_optimize = Some(selection);
        self
    }

    pub fn build(self) -> Result<PlacementConfig, ConfigError> {
        let scoring = ScoringConfig {
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
        };
        Ok(PlacementConfig {
            scoring,
            optimization,
            residues_to_optimize: self
                .residues_to_optimize
                .ok_or(ConfigError::MissingParameter("residues_to_optimize"))?,
        })
    }
}

pub type DesignSpec = HashMap<ResidueSpecifier, Vec<ResidueType>>;

#[derive(Debug, Clone, PartialEq)]
pub struct DesignConfig {
    pub scoring: ScoringConfig,
    pub optimization: OptimizationConfig,
    pub design_spec: DesignSpec,
    pub neighbors_to_repack: ResidueSelection,
}

#[derive(Default)]
pub struct DesignConfigBuilder {
    forcefield_path: Option<PathBuf>,
    rotamer_library_path: Option<PathBuf>,
    delta_params_path: Option<PathBuf>,
    s_factor: Option<f64>,
    max_iterations: Option<usize>,
    convergence_threshold: Option<f64>,
    num_solutions: Option<usize>,
    design_spec: Option<DesignSpec>,
    neighbors_to_repack: Option<ResidueSelection>,
}

impl DesignConfigBuilder {
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
    pub fn design_spec(mut self, spec: DesignSpec) -> Self {
        self.design_spec = Some(spec);
        self
    }
    pub fn neighbors_to_repack(mut self, selection: ResidueSelection) -> Self {
        self.neighbors_to_repack = Some(selection);
        self
    }

    pub fn build(self) -> Result<DesignConfig, ConfigError> {
        let scoring = ScoringConfig {
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
        };
        Ok(DesignConfig {
            scoring,
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
    pub scoring: ScoringConfig,
    pub analysis_type: AnalysisType,
}

#[derive(Default)]
pub struct AnalyzeConfigBuilder {
    forcefield_path: Option<PathBuf>,
    rotamer_library_path: Option<PathBuf>,
    delta_params_path: Option<PathBuf>,
    s_factor: Option<f64>,
    analysis_type: Option<AnalysisType>,
}

impl AnalyzeConfigBuilder {
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
    pub fn analysis_type(mut self, analysis: AnalysisType) -> Self {
        self.analysis_type = Some(analysis);
        self
    }

    pub fn build(self) -> Result<AnalyzeConfig, ConfigError> {
        let scoring = ScoringConfig {
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
        Ok(AnalyzeConfig {
            scoring,
            analysis_type: self
                .analysis_type
                .ok_or(ConfigError::MissingParameter("analysis_type"))?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::residue::ResidueType;
    use std::path::PathBuf;

    fn common_builder_setup(
        builder: PlacementConfigBuilder,
    ) -> (PlacementConfigBuilder, PlacementConfig) {
        let forcefield_path = PathBuf::from("ff.dat");
        let rotamer_library_path = PathBuf::from("rot.lib");
        let delta_params_path = PathBuf::from("delta.params");
        let residues_to_optimize = ResidueSelection::All;

        let configured_builder = builder
            .forcefield_path(forcefield_path.clone())
            .rotamer_library_path(rotamer_library_path.clone())
            .delta_params_path(delta_params_path.clone())
            .s_factor(0.5)
            .max_iterations(100)
            .convergence_threshold(0.01)
            .num_solutions(10)
            .residues_to_optimize(residues_to_optimize.clone());

        let expected_config = PlacementConfig {
            scoring: ScoringConfig {
                forcefield_path,
                rotamer_library_path,
                delta_params_path,
                s_factor: 0.5,
            },
            optimization: OptimizationConfig {
                max_iterations: 100,
                convergence_threshold: 0.01,
                num_solutions: 10,
            },
            residues_to_optimize,
        };

        (configured_builder, expected_config)
    }

    #[test]
    fn placement_config_builder_succeeds_with_all_parameters() {
        let (builder, expected_config) = common_builder_setup(PlacementConfigBuilder::new());
        let config = builder.build().unwrap();
        assert_eq!(config, expected_config);
    }

    #[test]
    fn placement_config_builder_fails_on_missing_forcefield_path() {
        let (builder, _) = common_builder_setup(PlacementConfigBuilder::new());
        let incomplete_builder = PlacementConfigBuilder {
            forcefield_path: None,
            ..builder
        };
        let result = incomplete_builder.build();
        assert_eq!(
            result,
            Err(ConfigError::MissingParameter("forcefield_path"))
        );
    }

    #[test]
    fn placement_config_builder_fails_on_missing_residues_to_optimize() {
        let (builder, _) = common_builder_setup(PlacementConfigBuilder::new());
        let incomplete_builder = PlacementConfigBuilder {
            residues_to_optimize: None,
            ..builder
        };
        let result = incomplete_builder.build();
        assert_eq!(
            result,
            Err(ConfigError::MissingParameter("residues_to_optimize"))
        );
    }

    #[test]
    fn design_config_builder_succeeds_with_all_parameters() {
        let forcefield_path = PathBuf::from("ff.dat");
        let rotamer_library_path = PathBuf::from("rot.lib");
        let delta_params_path = PathBuf::from("delta.params");
        let neighbors_to_repack = ResidueSelection::All;
        let mut design_spec = DesignSpec::new();
        design_spec.insert(
            ResidueSpecifier {
                chain_id: 'A',
                residue_number: 10,
            },
            vec![ResidueType::Alanine, ResidueType::Glycine],
        );

        let builder = DesignConfigBuilder::new()
            .forcefield_path(forcefield_path.clone())
            .rotamer_library_path(rotamer_library_path.clone())
            .delta_params_path(delta_params_path.clone())
            .s_factor(0.5)
            .max_iterations(100)
            .convergence_threshold(0.01)
            .num_solutions(10)
            .design_spec(design_spec.clone())
            .neighbors_to_repack(neighbors_to_repack.clone());

        let config = builder.build().unwrap();

        let expected_config = DesignConfig {
            scoring: ScoringConfig {
                forcefield_path,
                rotamer_library_path,
                delta_params_path,
                s_factor: 0.5,
            },
            optimization: OptimizationConfig {
                max_iterations: 100,
                convergence_threshold: 0.01,
                num_solutions: 10,
            },
            design_spec,
            neighbors_to_repack,
        };

        assert_eq!(config, expected_config);
    }

    #[test]
    fn design_config_builder_fails_on_missing_design_spec() {
        let builder = DesignConfigBuilder::new()
            .forcefield_path(PathBuf::from("ff.dat"))
            .rotamer_library_path(PathBuf::from("rot.lib"))
            .delta_params_path(PathBuf::from("delta.params"))
            .s_factor(0.5)
            .max_iterations(100)
            .convergence_threshold(0.01)
            .num_solutions(10)
            .neighbors_to_repack(ResidueSelection::All);

        let result = builder.build();
        assert_eq!(result, Err(ConfigError::MissingParameter("design_spec")));
    }

    #[test]
    fn analyze_config_builder_succeeds_with_all_parameters() {
        let forcefield_path = PathBuf::from("ff.dat");
        let rotamer_library_path = PathBuf::from("rot.lib");
        let delta_params_path = PathBuf::from("delta.params");
        let analysis_type = AnalysisType::ClashDetection {
            threshold_kcal_mol: 2.5,
        };

        let builder = AnalyzeConfigBuilder::new()
            .forcefield_path(forcefield_path.clone())
            .rotamer_library_path(rotamer_library_path.clone())
            .delta_params_path(delta_params_path.clone())
            .s_factor(1.0)
            .analysis_type(analysis_type.clone());

        let config = builder.build().unwrap();

        let expected_config = AnalyzeConfig {
            scoring: ScoringConfig {
                forcefield_path,
                rotamer_library_path,
                delta_params_path,
                s_factor: 1.0,
            },
            analysis_type,
        };

        assert_eq!(config, expected_config);
    }

    #[test]
    fn analyze_config_builder_fails_on_missing_analysis_type() {
        let builder = AnalyzeConfigBuilder::new()
            .forcefield_path(PathBuf::from("ff.dat"))
            .rotamer_library_path(PathBuf::from("rot.lib"))
            .delta_params_path(PathBuf::from("delta.params"))
            .s_factor(1.0);

        let result = builder.build();
        assert_eq!(result, Err(ConfigError::MissingParameter("analysis_type")));
    }

    #[test]
    fn analyze_config_builder_fails_on_missing_s_factor() {
        let builder = AnalyzeConfigBuilder::new()
            .forcefield_path(PathBuf::from("ff.dat"))
            .rotamer_library_path(PathBuf::from("rot.lib"))
            .delta_params_path(PathBuf::from("delta.params"))
            .analysis_type(AnalysisType::ClashDetection {
                threshold_kcal_mol: 2.5,
            });

        let result = builder.build();
        assert_eq!(result, Err(ConfigError::MissingParameter("s_factor")));
    }
}
