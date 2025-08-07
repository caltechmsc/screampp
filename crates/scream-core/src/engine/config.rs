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
        ligand_residue: ResidueSpecifier,
        radius_angstroms: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct ConvergenceConfig {
    pub energy_threshold: f64,
    pub patience_iterations: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SimulatedAnnealingConfig {
    pub initial_temperature: f64,
    pub final_temperature: f64,
    pub cooling_rate: f64,
    pub steps_per_temperature: usize,
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
}

#[derive(Debug, Clone, PartialEq)]
pub struct OptimizationConfig {
    pub max_iterations: usize,
    pub num_solutions: usize,
    pub include_input_conformation: bool,
    pub convergence: ConvergenceConfig,
    pub simulated_annealing: Option<SimulatedAnnealingConfig>,
    pub final_refinement_iterations: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PlacementConfig {
    pub forcefield: ForcefieldConfig,
    pub sampling: SamplingConfig,
    pub optimization: OptimizationConfig,
    pub residues_to_optimize: ResidueSelection,
    pub topology_registry_path: PathBuf,
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
    pub topology_registry_path: PathBuf,
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
    pub topology_registry_path: PathBuf,
}

#[derive(Default)]
pub struct PlacementConfigBuilder {
    forcefield_path: Option<PathBuf>,
    delta_params_path: Option<PathBuf>,
    s_factor: Option<f64>,
    rotamer_library_path: Option<PathBuf>,
    topology_registry_path: Option<PathBuf>,
    max_iterations: Option<usize>,
    num_solutions: Option<usize>,
    include_input_conformation: Option<bool>,
    convergence_config: Option<ConvergenceConfig>,
    simulated_annealing_config: Option<SimulatedAnnealingConfig>,
    final_refinement_iterations: Option<usize>,
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
    pub fn topology_registry_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.topology_registry_path = Some(path.into());
        self
    }
    pub fn max_iterations(mut self, iterations: usize) -> Self {
        self.max_iterations = Some(iterations);
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
    pub fn convergence_config(mut self, config: ConvergenceConfig) -> Self {
        self.convergence_config = Some(config);
        self
    }
    pub fn simulated_annealing_config(mut self, config: Option<SimulatedAnnealingConfig>) -> Self {
        self.simulated_annealing_config = config;
        self
    }
    pub fn final_refinement_iterations(mut self, iterations: usize) -> Self {
        self.final_refinement_iterations = Some(iterations);
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
        };
        let optimization = OptimizationConfig {
            max_iterations: self
                .max_iterations
                .ok_or(ConfigError::MissingParameter("max_iterations"))?,
            num_solutions: self
                .num_solutions
                .ok_or(ConfigError::MissingParameter("num_solutions"))?,
            include_input_conformation: self.include_input_conformation.unwrap_or(false),
            convergence: self
                .convergence_config
                .ok_or(ConfigError::MissingParameter("convergence_config"))?,
            simulated_annealing: self.simulated_annealing_config,
            final_refinement_iterations: self
                .final_refinement_iterations
                .ok_or(ConfigError::MissingParameter("final_refinement_iterations"))?,
        };

        Ok(PlacementConfig {
            forcefield,
            sampling,
            optimization,
            residues_to_optimize: self
                .residues_to_optimize
                .ok_or(ConfigError::MissingParameter("residues_to_optimize"))?,
            topology_registry_path: self
                .topology_registry_path
                .ok_or(ConfigError::MissingParameter("topology_registry_path"))?,
        })
    }
}

#[derive(Default)]
pub struct DesignConfigBuilder {
    forcefield_path: Option<PathBuf>,
    delta_params_path: Option<PathBuf>,
    s_factor: Option<f64>,
    rotamer_library_path: Option<PathBuf>,
    topology_registry_path: Option<PathBuf>,
    max_iterations: Option<usize>,
    num_solutions: Option<usize>,
    include_input_conformation: Option<bool>,
    convergence_config: Option<ConvergenceConfig>,
    simulated_annealing_config: Option<SimulatedAnnealingConfig>,
    final_refinement_iterations: Option<usize>,
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
    pub fn topology_registry_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.topology_registry_path = Some(path.into());
        self
    }

    pub fn max_iterations(mut self, iterations: usize) -> Self {
        self.max_iterations = Some(iterations);
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
    pub fn convergence_config(mut self, config: ConvergenceConfig) -> Self {
        self.convergence_config = Some(config);
        self
    }
    pub fn simulated_annealing_config(mut self, config: SimulatedAnnealingConfig) -> Self {
        self.simulated_annealing_config = Some(config);
        self
    }
    pub fn final_refinement_iterations(mut self, iterations: usize) -> Self {
        self.final_refinement_iterations = Some(iterations);
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
        };
        let optimization = OptimizationConfig {
            max_iterations: self
                .max_iterations
                .ok_or(ConfigError::MissingParameter("max_iterations"))?,
            num_solutions: self
                .num_solutions
                .ok_or(ConfigError::MissingParameter("num_solutions"))?,
            include_input_conformation: self
                .include_input_conformation
                .ok_or(ConfigError::MissingParameter("include_input_conformation"))?,
            convergence: self
                .convergence_config
                .ok_or(ConfigError::MissingParameter("convergence_config"))?,
            simulated_annealing: self.simulated_annealing_config,
            final_refinement_iterations: self
                .final_refinement_iterations
                .ok_or(ConfigError::MissingParameter("final_refinement_iterations"))?,
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
            topology_registry_path: self
                .topology_registry_path
                .ok_or(ConfigError::MissingParameter("topology_registry_path"))?,
        })
    }
}

#[derive(Default)]
pub struct AnalyzeConfigBuilder {
    forcefield_path: Option<PathBuf>,
    delta_params_path: Option<PathBuf>,
    s_factor: Option<f64>,
    topology_registry_path: Option<PathBuf>,
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
    pub fn topology_registry_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.topology_registry_path = Some(path.into());
        self
    }
    pub fn analysis_type(mut self, analysis_type: AnalysisType) -> Self {
        self.analysis_type = Some(analysis_type);
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
            topology_registry_path: self
                .topology_registry_path
                .ok_or(ConfigError::MissingParameter("topology_registry_path"))?,
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
        let convergence = ConvergenceConfig {
            energy_threshold: 0.01,
            patience_iterations: 10,
        };
        let simulated_annealing = SimulatedAnnealingConfig {
            initial_temperature: 100.0,
            final_temperature: 1.0,
            cooling_rate: 0.95,
            steps_per_temperature: 50,
        };
        let builder = PlacementConfigBuilder::new()
            .forcefield_path("ff.dat")
            .delta_params_path("delta.dat")
            .s_factor(0.5)
            .rotamer_library_path("rot.lib")
            .placement_registry_path("reg.json")
            .max_iterations(100)
            .num_solutions(10)
            .include_input_conformation(true)
            .convergence_config(convergence.clone())
            .simulated_annealing_config(Some(simulated_annealing.clone()))
            .final_refinement_iterations(5)
            .residues_to_optimize(ResidueSelection::All);

        let config = builder.build().unwrap();
        assert_eq!(config.optimization.convergence, convergence);
        assert_eq!(
            config.optimization.simulated_annealing,
            Some(simulated_annealing)
        );
        assert_eq!(config.optimization.final_refinement_iterations, 5);
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
            .num_solutions(10)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.01,
                patience_iterations: 10,
            })
            .final_refinement_iterations(5);

        assert_eq!(
            builder_missing_residues.build().unwrap_err(),
            ConfigError::MissingParameter("residues_to_optimize")
        );
    }

    #[test]
    fn placement_config_builder_uses_default_values_for_optional_fields() {
        let convergence = ConvergenceConfig {
            energy_threshold: 0.01,
            patience_iterations: 10,
        };
        let builder = PlacementConfigBuilder::new()
            .forcefield_path("ff.dat")
            .delta_params_path("delta.dat")
            .s_factor(0.5)
            .rotamer_library_path("rot.lib")
            .placement_registry_path("reg.json")
            .max_iterations(100)
            .num_solutions(10)
            .convergence_config(convergence.clone())
            .final_refinement_iterations(0)
            .residues_to_optimize(ResidueSelection::All);

        let config = builder.build().unwrap();
        assert_eq!(config.optimization.simulated_annealing, None);
        assert_eq!(config.optimization.final_refinement_iterations, 0);
    }

    #[test]
    fn design_config_builder_builds_successfully_with_all_parameters() {
        let convergence = ConvergenceConfig {
            energy_threshold: 0.01,
            patience_iterations: 10,
        };
        let simulated_annealing = SimulatedAnnealingConfig {
            initial_temperature: 100.0,
            final_temperature: 1.0,
            cooling_rate: 0.95,
            steps_per_temperature: 50,
        };
        let design_spec = DesignSpec::new();
        let builder = DesignConfigBuilder::new()
            .forcefield_path("ff.dat")
            .delta_params_path("delta.dat")
            .s_factor(0.5)
            .rotamer_library_path("rot.lib")
            .placement_registry_path("reg.json")
            .max_iterations(100)
            .num_solutions(10)
            .include_input_conformation(true)
            .convergence_config(convergence.clone())
            .simulated_annealing_config(simulated_annealing.clone())
            .final_refinement_iterations(5)
            .design_spec(design_spec)
            .neighbors_to_repack(ResidueSelection::All);

        let config = builder.build().unwrap();
        assert_eq!(config.optimization.convergence, convergence);
        assert_eq!(
            config.optimization.simulated_annealing,
            Some(simulated_annealing)
        );
        assert_eq!(config.optimization.final_refinement_iterations, 5);
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
            .num_solutions(10)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.01,
                patience_iterations: 10,
            })
            .final_refinement_iterations(5)
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
