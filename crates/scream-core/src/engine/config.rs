use crate::core::forcefield::params::EnergyWeights;
use crate::core::models::residue::ResidueType;
use std::collections::HashMap;
use std::path::PathBuf;
use thiserror::Error;

/// Errors that can occur during configuration building and validation.
///
/// This enum defines the possible errors that may arise when constructing
/// configuration objects for placement, design, or analysis operations in SCREAM++.
/// Each variant provides specific information about what went wrong to aid in
/// debugging configuration issues.
#[derive(Debug, Error, PartialEq, Eq, Clone)]
pub enum ConfigError {
    /// A required configuration parameter was not provided.
    ///
    /// This error is returned when attempting to build a configuration without
    /// supplying all mandatory parameters. The parameter name is included in
    /// the error message to help identify what needs to be specified.
    #[error("Missing required parameter: {0}")]
    MissingParameter(&'static str),
}

/// Uniquely identifies a residue within a protein structure.
///
/// This struct serves as a key to reference specific amino acid residues in
/// molecular structures. It combines the chain identifier with the residue
/// sequence number to provide unambiguous targeting for various computational
/// operations such as side-chain placement, design, or analysis.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ResidueSpecifier {
    /// The character identifier of the protein chain containing this residue.
    pub chain_id: char,
    /// The sequential position number of the residue within its chain.
    pub residue_number: isize,
}

/// Defines how to select residues for computational operations.
///
/// This enum provides flexible ways to specify which residues should be
/// included or excluded from operations like side-chain optimization or
/// protein design. It supports selecting all residues, specific lists,
/// or regions around a ligand binding site.
#[derive(Debug, Clone, PartialEq)]
pub enum ResidueSelection {
    /// Select all residues in the structure.
    All,
    /// Select specific residues with optional exclusions.
    ///
    /// # Arguments
    ///
    /// * `include` - List of residues to include in the selection.
    /// * `exclude` - List of residues to exclude from the selection.
    List {
        include: Vec<ResidueSpecifier>,
        exclude: Vec<ResidueSpecifier>,
    },
    /// Select residues within a specified radius of a ligand.
    ///
    /// # Arguments
    ///
    /// * `ligand_residue` - The residue identifier of the ligand.
    /// * `radius_angstroms` - The radius in Angstroms around the ligand.
    LigandBindingSite {
        ligand_residue: ResidueSpecifier,
        radius_angstroms: f64,
    },
}

/// Configuration parameters for optimization convergence criteria.
///
/// This struct defines the thresholds and patience settings used to determine
/// when an optimization algorithm has converged to a stable solution. It helps
/// balance computational efficiency with solution quality in iterative
/// refinement processes.
#[derive(Debug, Clone, PartialEq)]
pub struct ConvergenceConfig {
    /// The minimum energy difference threshold for convergence.
    ///
    /// Optimization stops when the energy change between iterations
    /// falls below this value (in kcal/mol).
    pub energy_threshold: f64,
    /// Number of iterations to wait before checking convergence.
    ///
    /// This prevents premature termination due to temporary energy fluctuations
    /// during the early stages of optimization.
    pub patience_iterations: usize,
}

/// Configuration for simulated annealing optimization.
///
/// This struct contains the parameters that control the simulated annealing
/// algorithm used for global optimization of side-chain conformations. The
/// algorithm uses a temperature schedule to explore the conformational space
/// and gradually converge to optimal solutions.
#[derive(Debug, Clone, PartialEq)]
pub struct SimulatedAnnealingConfig {
    /// The starting temperature for the annealing process.
    ///
    /// Higher values allow more exploration of the conformational space
    /// at the beginning of optimization.
    pub initial_temperature: f64,
    /// The final temperature at which annealing terminates.
    ///
    /// Lower values focus the search on local optima near the end.
    pub final_temperature: f64,
    /// The factor by which temperature decreases each step.
    ///
    /// Values between 0.8 and 0.99 are typical, with lower values
    /// providing faster cooling.
    pub cooling_rate: f64,
    /// Number of optimization steps performed at each temperature level.
    pub steps_per_temperature: usize,
}

/// Configuration for force field parameters and energy calculations.
///
/// This struct specifies the force field files and parameters used for
/// molecular mechanics energy calculations. It includes the main force field
/// parameters, delta parameters for side-chain corrections, and energy
/// weighting factors.
#[derive(Debug, Clone, PartialEq)]
pub struct ForcefieldConfig {
    /// Path to the main force field parameter file.
    pub forcefield_path: PathBuf,
    /// Path to the delta parameter file for side-chain corrections.
    pub delta_params_path: PathBuf,
    /// Scaling factor for non-bonded interactions.
    ///
    /// This parameter adjusts the strength of van der Waals and electrostatic
    /// interactions in the force field.
    pub s_factor: f64,
    /// Relative weights for different energy components.
    pub energy_weights: EnergyWeights,
}

/// Configuration for rotamer sampling during optimization.
///
/// This struct specifies the rotamer library used to sample side-chain
/// conformations during placement and design operations. Rotamer libraries
/// contain pre-computed, energetically favorable side-chain orientations.
#[derive(Debug, Clone, PartialEq)]
pub struct SamplingConfig {
    /// Path to the rotamer library file.
    pub rotamer_library_path: PathBuf,
}

/// Configuration for the overall optimization process.
///
/// This struct defines the parameters that control the optimization algorithm,
/// including iteration limits, solution generation, convergence criteria,
/// and optional simulated annealing settings.
#[derive(Debug, Clone, PartialEq)]
pub struct OptimizationConfig {
    /// Maximum number of optimization iterations.
    pub max_iterations: usize,
    /// Number of solution candidates to generate.
    pub num_solutions: usize,
    /// Whether to include the input conformation as a candidate.
    pub include_input_conformation: bool,
    /// Convergence criteria for optimization termination.
    pub convergence: ConvergenceConfig,
    /// Optional simulated annealing configuration.
    ///
    /// If None, standard optimization without annealing is used.
    pub simulated_annealing: Option<SimulatedAnnealingConfig>,
    /// Number of final refinement iterations after convergence.
    pub final_refinement_iterations: usize,
}

/// Complete configuration for side-chain placement operations.
///
/// This struct encapsulates all parameters needed to perform automated
/// side-chain placement on a protein structure. It includes force field
/// settings, sampling parameters, optimization controls, and residue
/// selection criteria.
#[derive(Debug, Clone, PartialEq)]
pub struct PlacementConfig {
    /// Force field configuration for energy calculations.
    pub forcefield: ForcefieldConfig,
    /// Rotamer sampling configuration.
    pub sampling: SamplingConfig,
    /// Optimization algorithm parameters.
    pub optimization: OptimizationConfig,
    /// Selection of residues to optimize.
    pub residues_to_optimize: ResidueSelection,
    /// Path to the topology registry file.
    pub topology_registry_path: PathBuf,
}

/// Type alias for specifying amino acid mutations in protein design.
///
/// This type maps residue specifiers to lists of allowed amino acid types,
/// defining the design space for computational protein design operations.
/// Each entry specifies which residue positions can mutate to which amino acids.
pub type DesignSpec = HashMap<ResidueSpecifier, Vec<ResidueType>>;

/// Extension trait for DesignSpec providing convenient access methods.
///
/// This trait adds utility methods to the DesignSpec type alias to simplify
/// common operations when working with design specifications in protein
/// design workflows.
pub trait DesignSpecExt {
    /// Retrieves the allowed amino acid types for a specific residue.
    ///
    /// # Arguments
    ///
    /// * `chain_id` - The chain identifier of the residue.
    /// * `residue_number` - The sequence number of the residue.
    ///
    /// # Return
    ///
    /// Returns `Some(Vec<ResidueType>)` if the residue is in the design spec,
    /// or `None` if it is not specified.
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

/// Complete configuration for computational protein design operations.
///
/// This struct defines all parameters required for protein design, including
/// force field settings, rotamer sampling, optimization parameters, the design
/// specification, and neighbor repacking criteria.
#[derive(Debug, Clone, PartialEq)]
pub struct DesignConfig {
    /// Force field configuration for energy calculations.
    pub forcefield: ForcefieldConfig,
    /// Rotamer sampling configuration.
    pub sampling: SamplingConfig,
    /// Optimization algorithm parameters.
    pub optimization: OptimizationConfig,
    /// Specification of allowed mutations at each position.
    pub design_spec: DesignSpec,
    /// Selection of neighboring residues to repack during design.
    pub neighbors_to_repack: ResidueSelection,
    /// Path to the topology registry file.
    pub topology_registry_path: PathBuf,
}

/// Defines how to select atoms for analysis operations.
///
/// This enum provides ways to specify which atoms should be included in
/// molecular analysis calculations, such as interaction energy computations
/// or clash detection.
#[derive(Debug, Clone, PartialEq)]
pub enum AtomSelection {
    /// Select a single residue's atoms.
    Residue(ResidueSpecifier),
    /// Select all atoms in a specific chain.
    Chain(char),
    /// Select all atoms in the structure.
    All,
}

/// Specifies the type of molecular analysis to perform.
///
/// This enum defines the available analysis operations that can be performed
/// on molecular structures, including interaction energy calculations and
/// steric clash detection.
#[derive(Debug, Clone, PartialEq)]
pub enum AnalysisType {
    /// Calculate interaction energies between two atom groups.
    ///
    /// # Arguments
    ///
    /// * `group1` - First group of atoms for interaction calculation.
    /// * `group2` - Second group of atoms for interaction calculation.
    Interaction {
        group1: AtomSelection,
        group2: AtomSelection,
    },
    /// Detect steric clashes above a threshold.
    ///
    /// # Arguments
    ///
    /// * `threshold_kcal_mol` - Energy threshold for clash detection in kcal/mol.
    ClashDetection { threshold_kcal_mol: f64 },
}

/// Complete configuration for molecular analysis operations.
///
/// This struct contains all parameters needed to perform analysis on
/// molecular structures, including force field settings and the specific
/// type of analysis to conduct.
#[derive(Debug, Clone, PartialEq)]
pub struct AnalyzeConfig {
    /// Force field configuration for energy calculations.
    pub forcefield: ForcefieldConfig,
    /// The type of analysis to perform.
    pub analysis_type: AnalysisType,
    /// Path to the topology registry file.
    pub topology_registry_path: PathBuf,
}

/// Builder pattern implementation for constructing PlacementConfig.
///
/// This struct provides a fluent interface for building placement configurations
/// with validation. It ensures all required parameters are provided and
/// uses sensible defaults for optional fields.
#[derive(Default)]
pub struct PlacementConfigBuilder {
    forcefield_path: Option<PathBuf>,
    delta_params_path: Option<PathBuf>,
    s_factor: Option<f64>,
    energy_weights: Option<EnergyWeights>,
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
    /// Creates a new builder with default values.
    ///
    /// # Return
    ///
    /// Returns a new `PlacementConfigBuilder` instance with all fields unset.
    pub fn new() -> Self {
        Self::default()
    }
    /// Sets the force field parameter file path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the force field parameter file.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn forcefield_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.forcefield_path = Some(path.into());
        self
    }
    /// Sets the delta parameters file path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the delta parameters file.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn delta_params_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.delta_params_path = Some(path.into());
        self
    }
    /// Sets the scaling factor for non-bonded interactions.
    ///
    /// # Arguments
    ///
    /// * `factor` - The scaling factor value.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn s_factor(mut self, factor: f64) -> Self {
        self.s_factor = Some(factor);
        self
    }
    /// Sets the energy weights for different force field components.
    ///
    /// # Arguments
    ///
    /// * `weights` - The energy weights configuration.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn energy_weights(mut self, weights: EnergyWeights) -> Self {
        self.energy_weights = Some(weights);
        self
    }
    /// Sets the rotamer library file path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the rotamer library file.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn rotamer_library_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.rotamer_library_path = Some(path.into());
        self
    }
    /// Sets the topology registry file path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the topology registry file.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn topology_registry_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.topology_registry_path = Some(path.into());
        self
    }
    /// Sets the maximum number of optimization iterations.
    ///
    /// # Arguments
    ///
    /// * `iterations` - The maximum number of iterations.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn max_iterations(mut self, iterations: usize) -> Self {
        self.max_iterations = Some(iterations);
        self
    }
    /// Sets the number of solution candidates to generate.
    ///
    /// # Arguments
    ///
    /// * `n` - The number of solutions.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn num_solutions(mut self, n: usize) -> Self {
        self.num_solutions = Some(n);
        self
    }
    /// Sets whether to include the input conformation as a candidate.
    ///
    /// # Arguments
    ///
    /// * `include` - Whether to include the input conformation.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn include_input_conformation(mut self, include: bool) -> Self {
        self.include_input_conformation = Some(include);
        self
    }
    /// Sets the convergence configuration.
    ///
    /// # Arguments
    ///
    /// * `config` - The convergence configuration.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn convergence_config(mut self, config: ConvergenceConfig) -> Self {
        self.convergence_config = Some(config);
        self
    }
    /// Sets the simulated annealing configuration.
    ///
    /// # Arguments
    ///
    /// * `config` - The simulated annealing configuration, or None to disable.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn simulated_annealing_config(mut self, config: Option<SimulatedAnnealingConfig>) -> Self {
        self.simulated_annealing_config = config;
        self
    }
    /// Sets the number of final refinement iterations.
    ///
    /// # Arguments
    ///
    /// * `iterations` - The number of refinement iterations.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn final_refinement_iterations(mut self, iterations: usize) -> Self {
        self.final_refinement_iterations = Some(iterations);
        self
    }
    /// Sets the residue selection for optimization.
    ///
    /// # Arguments
    ///
    /// * `selection` - The residue selection criteria.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn residues_to_optimize(mut self, selection: ResidueSelection) -> Self {
        self.residues_to_optimize = Some(selection);
        self
    }

    /// Builds the PlacementConfig from the current builder state.
    ///
    /// # Return
    ///
    /// Returns `Ok(PlacementConfig)` if all required parameters are set,
    /// or `Err(ConfigError)` if any required parameter is missing.
    ///
    /// # Errors
    ///
    /// Returns `ConfigError::MissingParameter` if any required field is not set.
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
            energy_weights: self.energy_weights.unwrap_or_default(),
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

/// Builder pattern implementation for constructing DesignConfig.
///
/// This struct provides a fluent interface for building design configurations
/// with validation. It ensures all required parameters are provided and
/// uses sensible defaults for optional fields.
#[derive(Default)]
pub struct DesignConfigBuilder {
    forcefield_path: Option<PathBuf>,
    delta_params_path: Option<PathBuf>,
    s_factor: Option<f64>,
    energy_weights: Option<EnergyWeights>,
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
    /// Creates a new builder with default values.
    ///
    /// # Return
    ///
    /// Returns a new `DesignConfigBuilder` instance with all fields unset.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the force field parameter file path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the force field parameter file.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn forcefield_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.forcefield_path = Some(path.into());
        self
    }
    /// Sets the delta parameters file path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the delta parameters file.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn delta_params_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.delta_params_path = Some(path.into());
        self
    }
    /// Sets the scaling factor for non-bonded interactions.
    ///
    /// # Arguments
    ///
    /// * `factor` - The scaling factor value.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn s_factor(mut self, factor: f64) -> Self {
        self.s_factor = Some(factor);
        self
    }
    /// Sets the energy weights for different force field components.
    ///
    /// # Arguments
    ///
    /// * `weights` - The energy weights configuration.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn energy_weights(mut self, weights: EnergyWeights) -> Self {
        self.energy_weights = Some(weights);
        self
    }

    /// Sets the rotamer library file path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the rotamer library file.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn rotamer_library_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.rotamer_library_path = Some(path.into());
        self
    }
    /// Sets the topology registry file path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the topology registry file.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn topology_registry_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.topology_registry_path = Some(path.into());
        self
    }

    /// Sets the maximum number of optimization iterations.
    ///
    /// # Arguments
    ///
    /// * `iterations` - The maximum number of iterations.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn max_iterations(mut self, iterations: usize) -> Self {
        self.max_iterations = Some(iterations);
        self
    }
    /// Sets the number of solution candidates to generate.
    ///
    /// # Arguments
    ///
    /// * `n` - The number of solutions.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn num_solutions(mut self, n: usize) -> Self {
        self.num_solutions = Some(n);
        self
    }
    /// Sets whether to include the input conformation as a candidate.
    ///
    /// # Arguments
    ///
    /// * `include` - Whether to include the input conformation.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn include_input_conformation(mut self, include: bool) -> Self {
        self.include_input_conformation = Some(include);
        self
    }
    /// Sets the convergence configuration.
    ///
    /// # Arguments
    ///
    /// * `config` - The convergence configuration.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn convergence_config(mut self, config: ConvergenceConfig) -> Self {
        self.convergence_config = Some(config);
        self
    }
    /// Sets the simulated annealing configuration.
    ///
    /// # Arguments
    ///
    /// * `config` - The simulated annealing configuration.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn simulated_annealing_config(mut self, config: SimulatedAnnealingConfig) -> Self {
        self.simulated_annealing_config = Some(config);
        self
    }
    /// Sets the number of final refinement iterations.
    ///
    /// # Arguments
    ///
    /// * `iterations` - The number of refinement iterations.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn final_refinement_iterations(mut self, iterations: usize) -> Self {
        self.final_refinement_iterations = Some(iterations);
        self
    }

    /// Sets the design specification.
    ///
    /// # Arguments
    ///
    /// * `spec` - The design specification mapping residues to allowed types.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn design_spec(mut self, spec: DesignSpec) -> Self {
        self.design_spec = Some(spec);
        self
    }
    /// Sets the residue selection for neighbors to repack.
    ///
    /// # Arguments
    ///
    /// * `selection` - The residue selection criteria for neighbors.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn neighbors_to_repack(mut self, selection: ResidueSelection) -> Self {
        self.neighbors_to_repack = Some(selection);
        self
    }

    /// Builds the DesignConfig from the current builder state.
    ///
    /// # Return
    ///
    /// Returns `Ok(DesignConfig)` if all required parameters are set,
    /// or `Err(ConfigError)` if any required parameter is missing.
    ///
    /// # Errors
    ///
    /// Returns `ConfigError::MissingParameter` if any required field is not set.
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
            energy_weights: self.energy_weights.unwrap_or_default(),
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

/// Builder pattern implementation for constructing AnalyzeConfig.
///
/// This struct provides a fluent interface for building analysis configurations
/// with validation. It ensures all required parameters are provided and
/// uses sensible defaults for optional fields.
#[derive(Default)]
pub struct AnalyzeConfigBuilder {
    forcefield_path: Option<PathBuf>,
    delta_params_path: Option<PathBuf>,
    s_factor: Option<f64>,
    energy_weights: Option<EnergyWeights>,
    topology_registry_path: Option<PathBuf>,
    analysis_type: Option<AnalysisType>,
}

impl AnalyzeConfigBuilder {
    /// Creates a new builder with default values.
    ///
    /// # Return
    ///
    /// Returns a new `AnalyzeConfigBuilder` instance with all fields unset.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the force field parameter file path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the force field parameter file.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn forcefield_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.forcefield_path = Some(path.into());
        self
    }
    /// Sets the delta parameters file path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the delta parameters file.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn delta_params_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.delta_params_path = Some(path.into());
        self
    }
    /// Sets the scaling factor for non-bonded interactions.
    ///
    /// # Arguments
    ///
    /// * `factor` - The scaling factor value.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn s_factor(mut self, factor: f64) -> Self {
        self.s_factor = Some(factor);
        self
    }
    /// Sets the energy weights for different force field components.
    ///
    /// # Arguments
    ///
    /// * `weights` - The energy weights configuration.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn energy_weights(mut self, weights: EnergyWeights) -> Self {
        self.energy_weights = Some(weights);
        self
    }
    /// Sets the topology registry file path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the topology registry file.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn topology_registry_path(mut self, path: impl Into<PathBuf>) -> Self {
        self.topology_registry_path = Some(path.into());
        self
    }
    /// Sets the type of analysis to perform.
    ///
    /// # Arguments
    ///
    /// * `analysis_type` - The analysis type configuration.
    ///
    /// # Return
    ///
    /// Returns the builder for method chaining.
    pub fn analysis_type(mut self, analysis_type: AnalysisType) -> Self {
        self.analysis_type = Some(analysis_type);
        self
    }

    /// Builds the AnalyzeConfig from the current builder state.
    ///
    /// # Return
    ///
    /// Returns `Ok(AnalyzeConfig)` if all required parameters are set,
    /// or `Err(ConfigError)` if any required parameter is missing.
    ///
    /// # Errors
    ///
    /// Returns `ConfigError::MissingParameter` if any required field is not set.
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
            energy_weights: self.energy_weights.unwrap_or_default(),
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
            .topology_registry_path("topology.dat")
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
            .topology_registry_path("topology.dat")
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
            .topology_registry_path("topology.dat")
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
            .topology_registry_path("topology.dat")
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
            .topology_registry_path("topology.dat")
            .max_iterations(100)
            .num_solutions(10)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.01,
                patience_iterations: 10,
            })
            .final_refinement_iterations(5)
            .neighbors_to_repack(ResidueSelection::All)
            .include_input_conformation(false);

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
            .topology_registry_path("topology.dat")
            .analysis_type(AnalysisType::ClashDetection {
                threshold_kcal_mol: 1.0,
            });

        let config = builder.build().unwrap();

        assert_eq!(config.forcefield.forcefield_path, Path::new("ff.dat"));
        assert_eq!(config.topology_registry_path, Path::new("topology.dat"));
        assert!(matches!(
            config.analysis_type,
            AnalysisType::ClashDetection {
                threshold_kcal_mol: 1.0
            }
        ));
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
