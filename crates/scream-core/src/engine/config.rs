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

impl Default for PlacementSelection {
    fn default() -> Self {
        PlacementSelection::All
    }
}

impl Default for PlacementMode {
    fn default() -> Self {
        PlacementMode::Standard(PlacementSelection::default())
    }
}

impl Default for EngineConfig {
    fn default() -> Self {
        Self {
            mode: PlacementMode::default(),
            s_factor: 1.0,
            max_iterations: 10,
            convergence_tolerance: 0.01,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn placement_selection_defaults_to_all() {
        assert_eq!(PlacementSelection::default(), PlacementSelection::All);
    }

    #[test]
    fn placement_mode_defaults_to_standard_with_default_selection() {
        let default_mode = PlacementMode::default();
        let expected_mode = PlacementMode::Standard(PlacementSelection::default());

        match (default_mode, expected_mode) {
            (PlacementMode::Standard(s1), PlacementMode::Standard(s2)) => assert_eq!(s1, s2),
            _ => panic!("Default PlacementMode is not Standard"),
        }
    }

    #[test]
    fn engine_config_defaults_to_sensible_values() {
        let default_config = EngineConfig::default();

        assert_eq!(default_config.s_factor, 1.0);
        assert_eq!(default_config.max_iterations, 10);
        assert_eq!(default_config.convergence_tolerance, 0.01);

        match default_config.mode {
            PlacementMode::Standard(selection) => assert_eq!(selection, PlacementSelection::All),
            _ => panic!("Default EngineConfig mode is not Standard(All)"),
        }
    }
}
