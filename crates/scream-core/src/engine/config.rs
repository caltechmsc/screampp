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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::residue::ResidueType;
    use std::collections::HashMap;
    use std::path::PathBuf;

    #[test]
    fn input_files_can_be_constructed() {
        let files = InputFiles {
            structure: PathBuf::from("input.pdb"),
            non_bonded_params: PathBuf::from("non_bonded.toml"),
            topology_params: PathBuf::from("topology.toml"),
            placement_params: PathBuf::from("placement.toml"),
            charge_params: PathBuf::from("charges.toml"),
            delta_params: PathBuf::from("deltas.toml"),
            rotamer_library: PathBuf::from("rotamers.toml"),
        };
        assert_eq!(files.structure.to_str(), Some("input.pdb"));
        assert_eq!(files.non_bonded_params.to_str(), Some("non_bonded.toml"));
        assert_eq!(files.topology_params.to_str(), Some("topology.toml"));
        assert_eq!(files.placement_params.to_str(), Some("placement.toml"));
        assert_eq!(files.charge_params.to_str(), Some("charges.toml"));
        assert_eq!(files.delta_params.to_str(), Some("deltas.toml"));
        assert_eq!(files.rotamer_library.to_str(), Some("rotamers.toml"));
    }

    #[test]
    fn target_residues_all_is_created_correctly() {
        let targets = TargetResidues::All;
        assert_eq!(targets, TargetResidues::All);
    }

    #[test]
    fn target_residues_all_except_is_created_correctly() {
        let exceptions = vec![('A', 10), ('B', 25)];
        let targets = TargetResidues::AllExcept(exceptions.clone());
        assert_eq!(targets, TargetResidues::AllExcept(exceptions));
    }

    #[test]
    fn target_residues_explicit_is_created_correctly() {
        let explicit_list = vec![('C', 100)];
        let targets = TargetResidues::Explicit(explicit_list.clone());
        assert_eq!(targets, TargetResidues::Explicit(explicit_list));
    }

    #[test]
    fn design_task_can_be_constructed() {
        let mut positions = HashMap::new();
        positions.insert(('A', 5), vec![ResidueType::Alanine, ResidueType::Leucine]);
        positions.insert(('B', 20), vec![ResidueType::Tryptophan]);

        let task = DesignTask {
            positions: positions.clone(),
        };

        assert_eq!(task.positions.len(), 2);
        assert_eq!(task.positions.get(&('A', 5)).unwrap().len(), 2);
        assert_eq!(
            task.positions.get(&('B', 20)).unwrap(),
            &vec![ResidueType::Tryptophan]
        );
    }

    #[test]
    fn interaction_analysis_task_can_be_constructed() {
        let task = InteractionAnalysisTask {
            group1_selector: "chain A and resid 10".to_string(),
            group2_selector: "chain B and resname LIG".to_string(),
        };
        assert_eq!(task.group1_selector, "chain A and resid 10");
        assert_eq!(task.group2_selector, "chain B and resname LIG");
    }

    #[test]
    fn scream_task_can_be_place() {
        let task = ScreamTask::Place(TargetResidues::All);
        match task {
            ScreamTask::Place(TargetResidues::All) => (),
            _ => panic!("Expected ScreamTask::Place"),
        }
    }

    #[test]
    fn scream_task_can_be_design() {
        let design_task = DesignTask {
            positions: HashMap::new(),
        };
        let task = ScreamTask::Design(design_task);
        match task {
            ScreamTask::Design(_) => (),
            _ => panic!("Expected ScreamTask::Design"),
        }
    }

    #[test]
    fn scream_task_can_be_analyze() {
        let analysis_task = InteractionAnalysisTask {
            group1_selector: "group1".to_string(),
            group2_selector: "group2".to_string(),
        };
        let task = ScreamTask::Analyze(analysis_task);
        match task {
            ScreamTask::Analyze(_) => (),
            _ => panic!("Expected ScreamTask::Analyze"),
        }
    }

    #[test]
    fn algorithm_config_defaults_are_sensible() {
        let config = AlgorithmConfig::default();
        assert_eq!(config.s_factor, 1.0);
        assert_eq!(config.max_iterations, 10);
        assert_eq!(config.convergence_tolerance, 0.01);
    }
}
