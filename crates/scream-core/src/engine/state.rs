use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use std::collections::{BinaryHeap, HashMap};

#[derive(Debug, Clone)]
pub struct InitialState {
    pub system: MolecularSystem,
    pub total_energy: f64,
    pub optimization_score: f64,
}

#[derive(Debug, Clone)]
pub struct SolutionState {
    pub system: MolecularSystem,
    pub rotamers: HashMap<ResidueId, usize>,
}

#[derive(Debug, Clone)]
pub struct Solution {
    pub total_energy: f64,
    pub optimization_score: f64,
    pub state: SolutionState,
}

impl PartialEq for Solution {
    fn eq(&self, other: &Self) -> bool {
        self.optimization_score == other.optimization_score
    }
}
impl Eq for Solution {}

impl PartialOrd for Solution {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.optimization_score
            .partial_cmp(&other.optimization_score)
    }
}

impl Ord for Solution {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap_or(std::cmp::Ordering::Equal)
    }
}

#[derive(Debug, Clone)]
pub struct OptimizationState {
    pub working_state: SolutionState,
    pub current_optimization_score: f64,
    solutions: BinaryHeap<Solution>,
    max_solutions: usize,
}

impl OptimizationState {
    pub fn new(
        initial_system: MolecularSystem,
        initial_rotamers: HashMap<ResidueId, usize>,
        initial_optimization_score: f64,
        max_solutions: usize,
    ) -> Self {
        let max_s = if max_solutions == 0 { 1 } else { max_solutions };

        let initial_state = SolutionState {
            system: initial_system,
            rotamers: initial_rotamers,
        };

        let mut solutions = BinaryHeap::with_capacity(max_s);
        solutions.push(Solution {
            total_energy: f64::NAN,
            optimization_score: initial_optimization_score,
            state: initial_state.clone(),
        });

        Self {
            working_state: initial_state,
            current_optimization_score: initial_optimization_score,
            solutions,
            max_solutions: max_s,
        }
    }

    pub fn submit_current_solution(&mut self) {
        if self.solutions.len() < self.max_solutions {
            self.solutions.push(Solution {
                total_energy: f64::NAN,
                optimization_score: self.current_optimization_score,
                state: self.working_state.clone(),
            });
            return;
        }

        if let Some(worst_of_the_best) = self.solutions.peek() {
            if self.current_optimization_score < worst_of_the_best.optimization_score {
                let mut worst_solution_placeholder = self.solutions.pop().unwrap();
                worst_solution_placeholder.optimization_score = self.current_optimization_score;
                worst_solution_placeholder.state = self.working_state.clone();
                self.solutions.push(worst_solution_placeholder);
            }
        }
    }

    pub fn into_sorted_solutions(self) -> Vec<Solution> {
        self.solutions.into_sorted_vec()
    }

    pub fn best_energy(&self) -> f64 {
        self.solutions
            .iter()
            .map(|s| s.optimization_score)
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or(f64::INFINITY)
    }
}
