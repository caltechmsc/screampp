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
pub(crate) struct OptimizationState {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::ids::ResidueId;
    use crate::core::models::system::MolecularSystem;
    use slotmap::KeyData;
    use std::collections::HashMap;

    fn create_dummy_solution_state(residue_id_val: u64, rotamer_idx: usize) -> SolutionState {
        let mut rotamers = HashMap::new();
        let residue_id = ResidueId::from(KeyData::from_ffi(residue_id_val));
        rotamers.insert(residue_id, rotamer_idx);
        SolutionState {
            system: MolecularSystem::new(),
            rotamers,
        }
    }

    #[test]
    fn new_initializes_correctly() {
        let state = OptimizationState::new(MolecularSystem::new(), HashMap::new(), -100.0, 3);

        assert_eq!(state.max_solutions, 3);
        assert_eq!(state.solutions.len(), 1);
        assert_eq!(state.current_optimization_score, -100.0);
        assert_eq!(state.best_energy(), -100.0);
    }

    #[test]
    fn new_handles_zero_max_solutions() {
        let state = OptimizationState::new(MolecularSystem::new(), HashMap::new(), -100.0, 0);
        assert_eq!(state.max_solutions, 1);
        assert_eq!(state.solutions.len(), 1);
    }

    #[test]
    fn submit_solution_fills_heap_until_capacity() {
        let mut state = OptimizationState::new(MolecularSystem::new(), HashMap::new(), -10.0, 3);

        state.working_state = create_dummy_solution_state(1, 1);
        state.current_optimization_score = -20.0;
        state.submit_current_solution();

        state.working_state = create_dummy_solution_state(2, 2);
        state.current_optimization_score = -5.0;
        state.submit_current_solution();

        assert_eq!(state.solutions.len(), 3);
        assert_eq!(state.solutions.peek().unwrap().optimization_score, -5.0);
    }

    #[test]
    fn submit_better_solution_replaces_worst_when_full() {
        let mut state = OptimizationState::new(MolecularSystem::new(), HashMap::new(), -10.0, 2);
        state.current_optimization_score = -20.0;
        state.submit_current_solution();

        assert_eq!(state.best_energy(), -20.0);
        assert_eq!(state.solutions.peek().unwrap().optimization_score, -10.0);

        state.working_state = create_dummy_solution_state(3, 3);
        state.current_optimization_score = -30.0;
        state.submit_current_solution(); // Heap: [-20, -30]. Top is -20.

        assert_eq!(state.solutions.len(), 2);
        assert_eq!(state.best_energy(), -30.0);
        assert_eq!(state.solutions.peek().unwrap().optimization_score, -20.0);
    }

    #[test]
    fn submit_worse_solution_is_ignored_when_full() {
        let mut state = OptimizationState::new(MolecularSystem::new(), HashMap::new(), -20.0, 2);
        state.current_optimization_score = -30.0;
        state.submit_current_solution(); // Heap: [-20, -30] -> Top is -20

        state.working_state = create_dummy_solution_state(4, 4);
        state.current_optimization_score = -10.0;
        state.submit_current_solution();

        assert_eq!(state.solutions.len(), 2);
        let sorted_solutions = state.into_sorted_solutions();
        assert_eq!(sorted_solutions[0].optimization_score, -30.0);
        assert_eq!(sorted_solutions[1].optimization_score, -20.0);
    }

    #[test]
    fn submit_equal_solution_is_ignored_when_full() {
        let mut state = OptimizationState::new(MolecularSystem::new(), HashMap::new(), -20.0, 2);
        state.current_optimization_score = -30.0;
        state.submit_current_solution(); // Heap: [-20, -30] -> Top is -20

        state.working_state = create_dummy_solution_state(5, 5);
        state.current_optimization_score = -20.0;
        state.submit_current_solution();

        assert_eq!(state.solutions.len(), 2);
        let sorted_solutions = state.into_sorted_solutions();
        assert_eq!(sorted_solutions[0].optimization_score, -30.0);
        assert_eq!(sorted_solutions[1].optimization_score, -20.0);
    }

    #[test]
    fn into_sorted_solutions_returns_correctly_ordered_vec() {
        let mut state = OptimizationState::new(MolecularSystem::new(), HashMap::new(), -15.0, 5);
        state.current_optimization_score = -25.0;
        state.submit_current_solution();
        state.current_optimization_score = -5.0;
        state.submit_current_solution();
        state.current_optimization_score = -35.0;
        state.submit_current_solution();

        let sorted = state.into_sorted_solutions();
        assert_eq!(sorted.len(), 4);
        assert_eq!(sorted[0].optimization_score, -35.0);
        assert_eq!(sorted[1].optimization_score, -25.0);
        assert_eq!(sorted[2].optimization_score, -15.0);
        assert_eq!(sorted[3].optimization_score, -5.0);
    }

    #[test]
    fn best_energy_updates_correctly() {
        let mut state = OptimizationState::new(MolecularSystem::new(), HashMap::new(), -10.0, 3);
        assert_eq!(state.best_energy(), -10.0);

        state.current_optimization_score = -5.0;
        state.submit_current_solution();
        assert_eq!(state.best_energy(), -10.0);

        state.current_optimization_score = -20.0;
        state.submit_current_solution();
        assert_eq!(state.best_energy(), -20.0);
    }

    #[test]
    fn submitted_solution_state_is_a_deep_clone() {
        let mut state = OptimizationState::new(MolecularSystem::new(), HashMap::new(), -10.0, 2);

        state.working_state = create_dummy_solution_state(1, 1);
        state.current_optimization_score = -20.0;
        state.submit_current_solution();

        state.working_state = create_dummy_solution_state(2, 2);
        state.current_optimization_score = -30.0;

        let solutions = state.clone().into_sorted_solutions();
        let solution_minus_20 = solutions
            .iter()
            .find(|s| (s.optimization_score - (-20.0)).abs() < 1e-9)
            .unwrap();

        let original_residue_id = ResidueId::from(KeyData::from_ffi(1));
        assert_eq!(
            solution_minus_20.state.rotamers.get(&original_residue_id),
            Some(&1)
        );
        assert_eq!(solution_minus_20.state.rotamers.len(), 1);

        state.submit_current_solution();
        let final_solutions = state.into_sorted_solutions();
        assert_eq!(final_solutions[0].optimization_score, -30.0);
        assert_eq!(final_solutions[1].optimization_score, -20.0);

        let final_residue_id = ResidueId::from(KeyData::from_ffi(2));
        assert_eq!(
            final_solutions[0].state.rotamers.get(&final_residue_id),
            Some(&2)
        );
    }
}
