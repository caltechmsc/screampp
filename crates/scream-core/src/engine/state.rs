use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use std::collections::{BinaryHeap, HashMap};

#[derive(Debug, Clone)]
pub struct SolutionState {
    pub system: MolecularSystem,
    pub rotamers: HashMap<ResidueId, usize>,
}

#[derive(Debug, Clone)]
pub struct Solution {
    pub total_energy: f64, // TotalSystemEnergy = optimization_score + energy_offset_constant
    pub optimization_score: f64, // OptimizationScore = ΣE_EL(Sc_A) + E_inter(Sc_A, Sc_A)
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
    pub current_total_energy: f64,
    pub current_optimization_score: f64,
    top_solutions: BinaryHeap<Solution>,
    max_solutions: usize,
}

impl OptimizationState {
    pub fn new(
        initial_system: MolecularSystem,
        initial_rotamers: HashMap<ResidueId, usize>,
        initial_total_energy: f64,
        initial_optimization_score: f64,
        max_solutions: usize,
    ) -> Self {
        let max_s = if max_solutions == 0 { 1 } else { max_solutions };

        let initial_state = SolutionState {
            system: initial_system,
            rotamers: initial_rotamers,
        };

        let mut top_solutions = BinaryHeap::with_capacity(max_s);
        top_solutions.push(Solution {
            total_energy: initial_total_energy,
            optimization_score: initial_optimization_score,
            state: initial_state.clone(),
        });

        Self {
            working_state: initial_state,
            current_total_energy: initial_total_energy,
            current_optimization_score: initial_optimization_score,
            top_solutions,
            max_solutions: max_s,
        }
    }

    pub fn submit_current_solution(&mut self) {
        if self.top_solutions.len() < self.max_solutions {
            self.top_solutions.push(Solution {
                total_energy: self.current_total_energy,
                optimization_score: self.current_optimization_score,
                state: self.working_state.clone(),
            });
            return;
        }

        if let Some(worst_of_the_best) = self.top_solutions.peek() {
            if self.current_optimization_score < worst_of_the_best.optimization_score {
                if let Some(mut placeholder) = self.top_solutions.pop() {
                    placeholder.total_energy = self.current_total_energy;
                    placeholder.optimization_score = self.current_optimization_score;

                    std::mem::swap(&mut placeholder.state, &mut self.working_state);
                    self.top_solutions.push(placeholder);

                    if let Some(mut new_worst) = self.top_solutions.peek_mut() {
                        std::mem::swap(&mut new_worst.state, &mut self.working_state);
                    }
                }
            }
        }
    }

    pub fn into_sorted_solutions(self) -> Vec<Solution> {
        let mut solutions = self.top_solutions.into_vec();
        solutions.sort_unstable_by(|a, b| {
            a.optimization_score
                .partial_cmp(&b.optimization_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        solutions
    }

    pub fn best_solution(&self) -> Option<&Solution> {
        self.top_solutions.iter().min_by(|a, b| {
            a.optimization_score
                .partial_cmp(&b.optimization_score)
                .unwrap()
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn state_initializes_correctly_with_one_solution() {
        let system = MolecularSystem::new();
        let rotamers = HashMap::new();
        let state = OptimizationState::new(system, rotamers, -50.0, 5);

        assert_eq!(state.current_energy, -50.0);
        assert_eq!(state.max_solutions, 5);
        assert_eq!(state.top_solutions.len(), 1);
        assert_eq!(state.best_solution().unwrap().energy, -50.0);
    }

    #[test]
    fn state_handles_zero_max_solutions_by_setting_it_to_one() {
        let system = MolecularSystem::new();
        let rotamers = HashMap::new();
        let state = OptimizationState::new(system, rotamers, -50.0, 0);
        assert_eq!(state.max_solutions, 1);
        assert_eq!(state.top_solutions.len(), 1);
    }

    #[test]
    fn submit_better_solution_replaces_worst_when_heap_is_full() {
        let system = MolecularSystem::new();
        let rotamers = HashMap::new();
        let mut state = OptimizationState::new(system, rotamers, -50.0, 1);

        state.current_energy = -100.0;
        state.submit_current_solution();

        assert_eq!(state.top_solutions.len(), 1);
        assert_eq!(state.best_solution().unwrap().energy, -100.0);
    }

    #[test]
    fn submit_worse_solution_is_ignored_when_heap_is_full() {
        let system = MolecularSystem::new();
        let rotamers = HashMap::new();
        let mut state = OptimizationState::new(system, rotamers, -100.0, 1);

        state.current_energy = -50.0;
        state.submit_current_solution();

        assert_eq!(state.top_solutions.len(), 1);
        assert_eq!(state.best_solution().unwrap().energy, -100.0);
    }

    #[test]
    fn submit_solution_fills_heap_until_capacity() {
        let system = MolecularSystem::new();
        let rotamers = HashMap::new();
        let mut state = OptimizationState::new(system.clone(), rotamers.clone(), -10.0, 3);

        state.current_energy = -20.0;
        state.working_state.system = system.clone();
        state.submit_current_solution();

        state.current_energy = -5.0;
        state.working_state.system = system.clone();
        state.submit_current_solution();

        assert_eq!(state.top_solutions.len(), 3);
        let solutions = state.clone().into_sorted_solutions();
        assert_eq!(solutions[0].energy, -20.0);
        assert_eq!(solutions[1].energy, -10.0);
        assert_eq!(solutions[2].energy, -5.0);
    }

    #[test]
    fn submit_solutions_maintains_correct_top_n() {
        let system = MolecularSystem::new();
        let rotamers = HashMap::new();
        let mut state = OptimizationState::new(system.clone(), rotamers.clone(), -10.0, 2);

        state.current_energy = -20.0;
        state.working_state.system = system.clone();
        state.submit_current_solution(); // Heap: [-20, -10]

        state.current_energy = -5.0;
        state.working_state.system = system.clone();
        state.submit_current_solution(); // Heap: [-20, -10], -5 is worse than -10, ignored

        state.current_energy = -30.0;
        state.working_state.system = system.clone();
        state.submit_current_solution(); // Heap: [-30, -20], -10 is popped

        let solutions = state.into_sorted_solutions();
        assert_eq!(solutions.len(), 2);
        assert_eq!(solutions[0].energy, -30.0);
        assert_eq!(solutions[1].energy, -20.0);
    }

    #[test]
    fn into_sorted_solutions_returns_correctly_ordered_list() {
        let system = MolecularSystem::new();
        let rotamers = HashMap::new();
        let mut state = OptimizationState::new(system.clone(), rotamers, -10.0, 5);

        state.current_energy = -30.0;
        state.submit_current_solution();
        state.current_energy = -20.0;
        state.submit_current_solution();

        let solutions = state.into_sorted_solutions();
        assert_eq!(solutions.len(), 3);
        assert_eq!(solutions[0].energy, -30.0);
        assert_eq!(solutions[1].energy, -20.0);
        assert_eq!(solutions[2].energy, -10.0);
    }

    #[test]
    fn best_solution_returns_minimum_energy_solution() {
        let system = MolecularSystem::new();
        let rotamers = HashMap::new();
        let mut state = OptimizationState::new(system.clone(), rotamers, -10.0, 5);

        state.current_energy = -30.0;
        state.submit_current_solution();
        state.current_energy = -5.0;
        state.submit_current_solution();

        let best = state.best_solution().unwrap();
        assert_eq!(best.energy, -30.0);
    }

    #[test]
    fn submit_solution_with_equal_energy_is_ignored_when_full() {
        let system = MolecularSystem::new();
        let rotamers = HashMap::new();
        let mut state = OptimizationState::new(system.clone(), rotamers, -10.0, 2);
        state.current_energy = -20.0;
        state.submit_current_solution(); // Heap: [-20, -10]

        // Submit with energy equal to the worst in the heap
        state.current_energy = -10.0;
        state.submit_current_solution();

        let solutions = state.into_sorted_solutions();
        assert_eq!(solutions.len(), 2);
        assert_eq!(solutions[0].energy, -20.0);
        assert_eq!(solutions[1].energy, -10.0);
    }
}
