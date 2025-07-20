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
    pub energy: f64,
    pub state: SolutionState,
}

impl PartialEq for Solution {
    fn eq(&self, other: &Self) -> bool {
        self.energy == other.energy
    }
}
impl Eq for Solution {}

impl PartialOrd for Solution {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.energy.partial_cmp(&other.energy)
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
    pub current_energy: f64,
    top_solutions: BinaryHeap<Solution>,
    max_solutions: usize,
}

impl OptimizationState {
    pub fn new(
        initial_system: MolecularSystem,
        initial_rotamers: HashMap<ResidueId, usize>,
        initial_energy: f64,
        max_solutions: usize,
    ) -> Self {
        let max_s = if max_solutions == 0 { 1 } else { max_solutions };

        let initial_state = SolutionState {
            system: initial_system,
            rotamers: initial_rotamers,
        };

        let mut top_solutions = BinaryHeap::with_capacity(max_s);
        top_solutions.push(Solution {
            energy: initial_energy,
            state: initial_state.clone(),
        });

        Self {
            working_state: initial_state,
            current_energy: initial_energy,
            top_solutions,
            max_solutions: max_s,
        }
    }

    pub fn submit_current_solution(&mut self) {
        if self.top_solutions.len() < self.max_solutions {
            let new_solution = Solution {
                energy: self.current_energy,
                state: self.working_state.clone(),
            };
            self.top_solutions.push(new_solution);
        } else {
            let worst_of_the_best = self.top_solutions.peek().unwrap();

            if self.current_energy < worst_of_the_best.energy {
                let mut new_solution = Solution {
                    energy: self.current_energy,
                    state: self.working_state.clone(),
                };

                if let Some(mut placeholder) = self.top_solutions.pop() {
                    std::mem::swap(&mut placeholder.state, &mut new_solution.state);
                    placeholder.energy = new_solution.energy;
                    self.top_solutions.push(placeholder);
                }
            }
        }
    }

    pub fn into_sorted_solutions(self) -> Vec<Solution> {
        let mut solutions = self.top_solutions.into_vec();
        solutions.sort_unstable();
        solutions
    }

    pub fn best_solution(&self) -> Option<&Solution> {
        self.top_solutions.iter().min()
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
