use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use std::collections::{BinaryHeap, HashMap};

/// Represents the initial state of a molecular system before optimization.
///
/// This struct captures the starting point of an optimization process, including
/// the molecular system configuration and its associated energy values. It serves
/// as a baseline for tracking optimization progress and comparing final results.
#[derive(Debug, Clone)]
pub struct InitialState {
    /// The molecular system in its initial configuration.
    pub system: MolecularSystem,
    /// The total energy of the initial system configuration.
    pub total_energy: f64,
    /// The optimization score used to evaluate the initial state.
    pub optimization_score: f64,
}

/// Represents a specific state of the molecular system with assigned rotamers.
///
/// This struct stores a snapshot of the molecular system along with the rotamer
/// assignments for each residue. It captures both the structural configuration
/// and the conformational choices made during optimization.
#[derive(Debug, Clone)]
pub struct SolutionState {
    /// The molecular system with its current atomic coordinates and structure.
    pub system: MolecularSystem,
    /// Mapping of residue IDs to their assigned rotamer indices.
    pub rotamers: HashMap<ResidueId, usize>,
}

/// Represents a complete solution from the optimization process.
///
/// This struct encapsulates a solution found during optimization, including its
/// energy values, optimization score, and the complete system state with rotamer
/// assignments. Solutions are comparable based on their optimization scores.
#[derive(Debug, Clone)]
pub struct Solution {
    /// The total energy of this solution configuration.
    pub total_energy: f64,
    /// The optimization score used for ranking and comparison of solutions.
    pub optimization_score: f64,
    /// The complete state of the molecular system and rotamer assignments.
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
    // Solutions are ordered by optimization score, with lower scores being better
    // This makes the BinaryHeap act as a max-heap where the "largest" (worst) solution is at the top
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap_or(std::cmp::Ordering::Equal)
    }
}

/// Manages the state of an ongoing optimization process.
///
/// This struct tracks the current working state of the molecular system during
/// optimization, maintains a collection of the best solutions found so far, and
/// provides methods for submitting new solutions and retrieving results. It uses
/// a binary heap to efficiently maintain the top N solutions by optimization score.
#[derive(Debug, Clone)]
pub(crate) struct OptimizationState {
    /// The current working state of the molecular system being optimized.
    pub working_state: SolutionState,
    /// The optimization score of the current working state.
    pub current_optimization_score: f64,
    /// Binary heap storing the best solutions found, ordered by optimization score.
    solutions: BinaryHeap<Solution>,
    /// Maximum number of solutions to maintain in the collection.
    max_solutions: usize,
}

impl OptimizationState {
    /// Creates a new optimization state with the initial system configuration.
    ///
    /// # Arguments
    ///
    /// * `initial_system` - The starting molecular system configuration.
    /// * `initial_rotamers` - Initial rotamer assignments for residues.
    /// * `initial_optimization_score` - Starting optimization score.
    /// * `max_solutions` - Maximum number of solutions to track (minimum 1).
    ///
    /// # Return
    ///
    /// Returns a new `OptimizationState` instance initialized with the provided parameters.
    pub fn new(
        initial_system: MolecularSystem,
        initial_rotamers: HashMap<ResidueId, usize>,
        initial_optimization_score: f64,
        max_solutions: usize,
    ) -> Self {
        // Ensure at least one solution is tracked to avoid empty collections
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

    /// Submits the current working state as a potential solution.
    ///
    /// This method evaluates the current working state and adds it to the solution
    /// collection if it improves upon the existing solutions. If the collection is
    /// at capacity, it replaces the worst solution if the current state is better.
    /// The optimization score determines solution quality, with lower scores being better.
    pub fn submit_current_solution(&mut self) {
        if self.solutions.len() < self.max_solutions {
            self.solutions.push(Solution {
                total_energy: f64::NAN,
                optimization_score: self.current_optimization_score,
                state: self.working_state.clone(),
            });
            return;
        }

        // BinaryHeap is a max-heap, so peek() gives the worst solution (highest score)
        // Replace it only if current solution has a better (lower) score
        if let Some(worst_of_the_best) = self.solutions.peek() {
            if self.current_optimization_score < worst_of_the_best.optimization_score {
                let mut worst_solution_placeholder = self.solutions.pop().unwrap();
                worst_solution_placeholder.optimization_score = self.current_optimization_score;
                worst_solution_placeholder.state = self.working_state.clone();
                self.solutions.push(worst_solution_placeholder);
            }
        }
    }

    /// Consumes the optimization state and returns all solutions in sorted order.
    ///
    /// Solutions are sorted by optimization score in ascending order (best first).
    /// This method takes ownership of the state, so it should be called when
    /// optimization is complete.
    ///
    /// # Return
    ///
    /// Returns a `Vec<Solution>` containing all tracked solutions, sorted by quality.
    pub fn into_sorted_solutions(self) -> Vec<Solution> {
        self.solutions.into_sorted_vec()
    }

    /// Returns the best optimization score among all tracked solutions.
    ///
    /// # Return
    ///
    /// Returns the minimum optimization score (best value) from all solutions,
    /// or `f64::INFINITY` if no solutions are available.
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
