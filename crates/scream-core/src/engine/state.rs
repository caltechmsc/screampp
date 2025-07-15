use crate::core::models::system::MolecularSystem;
use std::collections::BinaryHeap;

#[derive(Debug, Clone)]
pub struct Solution {
    pub energy: f64,
    pub system: MolecularSystem,
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
    pub working_system: MolecularSystem,
    pub current_energy: f64,
    top_solutions: BinaryHeap<Solution>,
    max_solutions: usize,
}

impl OptimizationState {
    pub fn new(initial_system: MolecularSystem, initial_energy: f64, max_solutions: usize) -> Self {
        let max_s = if max_solutions == 0 { 1 } else { max_solutions };

        let mut top_solutions = BinaryHeap::with_capacity(max_s);

        top_solutions.push(Solution {
            energy: initial_energy,
            system: initial_system.clone(),
        });

        Self {
            working_system: initial_system,
            current_energy: initial_energy,
            top_solutions,
            max_solutions: max_s,
        }
    }

    pub fn submit_current_solution(&mut self) {
        if self.top_solutions.len() < self.max_solutions {
            let new_solution = Solution {
                energy: self.current_energy,
                system: self.working_system.clone(),
            };
            self.top_solutions.push(new_solution);
        } else {
            let worst_of_the_best = self.top_solutions.peek().unwrap();

            if self.current_energy < worst_of_the_best.energy {
                let mut new_solution = Solution {
                    energy: self.current_energy,
                    system: self.working_system.clone(),
                };

                if let Some(mut placeholder) = self.top_solutions.pop() {
                    std::mem::swap(&mut placeholder.system, &mut new_solution.system);
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
