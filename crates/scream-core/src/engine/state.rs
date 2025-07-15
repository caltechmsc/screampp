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
        other.energy.partial_cmp(&self.energy)
    }
}

impl Ord for Solution {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap_or(std::cmp::Ordering::Equal)
    }
}
