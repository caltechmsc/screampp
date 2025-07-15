use crate::core::models::system::MolecularSystem;
use std::collections::BinaryHeap;

#[derive(Debug, Clone)]
pub struct Solution {
    pub energy: f64,
    pub system: MolecularSystem,
}
