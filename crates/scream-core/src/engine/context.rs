use super::config::ResidueSelection;
use super::error::EngineError;
use crate::core::forcefield::params::Forcefield;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::core::rotamers::library::RotamerLibrary;
use std::collections::HashSet;

#[derive(Clone)]
pub struct Context<'a, C> {
    pub system: &'a MolecularSystem,
    pub forcefield: &'a Forcefield,
    pub config: &'a C,
}

impl<'a, C> Context<'a, C> {
    pub fn new(system: &'a MolecularSystem, forcefield: &'a Forcefield, config: &'a C) -> Self {
        Self {
            system,
            forcefield,
            config,
        }
    }
}
