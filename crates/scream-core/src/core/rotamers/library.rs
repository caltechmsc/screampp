use super::rotamer::{Rotamer, RotamerData};
use crate::core::forcefield::parameterization::{ParameterizationError, Parameterizer};
use crate::core::forcefield::params::Forcefield;
use crate::core::models::atom::Atom;
use crate::core::models::ids::ResidueId;
use crate::core::models::residue::ResidueType;
use nalgebra::Point3;
use std::collections::HashMap;
use std::path::Path;
use std::str::FromStr;
use thiserror::Error;

type RawRotamerFile = HashMap<String, Vec<RotamerData>>;

#[derive(Debug, Default)]
pub struct RotamerLibrary {
    pub rotamers: HashMap<ResidueType, Vec<Rotamer>>,
}
