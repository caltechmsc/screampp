use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Deserialize, Clone, PartialEq)]
#[serde(untagged)]
pub enum VdwParam {
    Buckingham {
        radius: f64,
        well_depth: f64,
        scale: f64,
    },
    LennardJones {
        radius: f64,
        well_depth: f64,
    },
}

#[derive(Debug, Deserialize, Clone, PartialEq)]
pub struct HBondParam {
    pub equilibrium_dist: f64,
    pub well_depth: f64,
}

#[derive(Debug, Deserialize, Clone)]
pub struct DeltaParam {
    pub residue_type: String,
    pub atom_name: String,
    pub mu: f64,
    pub sigma: f64,
}

#[derive(Debug, Deserialize, Clone)]
pub struct ChargeParam {
    pub residue_type: String,
    pub atom_name: String,
    pub partial_charge: f64,
}

#[derive(Debug, Deserialize, Clone)]
pub struct TopologyAtomParam {
    pub name: String,
    pub ff_type: String,
}

#[derive(Debug, Deserialize, Clone)]
pub struct TopologyResidueParams {
    pub atoms: Vec<TopologyAtomParam>,
    pub bonds: Vec<[String; 2]>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct TopologyParams {
    #[serde(flatten)]
    pub residues: HashMap<String, TopologyResidueParams>,
}
