use super::params::Forcefield;
use crate::core::{
    models::{
        atom::{Atom, AtomRole, CachedVdwParam},
        chain::ChainType,
        ids::{AtomId, ResidueId},
        system::MolecularSystem,
    },
    topology::registry::{ResidueTopology, TopologyRegistry},
};
use std::collections::{HashMap, HashSet};
use thiserror::Error;
use tracing::warn;

#[derive(Debug, Error, PartialEq, Eq)]
pub enum ParameterizationError {
    #[error(
        "Missing VDW parameter for force field type: '{ff_type}' in atom '{atom_name}' of residue {residue_name}"
    )]
    MissingVdwParams {
        ff_type: String,
        atom_name: String,
        residue_name: String,
    },
    #[error(
        "In residue '{residue_name}', a sidechain atom '{atom_name}' defined in the topology registry was not found."
    )]
    MissingSidechainAtom {
        residue_name: String,
        atom_name: String,
    },
    #[error(
        "In residue '{residue_name}', an anchor atom '{atom_name}' was incorrectly defined as a sidechain atom in the topology registry."
    )]
    AnchorAtomMisclassified {
        residue_name: String,
        atom_name: String,
    },
}

pub struct Parameterizer<'a> {
    forcefield: &'a Forcefield,
    topology_registry: &'a TopologyRegistry,
    delta_s_factor: f64,
}

impl<'a> Parameterizer<'a> {
    pub fn new(
        forcefield: &'a Forcefield,
        topology_registry: &'a TopologyRegistry,
        delta_s_factor: f64,
    ) -> Self {
        Self {
            forcefield,
            topology_registry,
            delta_s_factor,
        }
    }

    pub fn parameterize_system(
        &self,
        system: &mut MolecularSystem,
    ) -> Result<(), ParameterizationError> {
        let residue_ids: Vec<_> = system.residues_iter().map(|(id, _)| id).collect();

        // Pass 1: Assign Atom Roles for each residue based on topology.
        for residue_id in &residue_ids {
            self.assign_atom_roles_for_residue(*residue_id, system)?;
        }

        // Pass 2: Assign physicochemical parameters to each atom.
        for atom_id in system.atoms_iter().map(|(id, _)| id).collect::<Vec<_>>() {
            let residue_name = system
                .residue(system.atom(atom_id).unwrap().residue_id)
                .unwrap()
                .name
                .clone();
            let atom = system.atom_mut(atom_id).unwrap();
            self.assign_physicochemical_params(atom, &residue_name)?;
        }

        Ok(())
    }

    fn assign_atom_roles_for_residue(
        &self,
        residue_id: ResidueId,
        system: &mut MolecularSystem,
    ) -> Result<(), ParameterizationError> {
        let (chain_type, residue_name, atom_ids) = {
            let residue = system.residue(residue_id).unwrap();
            let chain = system.chain(residue.chain_id).unwrap();
            (
                chain.chain_type,
                residue.name.clone(),
                residue.atoms().to_vec(),
            )
        };

        match chain_type {
            ChainType::Protein | ChainType::DNA | ChainType::RNA => {
                if let Some(topology) = self.topology_registry.get(&residue_name) {
                    self.assign_protein_atom_roles(system, &atom_ids, &residue_name, topology)?
                } else {
                    warn!(
                        "Residue '{}' in a protein chain has no topology definition. Treating all its atoms as 'Other'.",
                        residue_name
                    );
                    for atom_id in atom_ids {
                        system.atom_mut(atom_id).unwrap().role = AtomRole::Other;
                    }
                }
            }
            ChainType::Ligand => {
                for atom_id in atom_ids {
                    system.atom_mut(atom_id).unwrap().role = AtomRole::Ligand;
                }
            }
            ChainType::Water => {
                for atom_id in atom_ids {
                    system.atom_mut(atom_id).unwrap().role = AtomRole::Water;
                }
            }
            ChainType::Other => {
                for atom_id in atom_ids {
                    system.atom_mut(atom_id).unwrap().role = AtomRole::Other;
                }
            }
        }
        Ok(())
    }

    fn assign_protein_atom_roles(
        &self,
        system: &mut MolecularSystem,
        atom_ids: &[AtomId],
        residue_name: &str,
        topology: &ResidueTopology,
    ) -> Result<(), ParameterizationError> {
        let mut atom_pool: HashMap<String, Vec<AtomId>> = HashMap::new();
        for &atom_id in atom_ids {
            let atom_name = system.atom(atom_id).unwrap().name.clone();
            atom_pool.entry(atom_name).or_default().push(atom_id);
        }

        for name in &topology.sidechain_atoms {
            let atom_id = atom_pool.get_mut(name).and_then(|ids| ids.pop()).ok_or(
                ParameterizationError::MissingSidechainAtom {
                    residue_name: residue_name.to_string(),
                    atom_name: name.clone(),
                },
            )?;
            system.atom_mut(atom_id).unwrap().role = AtomRole::Sidechain;
        }

        let mut backbone_atom_names = HashSet::new();
        for ids in atom_pool.values() {
            for &atom_id in ids {
                let atom = system.atom_mut(atom_id).unwrap();
                atom.role = AtomRole::Backbone;
                backbone_atom_names.insert(atom.name.clone());
            }
        }

        for name in &topology.anchor_atoms {
            if !backbone_atom_names.contains(name) {
                return Err(ParameterizationError::AnchorAtomMisclassified {
                    residue_name: residue_name.to_string(),
                    atom_name: name.clone(),
                });
            }
        }

        Ok(())
    }

    fn assign_physicochemical_params(
        &self,
        atom: &mut Atom,
        residue_name: &str,
    ) -> Result<(), ParameterizationError> {
        let delta_param = self
            .forcefield
            .deltas
            .get(&(residue_name.to_string(), atom.name.clone()));
        atom.delta = match delta_param {
            Some(p) => p.mu + self.delta_s_factor * p.sigma,
            None => 0.0,
        };

        let ff_type = &atom.force_field_type;
        if ff_type.is_empty() {
            return Ok(());
        }

        let vdw_param = self.forcefield.non_bonded.vdw.get(ff_type).ok_or_else(|| {
            ParameterizationError::MissingVdwParams {
                ff_type: ff_type.clone(),
                atom_name: atom.name.clone(),
                residue_name: residue_name.to_string(),
            }
        })?;

        atom.vdw_param = (*vdw_param).clone().into();

        atom.hbond_type_id = if self.forcefield.non_bonded.hbond.contains_key(ff_type) {
            if ff_type.starts_with('H') {
                0 // Donor Hydrogen
            } else {
                1 // Acceptor
            }
        } else {
            -1 // Not involved in H-bonding
        };

        Ok(())
    }

    pub fn parameterize_protein_atom(
        &self,
        atom: &mut Atom,
        residue_name: &str,
        topology: &ResidueTopology,
    ) -> Result<(), ParameterizationError> {
        // Step 1: Assign Role
        if topology.sidechain_atoms.contains(&atom.name) {
            atom.role = AtomRole::Sidechain;
        } else {
            atom.role = AtomRole::Backbone;
        }

        // Step 2: Assign physicochemical parameters
        self.assign_physicochemical_params(atom, residue_name)?;

        Ok(())
    }
}

impl From<super::params::VdwParam> for CachedVdwParam {
    fn from(param: super::params::VdwParam) -> Self {
        match param {
            super::params::VdwParam::LennardJones { radius, well_depth } => {
                CachedVdwParam::LennardJones { radius, well_depth }
            }
            super::params::VdwParam::Buckingham {
                radius,
                well_depth,
                scale,
            } => CachedVdwParam::Buckingham {
                radius,
                well_depth,
                scale,
            },
        }
    }
}
