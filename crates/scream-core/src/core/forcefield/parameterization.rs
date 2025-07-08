use crate::core::forcefield::params::Forcefield;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::core::models::topology::BondOrder;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ParameterizationError {
    #[error("Missing topology for residue type: {0}")]
    MissingTopology(String),
    #[error("Missing force field type for atom '{atom_name}' in residue '{res_type}'")]
    MissingForceFieldType { res_type: String, atom_name: String },
    #[error("Missing charge parameter for atom '{atom_name}' in residue '{res_type}'")]
    MissingCharge { res_type: String, atom_name: String },
    #[error("Missing VDW parameter for force field type: {0}")]
    MissingVdwParams(String),
    #[error("Missing delta parameter for atom '{atom_name}' in residue '{res_type}'")]
    MissingDelta { res_type: String, atom_name: String },
    #[error(
        "Inconsistent topology: Atom '{atom_name}' from input file not found in topology for residue '{res_type}'"
    )]
    AtomNotFoundInTopology { res_type: String, atom_name: String },
    #[error("Atom with name '{0}' not found in residue ID {1:?}")]
    AtomNameNotFound(String, ResidueId),
}

pub struct Parameterizer {
    forcefield: Forcefield,
}

impl Parameterizer {
    pub fn new(forcefield: Forcefield) -> Self {
        Self { forcefield }
    }

    pub fn parameterize_system(
        &self,
        system: &mut MolecularSystem,
    ) -> Result<(), ParameterizationError> {
        self.parameterize_topology(system)?;
        self.parameterize_charges(system)?;
        self.parameterize_non_bonded_properties(system)?;
        Ok(())
    }

    pub fn parameterize_topology(
        &self,
        system: &mut MolecularSystem,
    ) -> Result<(), ParameterizationError> {
        let mut modifications = Vec::new();
        let mut bonds_to_add = Vec::new();

        let residue_ids: Vec<_> = system.residues_iter().map(|(id, _)| id).collect();

        for res_id in &residue_ids {
            let res_id = *res_id;
            let residue = system.residue(res_id).unwrap();
            let res_type_str = &residue.name;

            let topology = self
                .forcefield
                .topology
                .get(res_type_str)
                .ok_or_else(|| ParameterizationError::MissingTopology(res_type_str.clone()))?;

            for atom_topo in &topology.atoms {
                if let Some(atom_id) = residue.get_atom_id_by_name(&atom_topo.name) {
                    modifications.push((atom_id, atom_topo.ff_type.clone()));
                } else {
                    return Err(ParameterizationError::AtomNotFoundInTopology {
                        res_type: res_type_str.clone(),
                        atom_name: atom_topo.name.clone(),
                    });
                }
            }

            for bond_pair in &topology.bonds {
                let atom1_id = residue.get_atom_id_by_name(&bond_pair[0]).ok_or_else(|| {
                    ParameterizationError::AtomNameNotFound(bond_pair[0].clone(), res_id)
                })?;
                let atom2_id = residue.get_atom_id_by_name(&bond_pair[1]).ok_or_else(|| {
                    ParameterizationError::AtomNameNotFound(bond_pair[1].clone(), res_id)
                })?;
                bonds_to_add.push((atom1_id, atom2_id));
            }
        }

        for (atom_id, ff_type) in modifications {
            if let Some(atom) = system.atom_mut(atom_id) {
                atom.force_field_type = ff_type;
            }
        }

        for (atom1_id, atom2_id) in bonds_to_add {
            let _ = system.add_bond(atom1_id, atom2_id, BondOrder::default());
        }

        self.build_peptide_bonds(system);

        Ok(())
    }

    pub fn parameterize_charges(
        &self,
        system: &mut MolecularSystem,
    ) -> Result<(), ParameterizationError> {
        let atom_ids: Vec<_> = system.atoms_iter().map(|(id, _)| id).collect();
        for atom_id in atom_ids {
            let (res_type_str, atom_name) = {
                let atom = system.atom(atom_id).unwrap();
                let residue = system.residue(atom.residue_id).unwrap();
                (residue.name.clone(), atom.name.clone())
            };

            let charge_param = self
                .forcefield
                .charges
                .get(&(res_type_str.clone(), atom_name.clone()))
                .ok_or_else(|| ParameterizationError::MissingCharge {
                    res_type: res_type_str,
                    atom_name,
                })?;

            if let Some(atom) = system.atom_mut(atom_id) {
                atom.partial_charge = charge_param.partial_charge;
            }
        }
        Ok(())
    }

    pub fn parameterize_non_bonded_properties(
        &self,
        system: &mut MolecularSystem,
    ) -> Result<(), ParameterizationError> {
        let atom_ids: Vec<_> = system.atoms_iter().map(|(id, _)| id).collect();
        for atom_id in atom_ids {
            let ff_type = system.atom(atom_id).unwrap().force_field_type.clone();

            if ff_type.is_empty() {
                return Err(ParameterizationError::MissingForceFieldType {
                    res_type: system
                        .residue(system.atom(atom_id).unwrap().residue_id)
                        .unwrap()
                        .name
                        .clone(),
                    atom_name: system.atom(atom_id).unwrap().name.clone(),
                });
            }

            let vdw_param = self
                .forcefield
                .non_bonded
                .vdw
                .get(&ff_type)
                .ok_or_else(|| ParameterizationError::MissingVdwParams(ff_type.clone()))?;

            let (vdw_radius, vdw_well_depth) = match vdw_param {
                super::params::VdwParam::LennardJones { radius, well_depth } => {
                    (*radius, *well_depth)
                }
                super::params::VdwParam::Buckingham {
                    radius, well_depth, ..
                } => (*radius, *well_depth),
            };

            let hbond_type_id = if ff_type == "H___A" {
                0
            } else {
                self.forcefield
                    .non_bonded
                    .hbond
                    .get(&ff_type)
                    .map_or(-1, |_| 1)
            };

            if let Some(atom) = system.atom_mut(atom_id) {
                atom.vdw_radius = vdw_radius;
                atom.vdw_well_depth = vdw_well_depth;
                atom.hbond_type_id = hbond_type_id;
            }
        }
        Ok(())
    }

    pub fn parameterize_deltas(
        &self,
        system: &mut MolecularSystem,
        _library_id: &str,
    ) -> Result<(), ParameterizationError> {
        let atom_ids: Vec<_> = system.atoms_iter().map(|(id, _)| id).collect();
        for atom_id in atom_ids {
            let (res_type_str, atom_name) = {
                let atom = system.atom(atom_id).unwrap();
                let residue = system.residue(atom.residue_id).unwrap();
                (residue.name.clone(), atom.name.clone())
            };

            let delta = if res_type_str == "GLY" {
                0.0
            } else {
                self.forcefield
                    .deltas
                    .get(&(res_type_str.clone(), atom_name.clone()))
                    .map_or(0.0, |p| p.mu)
            };

            if let Some(atom) = system.atom_mut(atom_id) {
                atom.delta = delta;
            }
        }
        Ok(())
    }

    fn build_peptide_bonds(&self, system: &mut MolecularSystem) {
        let mut peptide_bonds_to_add = Vec::new();

        for (_chain_id, chain) in system.chains_iter() {
            for residue_pair in chain.residues().windows(2) {
                let res1_id = residue_pair[0];
                let res2_id = residue_pair[1];

                let res1 = system.residue(res1_id).unwrap();
                let res2 = system.residue(res2_id).unwrap();

                if let (Some(c_id), Some(n_id)) =
                    (res1.get_atom_id_by_name("C"), res2.get_atom_id_by_name("N"))
                {
                    peptide_bonds_to_add.push((c_id, n_id));
                }
            }
        }

        for (atom1_id, atom2_id) in peptide_bonds_to_add {
            let _ = system.add_bond(atom1_id, atom2_id, BondOrder::Single);
        }
    }
}
