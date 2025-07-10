use super::params::Forcefield;
use crate::core::models::atom::Atom;
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
    delta_s_factor: f64,
}

impl Parameterizer {
    pub fn new(forcefield: Forcefield, delta_s_factor: f64) -> Self {
        Self {
            forcefield,
            delta_s_factor,
        }
    }

    pub fn parameterize_system(
        &self,
        system: &mut MolecularSystem,
    ) -> Result<(), ParameterizationError> {
        self.parameterize_topology(system)?;
        self.parameterize_charges(system)?;
        self.parameterize_deltas(system)?;
        self.parameterize_non_bonded_properties(system)?;
        Ok(())
    }

    pub fn parameterize_topology(
        &self,
        system: &mut MolecularSystem,
    ) -> Result<(), ParameterizationError> {
        let mut ff_type_modifications = Vec::new();
        let mut bonds_to_add = Vec::new();

        let residue_ids: Vec<_> = system.residues_iter().map(|(id, _)| id).collect();

        for res_id in residue_ids {
            let residue = system.residue(res_id).unwrap();
            let res_type_str = &residue.name;

            let topology = self
                .forcefield
                .topology
                .get(res_type_str)
                .ok_or_else(|| ParameterizationError::MissingTopology(res_type_str.clone()))?;

            for atom_id in residue.atoms() {
                let atom_name = &system.atom(*atom_id).unwrap().name;
                let atom_topo = topology
                    .atoms
                    .iter()
                    .find(|a| &a.name == atom_name)
                    .ok_or_else(|| ParameterizationError::AtomNotFoundInTopology {
                        res_type: res_type_str.clone(),
                        atom_name: atom_name.clone(),
                    })?;
                ff_type_modifications.push((*atom_id, atom_topo.ff_type.clone()));
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

        for (atom_id, ff_type) in ff_type_modifications {
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

    pub fn parameterize_deltas(
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

            let final_delta = if res_type_str == "GLY" {
                0.0
            } else {
                let delta_param = self
                    .forcefield
                    .deltas
                    .get(&(res_type_str.clone(), atom_name.clone()));

                match delta_param {
                    Some(p) => p.mu + self.delta_s_factor * p.sigma,
                    None => 0.0,
                }
            };

            if let Some(atom) = system.atom_mut(atom_id) {
                atom.delta = final_delta;
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

    pub fn parameterize_atom(
        &self,
        atom: &mut Atom,
        res_type_str: &str,
    ) -> Result<(), ParameterizationError> {
        let topology = self
            .forcefield
            .topology
            .get(res_type_str)
            .ok_or_else(|| ParameterizationError::MissingTopology(res_type_str.to_string()))?;

        let atom_topo = topology
            .atoms
            .iter()
            .find(|a| a.name == atom.name)
            .ok_or_else(|| ParameterizationError::AtomNotFoundInTopology {
                res_type: res_type_str.to_string(),
                atom_name: atom.name.clone(),
            })?;

        atom.force_field_type = atom_topo.ff_type.clone();

        atom.delta = if res_type_str == "GLY" {
            0.0
        } else {
            self.forcefield
                .deltas
                .get(&(res_type_str.to_string(), atom.name.clone()))
                .map_or(0.0, |p| p.mu + self.delta_s_factor * p.sigma)
        };

        let ff_type = &atom.force_field_type;
        let vdw_param = self
            .forcefield
            .non_bonded
            .vdw
            .get(ff_type)
            .ok_or_else(|| ParameterizationError::MissingVdwParams(ff_type.clone()))?;

        let (vdw_radius, vdw_well_depth) = match vdw_param {
            super::params::VdwParam::LennardJones { radius, well_depth } => (*radius, *well_depth),
            super::params::VdwParam::Buckingham {
                radius, well_depth, ..
            } => (*radius, *well_depth),
        };
        atom.vdw_radius = vdw_radius;
        atom.vdw_well_depth = vdw_well_depth;

        atom.hbond_type_id = if ff_type == "H___A" {
            0
        } else {
            self.forcefield
                .non_bonded
                .hbond
                .get(ff_type)
                .map_or(-1, |_| 1)
        };

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::params::{
        ChargeParam, DeltaParam, Forcefield, GlobalParams, HBondParam, NonBondedParams,
        TopologyAtomParam, TopologyResidueParams, VdwParam,
    };
    use crate::core::models::atom::Atom;
    use crate::core::models::chain::ChainType;
    use crate::core::models::system::MolecularSystem;
    use nalgebra::Point3;
    use std::collections::HashMap;

    fn create_dummy_forcefield() -> Forcefield {
        let topology = {
            let mut map = HashMap::new();
            map.insert(
                "ALA".to_string(),
                TopologyResidueParams {
                    atoms: vec![
                        TopologyAtomParam {
                            name: "N".to_string(),
                            ff_type: "N_AM".to_string(),
                        },
                        TopologyAtomParam {
                            name: "CA".to_string(),
                            ff_type: "C_SP3".to_string(),
                        },
                        TopologyAtomParam {
                            name: "C".to_string(),
                            ff_type: "C_CARB".to_string(),
                        },
                    ],
                    bonds: vec![["N".to_string(), "CA".to_string()]],
                },
            );
            map.insert(
                "GLY".to_string(),
                TopologyResidueParams {
                    atoms: vec![TopologyAtomParam {
                        name: "CA".to_string(),
                        ff_type: "C_GLY".to_string(),
                    }],
                    bonds: vec![],
                },
            );
            map
        };

        let charges = {
            let mut map = HashMap::new();
            map.insert(
                ("ALA".to_string(), "N".to_string()),
                ChargeParam {
                    res_type: "ALA".to_string(),
                    atom_name: "N".to_string(),
                    partial_charge: -0.5,
                },
            );
            map.insert(
                ("ALA".to_string(), "CA".to_string()),
                ChargeParam {
                    res_type: "ALA".to_string(),
                    atom_name: "CA".to_string(),
                    partial_charge: 0.1,
                },
            );
            map.insert(
                ("ALA".to_string(), "C".to_string()),
                ChargeParam {
                    res_type: "ALA".to_string(),
                    atom_name: "C".to_string(),
                    partial_charge: 0.4,
                },
            );
            map
        };

        let deltas = {
            let mut map = HashMap::new();
            map.insert(
                ("ALA".to_string(), "CA".to_string()),
                DeltaParam {
                    res_type: "ALA".to_string(),
                    atom_name: "CA".to_string(),
                    mu: 1.23,
                    sigma: 0.0,
                },
            );
            map
        };

        let non_bonded = {
            let mut vdw = HashMap::new();
            vdw.insert(
                "N_AM".to_string(),
                VdwParam::LennardJones {
                    radius: 1.8,
                    well_depth: 0.2,
                },
            );
            vdw.insert(
                "C_SP3".to_string(),
                VdwParam::Buckingham {
                    radius: 2.0,
                    well_depth: 0.1,
                    scale: 12.0,
                },
            );
            vdw.insert(
                "C_CARB".to_string(),
                VdwParam::LennardJones {
                    radius: 1.9,
                    well_depth: 0.08,
                },
            );
            vdw.insert(
                "H___A".to_string(),
                VdwParam::LennardJones {
                    radius: 1.0,
                    well_depth: 0.01,
                },
            );
            let mut hbond = HashMap::new();
            hbond.insert(
                "N_AM".to_string(),
                HBondParam {
                    equilibrium_dist: 2.7,
                    well_depth: 5.0,
                },
            );
            NonBondedParams {
                globals: GlobalParams {
                    dielectric_constant: 1.0,
                    potential_function: "lj".to_string(),
                },
                vdw,
                hbond,
            }
        };

        Forcefield {
            topology,
            charges,
            deltas,
            non_bonded,
        }
    }

    fn create_test_system() -> MolecularSystem {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);

        let res1_id = system.add_residue(chain_id, 1, "ALA", None).unwrap();
        let atom1_n = Atom::new(1, "N", res1_id, Point3::origin());
        let atom1_ca = Atom::new(2, "CA", res1_id, Point3::origin());
        let atom1_c = Atom::new(3, "C", res1_id, Point3::origin());
        system.add_atom_to_residue(res1_id, atom1_n).unwrap();
        system.add_atom_to_residue(res1_id, atom1_ca).unwrap();
        system.add_atom_to_residue(res1_id, atom1_c).unwrap();

        let res2_id = system.add_residue(chain_id, 2, "ALA", None).unwrap();
        let atom2_n = Atom::new(4, "N", res2_id, Point3::origin());
        let atom2_ca = Atom::new(5, "CA", res2_id, Point3::origin());
        let atom2_c = Atom::new(6, "C", res2_id, Point3::origin());
        system.add_atom_to_residue(res2_id, atom2_n).unwrap();
        system.add_atom_to_residue(res2_id, atom2_ca).unwrap();
        system.add_atom_to_residue(res2_id, atom2_c).unwrap();

        system
    }

    #[test]
    fn parameterize_topology_assigns_ff_types_and_bonds() {
        let ff = create_dummy_forcefield();
        let mut system = create_test_system();
        let parameterizer = Parameterizer::new(ff, 0.0);

        parameterizer.parameterize_topology(&mut system).unwrap();

        let atom_n_id = system.find_atom_by_serial(1).unwrap();
        let atom_ca_id = system.find_atom_by_serial(2).unwrap();
        let atom_n = system.atom(atom_n_id).unwrap();
        let atom_ca = system.atom(atom_ca_id).unwrap();

        assert_eq!(atom_n.force_field_type, "N_AM");
        assert_eq!(atom_ca.force_field_type, "C_SP3");
        assert_eq!(system.bonds().len(), 3);
        assert!(
            system
                .bonds()
                .iter()
                .any(|b| b.contains(atom_n_id) && b.contains(atom_ca_id))
        );
    }

    #[test]
    fn parameterize_topology_creates_peptide_bonds() {
        let ff = create_dummy_forcefield();
        let mut system = create_test_system();
        let parameterizer = Parameterizer::new(ff, 0.0);

        parameterizer.parameterize_topology(&mut system).unwrap();

        let atom_c1_id = system.find_atom_by_serial(3).unwrap();
        let atom_n2_id = system.find_atom_by_serial(4).unwrap();

        assert!(
            system
                .bonds()
                .iter()
                .any(|b| b.contains(atom_c1_id) && b.contains(atom_n2_id))
        );
    }

    #[test]
    fn parameterize_topology_fails_for_missing_residue_topology() {
        let ff = create_dummy_forcefield();
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        system.add_residue(chain_id, 1, "UNK", None).unwrap();
        let parameterizer = Parameterizer::new(ff, 0.0);

        let result = parameterizer.parameterize_topology(&mut system);
        assert!(matches!(
            result,
            Err(ParameterizationError::MissingTopology(res_type)) if res_type == "UNK"
        ));
    }

    #[test]
    fn parameterize_topology_fails_if_system_atom_is_missing_from_topology_definition() {
        let ff = create_dummy_forcefield();
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_id = system.add_residue(chain_id, 1, "ALA", None).unwrap();

        system
            .add_atom_to_residue(res_id, Atom::new(1, "N", res_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(res_id, Atom::new(2, "CA", res_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(res_id, Atom::new(3, "C", res_id, Point3::origin()))
            .unwrap();

        system
            .add_atom_to_residue(res_id, Atom::new(4, "H_EXTRA", res_id, Point3::origin()))
            .unwrap();
    }

    #[test]
    fn parameterize_charges_assigns_partial_charges() {
        let ff = create_dummy_forcefield();
        let mut system = create_test_system();
        let parameterizer = Parameterizer::new(ff, 0.0);

        parameterizer.parameterize_charges(&mut system).unwrap();

        let atom_n = system.atom(system.find_atom_by_serial(1).unwrap()).unwrap();
        let atom_ca = system.atom(system.find_atom_by_serial(2).unwrap()).unwrap();

        assert_eq!(atom_n.partial_charge, -0.5);
        assert_eq!(atom_ca.partial_charge, 0.1);
    }

    #[test]
    fn parameterize_charges_fails_for_missing_charge_parameter() {
        let ff = create_dummy_forcefield();
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_id = system.add_residue(chain_id, 1, "GLY", None).unwrap();
        let atom = Atom::new(1, "CA", res_id, Point3::origin());
        system.add_atom_to_residue(res_id, atom).unwrap();
        let parameterizer = Parameterizer::new(ff, 0.0);

        let result = parameterizer.parameterize_charges(&mut system);
        assert!(matches!(
            result,
            Err(ParameterizationError::MissingCharge { res_type, atom_name })
            if res_type == "GLY" && atom_name == "CA"
        ));
    }

    #[test]
    fn parameterize_deltas_assigns_values_correctly() {
        let ff = create_dummy_forcefield();
        let mut system = create_test_system();
        let parameterizer = Parameterizer::new(ff, 0.0);

        parameterizer.parameterize_system(&mut system).unwrap();

        let atom_n = system.atom(system.find_atom_by_serial(1).unwrap()).unwrap();
        let atom_ca = system.atom(system.find_atom_by_serial(2).unwrap()).unwrap();

        assert_eq!(atom_n.delta, 0.0);
        assert_eq!(atom_ca.delta, 1.23);
    }

    #[test]
    fn parameterize_deltas_assigns_zero_for_glycine() {
        let ff = create_dummy_forcefield();
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_id = system.add_residue(chain_id, 1, "GLY", None).unwrap();
        let atom_ca = Atom::new(1, "CA", res_id, Point3::origin());
        system.add_atom_to_residue(res_id, atom_ca).unwrap();
        let parameterizer = Parameterizer::new(ff, 0.0);

        parameterizer.parameterize_topology(&mut system).unwrap();
        parameterizer.parameterize_deltas(&mut system).unwrap();

        let atom_ca_id = system.find_atom_by_serial(1).unwrap();
        let atom = system.atom(atom_ca_id).unwrap();

        assert_eq!(atom.delta, 0.0);
    }

    #[test]
    fn parameterize_non_bonded_properties_assigns_vdw_and_hbond_params() {
        let ff = create_dummy_forcefield();
        let mut system = create_test_system();
        let parameterizer = Parameterizer::new(ff, 0.0);

        parameterizer.parameterize_topology(&mut system).unwrap();
        parameterizer
            .parameterize_non_bonded_properties(&mut system)
            .unwrap();

        let atom_n = system.atom(system.find_atom_by_serial(1).unwrap()).unwrap();
        let atom_ca = system.atom(system.find_atom_by_serial(2).unwrap()).unwrap();

        assert_eq!(atom_n.vdw_radius, 1.8);
        assert_eq!(atom_n.vdw_well_depth, 0.2);
        assert_eq!(atom_n.hbond_type_id, 1);

        assert_eq!(atom_ca.vdw_radius, 2.0);
        assert_eq!(atom_ca.vdw_well_depth, 0.1);
        assert_eq!(atom_ca.hbond_type_id, -1);
    }

    #[test]
    fn parameterize_non_bonded_properties_fails_for_missing_vdw_params() {
        let ff = create_dummy_forcefield();
        let mut system = create_test_system();
        let parameterizer = Parameterizer::new(ff, 0.0);

        let atom_id = system.find_atom_by_serial(1).unwrap();
        system.atom_mut(atom_id).unwrap().force_field_type = "UNKNOWN_TYPE".to_string();

        let result = parameterizer.parameterize_non_bonded_properties(&mut system);
        assert!(matches!(
            result,
            Err(ParameterizationError::MissingVdwParams(ff_type)) if ff_type == "UNKNOWN_TYPE"
        ));
    }

    #[test]
    fn parameterize_non_bonded_properties_fails_if_ff_type_is_not_set() {
        let ff = create_dummy_forcefield();
        let mut system = create_test_system();
        let parameterizer = Parameterizer::new(ff, 0.0);

        let result = parameterizer.parameterize_non_bonded_properties(&mut system);
        assert!(matches!(
            result,
            Err(ParameterizationError::MissingForceFieldType { .. })
        ));
    }

    #[test]
    fn full_parameterize_system_succeeds_on_valid_system_and_ff() {
        let ff = create_dummy_forcefield();
        let mut system = create_test_system();
        let parameterizer = Parameterizer::new(ff, 0.0);

        let result = parameterizer.parameterize_system(&mut system);
        assert!(result.is_ok());

        let atom_ca = system.atom(system.find_atom_by_serial(2).unwrap()).unwrap();
        assert_eq!(atom_ca.force_field_type, "C_SP3");
        assert_eq!(atom_ca.partial_charge, 0.1);
        assert_eq!(atom_ca.delta, 1.23);
        assert_eq!(atom_ca.vdw_radius, 2.0);
        assert_eq!(system.bonds().len(), 3);
    }
}
