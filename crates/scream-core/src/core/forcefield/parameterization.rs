use super::params::Forcefield;
use crate::core::models::atom::{Atom, CachedVdwParam};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ParameterizationError {
    #[error("Missing VDW parameter for force field type: {0}")]
    MissingVdwParams(String),
    #[error("Missing delta parameter for atom '{atom_name}' in residue '{residue_type}'")]
    MissingDelta {
        residue_type: String,
        atom_name: String,
    },
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
        system: &mut crate::core::models::system::MolecularSystem,
    ) -> Result<(), ParameterizationError> {
        self.parameterize_deltas(system)?;
        self.parameterize_non_bonded_properties(system)?;
        Ok(())
    }

    pub fn parameterize_deltas(
        &self,
        system: &mut crate::core::models::system::MolecularSystem,
    ) -> Result<(), ParameterizationError> {
        let atom_ids: Vec<_> = system.atoms_iter().map(|(id, _)| id).collect();

        for atom_id in atom_ids {
            let (residue_type_str, atom_name) = {
                let atom = system.atom(atom_id).unwrap();
                let residue = system.residue(atom.residue_id).unwrap();
                (residue.name.clone(), atom.name.clone())
            };

            let final_delta = if residue_type_str == "GLY" {
                0.0
            } else {
                let delta_param = self
                    .forcefield
                    .deltas
                    .get(&(residue_type_str.clone(), atom_name.clone()));

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
        system: &mut crate::core::models::system::MolecularSystem,
    ) -> Result<(), ParameterizationError> {
        let atom_ids: Vec<_> = system.atoms_iter().map(|(id, _)| id).collect();
        for atom_id in atom_ids {
            let ff_type = system.atom(atom_id).unwrap().force_field_type.clone();

            if ff_type.is_empty() {
                return Err(ParameterizationError::MissingVdwParams(ff_type));
            }

            let vdw_param = self
                .forcefield
                .non_bonded
                .vdw
                .get(&ff_type)
                .ok_or_else(|| ParameterizationError::MissingVdwParams(ff_type.clone()))?;

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
                atom.vdw_param = match vdw_param {
                    super::params::VdwParam::LennardJones { radius, well_depth } => {
                        CachedVdwParam::LennardJones {
                            radius: *radius,
                            well_depth: *well_depth,
                        }
                    }
                    super::params::VdwParam::Buckingham {
                        radius,
                        well_depth,
                        scale,
                    } => CachedVdwParam::Buckingham {
                        radius: *radius,
                        well_depth: *well_depth,
                        scale: *scale,
                    },
                };
                atom.hbond_type_id = hbond_type_id;
            }
        }
        Ok(())
    }

    pub fn parameterize_atom(
        &self,
        atom: &mut Atom,
        residue_type_str: &str,
    ) -> Result<(), ParameterizationError> {
        atom.delta = if residue_type_str == "GLY" {
            0.0
        } else {
            self.forcefield
                .deltas
                .get(&(residue_type_str.to_string(), atom.name.clone()))
                .map_or(0.0, |p| p.mu + self.delta_s_factor * p.sigma)
        };

        let ff_type = &atom.force_field_type;
        let vdw_param = self
            .forcefield
            .non_bonded
            .vdw
            .get(ff_type)
            .ok_or_else(|| ParameterizationError::MissingVdwParams(ff_type.clone()))?;

        atom.vdw_param = match vdw_param {
            super::params::VdwParam::LennardJones { radius, well_depth } => {
                CachedVdwParam::LennardJones {
                    radius: *radius,
                    well_depth: *well_depth,
                }
            }
            super::params::VdwParam::Buckingham {
                radius,
                well_depth,
                scale,
            } => CachedVdwParam::Buckingham {
                radius: *radius,
                well_depth: *well_depth,
                scale: *scale,
            },
        };

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
        DeltaParam, Forcefield, GlobalParams, NonBondedParams, VdwParam,
    };
    use crate::core::models::atom::Atom;
    use crate::core::models::chain::ChainType;
    use crate::core::models::ids::AtomId;
    use crate::core::models::system::MolecularSystem;
    use nalgebra::Point3;
    use std::collections::HashMap;

    fn create_dummy_forcefield() -> Forcefield {
        let deltas = {
            let mut map = HashMap::new();
            map.insert(
                ("ALA".to_string(), "CA".to_string()),
                DeltaParam {
                    residue_type: "ALA".to_string(),
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
                "C_SP3".to_string(),
                VdwParam::Buckingham {
                    radius: 2.0,
                    well_depth: 0.1,
                    scale: 12.0,
                },
            );
            vdw.insert(
                "H___A".to_string(),
                VdwParam::LennardJones {
                    radius: 1.0,
                    well_depth: 0.01,
                },
            );
            let hbond = HashMap::new();
            NonBondedParams {
                globals: GlobalParams {
                    dielectric_constant: 1.0,
                    potential_function: "lennard-jones-12-6".to_string(),
                },
                vdw,
                hbond,
            }
        };

        Forcefield { deltas, non_bonded }
    }

    fn create_test_system() -> (MolecularSystem, AtomId) {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res1_id = system.add_residue(chain_id, 1, "ALA", None).unwrap();
        let mut atom1_ca = Atom::new("CA", res1_id, Point3::origin());
        atom1_ca.force_field_type = "C_SP3".to_string();
        let atom_id = system.add_atom_to_residue(res1_id, atom1_ca).unwrap();
        (system, atom_id)
    }

    fn create_dummy_forcefield_for_atom_test() -> Forcefield {
        let mut ff = create_dummy_forcefield();
        ff.deltas.insert(
            ("ALA".to_string(), "CA".to_string()),
            DeltaParam {
                residue_type: "ALA".to_string(),
                atom_name: "CA".to_string(),
                mu: 1.23,
                sigma: 0.2,
            },
        );
        ff
    }

    #[test]
    fn parameterize_deltas_assigns_values_correctly() {
        let ff = create_dummy_forcefield();
        let (mut system, atom_ca_id) = create_test_system();
        let parameterizer = Parameterizer::new(ff, 0.0);

        parameterizer.parameterize_system(&mut system).unwrap();

        let atom_ca = system.atom(atom_ca_id).unwrap();
        assert_eq!(atom_ca.delta, 1.23);
    }

    #[test]
    fn parameterize_deltas_assigns_zero_for_glycine() {
        let ff = create_dummy_forcefield();
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_id = system.add_residue(chain_id, 1, "GLY", None).unwrap();
        let mut atom_ca = Atom::new("CA", res_id, Point3::origin());
        atom_ca.force_field_type = "C_SP3".to_string();
        let atom_ca_id = system.add_atom_to_residue(res_id, atom_ca).unwrap();
        let parameterizer = Parameterizer::new(ff, 0.0);

        parameterizer.parameterize_deltas(&mut system).unwrap();

        let atom = system.atom(atom_ca_id).unwrap();
        assert_eq!(atom.delta, 0.0);
    }

    #[test]
    fn parameterize_non_bonded_properties_assigns_vdw_and_hbond_params() {
        let ff = create_dummy_forcefield();
        let (mut system, atom_ca_id) = create_test_system();
        let parameterizer = Parameterizer::new(ff, 0.0);

        parameterizer
            .parameterize_non_bonded_properties(&mut system)
            .unwrap();

        let atom_ca = system.atom(atom_ca_id).unwrap();
        match atom_ca.vdw_param {
            CachedVdwParam::Buckingham {
                radius, well_depth, ..
            } => {
                assert_eq!(radius, 2.0);
                assert_eq!(well_depth, 0.1);
            }
            _ => panic!("Expected Buckingham variant for CA atom"),
        }
        assert_eq!(atom_ca.hbond_type_id, -1);
    }

    #[test]
    fn parameterize_non_bonded_properties_fails_for_missing_vdw_params() {
        let ff = create_dummy_forcefield();
        let (mut system, atom_id) = create_test_system();
        let parameterizer = Parameterizer::new(ff, 0.0);

        system.atom_mut(atom_id).unwrap().force_field_type = "UNKNOWN_TYPE".to_string();

        let result = parameterizer.parameterize_non_bonded_properties(&mut system);
        assert!(matches!(
            result,
            Err(ParameterizationError::MissingVdwParams(ff_type)) if ff_type == "UNKNOWN_TYPE"
        ));
    }

    #[test]
    fn parameterize_atom_updates_atom_correctly() {
        let ff = create_dummy_forcefield_for_atom_test();
        let parameterizer = Parameterizer::new(ff, 1.0);
        let residue_id = crate::core::models::ids::ResidueId::default();

        let mut atom = Atom::new("CA", residue_id, Point3::origin());
        atom.partial_charge = -99.9;
        atom.force_field_type = "C_SP3".to_string();

        let result = parameterizer.parameterize_atom(&mut atom, "ALA");

        assert!(result.is_ok(), "Parameterization should succeed");
        assert_eq!(
            atom.partial_charge, -99.9,
            "Partial charge should NOT be modified by parameterize_atom"
        );
        assert!(
            (atom.delta - 1.43).abs() < 1e-9,
            "Delta value should be mu + s * sigma"
        );
        match atom.vdw_param {
            CachedVdwParam::Buckingham {
                radius, well_depth, ..
            } => {
                assert_eq!(radius, 2.0, "VDW radius should be set from forcefield");
                assert_eq!(
                    well_depth, 0.1,
                    "VDW well depth should be set from forcefield"
                );
            }
            _ => panic!("Expected Buckingham variant for CA atom"),
        }
        assert_eq!(
            atom.hbond_type_id, -1,
            "HBond type ID for non-HBond atom should be -1"
        );
    }

    #[test]
    fn parameterize_atom_handles_hbond_donor_hydrogen() {
        let residue_id = crate::core::models::ids::ResidueId::default();

        let mut atom = Atom::new("HN", residue_id, Point3::origin());
        atom.force_field_type = "H___A".to_string();
        atom.partial_charge = 0.35;

        let ff_for_hbond = create_dummy_forcefield_for_atom_test();
        let parameterizer_for_hbond = Parameterizer::new(ff_for_hbond, 1.0);
        let result = parameterizer_for_hbond.parameterize_atom(&mut atom, "ALA");

        assert!(result.is_ok());
        assert_eq!(atom.force_field_type, "H___A");
        assert_eq!(
            atom.hbond_type_id, 0,
            "Polar hydrogen should have hbond_type_id 0"
        );
    }
}
