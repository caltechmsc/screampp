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
use std::collections::HashMap;
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
        "Missing or misclassified anchor atom in residue '{residue_name}': Cannot find required anchor atom '{atom_name}', or it was incorrectly defined as a sidechain atom in the topology."
    )]
    InvalidAnchorAtom {
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
        // Step 1 & 2 remain the same: create pool and mark sidechains
        let mut atom_pool: HashMap<String, Vec<AtomId>> = HashMap::new();
        for &atom_id in atom_ids {
            let atom_name = system.atom(atom_id).unwrap().name.clone();
            atom_pool.entry(atom_name).or_default().push(atom_id);
        }

        for sidechain_name in &topology.sidechain_atoms {
            if let Some(ids) = atom_pool.get_mut(sidechain_name) {
                if let Some(atom_id_to_mark) = ids.pop() {
                    system.atom_mut(atom_id_to_mark).unwrap().role = AtomRole::Sidechain;
                }
            }
        }

        // Step 3 remains the same: mark remaining as backbone
        for ids in atom_pool.values() {
            for &atom_id in ids {
                system.atom_mut(atom_id).unwrap().role = AtomRole::Backbone;
            }
        }

        if atom_ids.is_empty() {
            if !topology.anchor_atoms.is_empty() {
                return Err(ParameterizationError::InvalidAnchorAtom {
                    residue_name: residue_name.to_string(),
                    atom_name: topology.anchor_atoms[0].clone(),
                });
            }
            return Ok(());
        }

        let parent_residue = system
            .residue(system.atom(atom_ids[0]).unwrap().residue_id)
            .unwrap();

        // Step 4: Validate that all required anchor atoms exist AND are classified as backbone.
        for anchor_name in &topology.anchor_atoms {
            match parent_residue.get_first_atom_id_by_name(anchor_name) {
                Some(atom_id) => {
                    let atom = system.atom(atom_id).unwrap();
                    if atom.role != AtomRole::Backbone {
                        return Err(ParameterizationError::InvalidAnchorAtom {
                            residue_name: residue_name.to_string(),
                            atom_name: anchor_name.clone(),
                        });
                    }
                }
                None => {
                    return Err(ParameterizationError::InvalidAnchorAtom {
                        residue_name: residue_name.to_string(),
                        atom_name: anchor_name.clone(),
                    });
                }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{
        forcefield::params::{DeltaParam, Forcefield, NonBondedParams, VdwParam},
        models::{chain::ChainType, residue::ResidueType, system::MolecularSystem},
        topology::registry::TopologyRegistry,
    };
    use nalgebra::Point3;
    use std::collections::HashMap;
    use std::io::Write;
    use tempfile::NamedTempFile;

    struct TestSetup {
        system: MolecularSystem,
        forcefield: Forcefield,
        topology_registry: TopologyRegistry,
    }

    fn setup() -> TestSetup {
        let mut system = MolecularSystem::new();

        let chain_a = system.add_chain('A', ChainType::Protein);
        let ala_id = system
            .add_residue(chain_a, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let unk_id = system.add_residue(chain_a, 2, "UNK", None).unwrap();

        let chain_l = system.add_chain('L', ChainType::Ligand);
        let lig_id = system.add_residue(chain_l, 101, "LIG", None).unwrap();

        let chain_w = system.add_chain('W', ChainType::Water);
        let hoh_id = system.add_residue(chain_w, 201, "HOH", None).unwrap();

        system
            .add_atom_to_residue(ala_id, Atom::new("N", ala_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(ala_id, Atom::new("CA", ala_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(ala_id, Atom::new("C", ala_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(ala_id, Atom::new("CB", ala_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(unk_id, Atom::new("X", unk_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(lig_id, Atom::new("C1", lig_id, Point3::origin()))
            .unwrap();
        system
            .add_atom_to_residue(hoh_id, Atom::new("O", hoh_id, Point3::origin()))
            .unwrap();

        let mut vdw = HashMap::new();
        vdw.insert(
            "C_BB".to_string(),
            VdwParam::LennardJones {
                radius: 1.8,
                well_depth: 0.1,
            },
        );
        vdw.insert(
            "C_SC".to_string(),
            VdwParam::LennardJones {
                radius: 1.9,
                well_depth: 0.12,
            },
        );
        vdw.insert(
            "O_H2O".to_string(),
            VdwParam::LennardJones {
                radius: 1.6,
                well_depth: 0.2,
            },
        );
        let mut deltas = HashMap::new();
        deltas.insert(
            ("ALA".to_string(), "CB".to_string()),
            DeltaParam {
                residue_type: "ALA".to_string(),
                atom_name: "CB".to_string(),
                mu: 0.1,
                sigma: 0.05,
            },
        );
        let forcefield = Forcefield {
            non_bonded: NonBondedParams {
                globals: Default::default(),
                vdw,
                hbond: HashMap::new(),
            },
            deltas,
        };

        let topo_content = r#"
[ALA]
anchor_atoms = ["N", "CA", "C"]
sidechain_atoms = ["CB"]
"#;
        let mut topo_file = NamedTempFile::new().unwrap();
        write!(topo_file, "{}", topo_content).unwrap();
        let topology_registry = TopologyRegistry::load(topo_file.path()).unwrap();

        TestSetup {
            system,
            forcefield,
            topology_registry,
        }
    }

    mod role_assignment {
        use super::*;

        #[test]
        fn assigns_roles_correctly_for_standard_protein_residue() {
            let TestSetup {
                mut system,
                forcefield,
                topology_registry,
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 1.0);

            parameterizer.parameterize_system(&mut system).unwrap();

            let ala_residue = system
                .residue(
                    system
                        .find_residue_by_id(system.find_chain_by_id('A').unwrap(), 1)
                        .unwrap(),
                )
                .unwrap();

            assert_eq!(
                system
                    .atom(ala_residue.get_first_atom_id_by_name("N").unwrap())
                    .unwrap()
                    .role,
                AtomRole::Backbone
            );
            assert_eq!(
                system
                    .atom(ala_residue.get_first_atom_id_by_name("CA").unwrap())
                    .unwrap()
                    .role,
                AtomRole::Backbone
            );
            assert_eq!(
                system
                    .atom(ala_residue.get_first_atom_id_by_name("C").unwrap())
                    .unwrap()
                    .role,
                AtomRole::Backbone
            );
            assert_eq!(
                system
                    .atom(ala_residue.get_first_atom_id_by_name("CB").unwrap())
                    .unwrap()
                    .role,
                AtomRole::Sidechain
            );
        }

        #[test]
        fn assigns_roles_correctly_for_ligand_and_water() {
            let TestSetup {
                mut system,
                forcefield,
                topology_registry,
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 1.0);

            parameterizer.parameterize_system(&mut system).unwrap();

            let lig_atom_id = system
                .residue(
                    system
                        .find_residue_by_id(system.find_chain_by_id('L').unwrap(), 101)
                        .unwrap(),
                )
                .unwrap()
                .get_first_atom_id_by_name("C1")
                .unwrap();
            assert_eq!(system.atom(lig_atom_id).unwrap().role, AtomRole::Ligand);

            let hoh_atom_id = system
                .residue(
                    system
                        .find_residue_by_id(system.find_chain_by_id('W').unwrap(), 201)
                        .unwrap(),
                )
                .unwrap()
                .get_first_atom_id_by_name("O")
                .unwrap();
            assert_eq!(system.atom(hoh_atom_id).unwrap().role, AtomRole::Water);
        }

        #[test]
        fn assigns_other_role_for_protein_residue_without_topology() {
            let TestSetup {
                mut system,
                forcefield,
                topology_registry,
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 1.0);

            parameterizer.parameterize_system(&mut system).unwrap();

            let unk_atom_id = system
                .residue(
                    system
                        .find_residue_by_id(system.find_chain_by_id('A').unwrap(), 2)
                        .unwrap(),
                )
                .unwrap()
                .get_first_atom_id_by_name("X")
                .unwrap();
            assert_eq!(system.atom(unk_atom_id).unwrap().role, AtomRole::Other);
        }
    }

    mod physicochemical_params {
        use super::*;

        #[test]
        fn assigns_physicochemical_params_correctly() {
            let TestSetup {
                mut system,
                forcefield,
                topology_registry,
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 1.0);

            let ala_cb_id = system
                .residue(
                    system
                        .find_residue_by_id(system.find_chain_by_id('A').unwrap(), 1)
                        .unwrap(),
                )
                .unwrap()
                .get_first_atom_id_by_name("CB")
                .unwrap();
            system.atom_mut(ala_cb_id).unwrap().force_field_type = "C_SC".to_string();

            parameterizer.parameterize_system(&mut system).unwrap();

            let cb_atom = system.atom(ala_cb_id).unwrap();

            assert!((cb_atom.delta - 0.15).abs() < 1e-9);

            match cb_atom.vdw_param {
                CachedVdwParam::LennardJones { radius, well_depth } => {
                    assert_eq!(radius, 1.9);
                    assert_eq!(well_depth, 0.12);
                }
                _ => panic!("Incorrect VDW param type cached"),
            }
            assert_eq!(cb_atom.hbond_type_id, -1);
        }
    }

    mod error_handling {
        use super::*;

        #[test]
        fn parameterizer_is_tolerant_of_missing_sidechain_atom_in_system() {
            let TestSetup {
                mut system,
                forcefield,
                topology_registry,
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 1.0);

            let ala_id = system
                .find_residue_by_id(system.find_chain_by_id('A').unwrap(), 1)
                .unwrap();

            let cb_id = system
                .residue(ala_id)
                .unwrap()
                .get_first_atom_id_by_name("CB")
                .unwrap();
            system.remove_atom(cb_id);

            let result = parameterizer.parameterize_system(&mut system);
            assert!(
                result.is_ok(),
                "Parameterization should succeed even with a missing sidechain atom, but it failed with: {:?}",
                result.err()
            );

            let n_id = system
                .residue(ala_id)
                .unwrap()
                .get_first_atom_id_by_name("N")
                .unwrap();
            assert_eq!(system.atom(n_id).unwrap().role, AtomRole::Backbone);
        }

        #[test]
        fn fails_when_anchor_atom_is_misclassified_as_sidechain_in_topology() {
            let TestSetup {
                mut system,
                forcefield,
                ..
            } = setup();

            let faulty_topo_content = r#"
[ALA]
anchor_atoms = ["N", "CA", "C"]
sidechain_atoms = ["C", "CB"]
"#;
            let mut faulty_topo_file = NamedTempFile::new().unwrap();
            write!(faulty_topo_file, "{}", faulty_topo_content).unwrap();
            let faulty_topology_registry = TopologyRegistry::load(faulty_topo_file.path()).unwrap();

            let parameterizer = Parameterizer::new(&forcefield, &faulty_topology_registry, 1.0);

            let result = parameterizer.parameterize_system(&mut system);

            assert_eq!(
                result.unwrap_err(),
                ParameterizationError::InvalidAnchorAtom {
                    residue_name: "ALA".to_string(),
                    atom_name: "C".to_string()
                }
            );
        }

        #[test]
        fn fails_on_missing_vdw_params() {
            let TestSetup {
                mut system,
                forcefield,
                topology_registry,
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 1.0);

            let lig_id = system
                .find_residue_by_id(system.find_chain_by_id('L').unwrap(), 101)
                .unwrap();
            let lig_atom_id = system
                .residue(lig_id)
                .unwrap()
                .get_first_atom_id_by_name("C1")
                .unwrap();
            system.atom_mut(lig_atom_id).unwrap().force_field_type = "UNKNOWN_TYPE".to_string();

            let result = parameterizer.parameterize_system(&mut system);
            assert!(
                matches!(result, Err(ParameterizationError::MissingVdwParams { ff_type, .. }) if ff_type == "UNKNOWN_TYPE")
            );
        }
    }

    mod single_atom_parameterization {
        use super::*;

        #[test]
        fn parameterize_protein_atom_works_correctly() {
            let TestSetup {
                forcefield,
                topology_registry,
                ..
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 1.0);

            let ala_topology = topology_registry.get("ALA").unwrap();

            let mut cb_atom = Atom::new("CB", ResidueId::default(), Point3::origin());
            cb_atom.force_field_type = "C_SC".to_string();
            parameterizer
                .parameterize_protein_atom(&mut cb_atom, "ALA", ala_topology)
                .unwrap();

            assert_eq!(cb_atom.role, AtomRole::Sidechain);
            assert!((cb_atom.delta - 0.15).abs() < 1e-9);

            let mut ca_atom = Atom::new("CA", ResidueId::default(), Point3::origin());
            ca_atom.force_field_type = "C_BB".to_string();
            parameterizer
                .parameterize_protein_atom(&mut ca_atom, "ALA", ala_topology)
                .unwrap();

            assert_eq!(ca_atom.role, AtomRole::Backbone);
            assert_eq!(ca_atom.delta, 0.0);
        }
    }
}
