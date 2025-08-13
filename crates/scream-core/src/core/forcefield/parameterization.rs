use super::params::Forcefield;
use crate::core::{
    models::{
        atom::{Atom, AtomRole, CachedVdwParam},
        chain::ChainType,
        ids::{AtomId, ResidueId},
        system::MolecularSystem,
    },
    rotamers::rotamer::Rotamer,
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

struct CalculatedAtomParams {
    delta: f64,
    vdw_param: CachedVdwParam,
    hbond_type_id: i8,
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
        // Pass 1: Assign Atom Roles for each residue.
        for residue_id in system.residues_iter().map(|(id, _)| id).collect::<Vec<_>>() {
            self.assign_atom_roles_for_residue(residue_id, system)?;
        }

        // Pass 2, Step A: "Collect" phase - compute all parameters, but do not modify the system yet
        let mut calculated_params = HashMap::new();
        for (atom_id, atom) in system.atoms_iter() {
            let residue = system.residue(atom.residue_id).unwrap();

            let (delta, vdw_param) = self.calculate_core_params(atom, &residue.name)?;

            let hbond_type_id = self.determine_hbond_role_for_system_atom(atom_id, system);

            calculated_params.insert(
                atom_id,
                CalculatedAtomParams {
                    delta,
                    vdw_param,
                    hbond_type_id,
                },
            );
        }

        // Pass 2, Step B: "Apply" phase - now it's safe to mutably borrow the system
        for (atom_id, params) in calculated_params {
            let atom = system.atom_mut(atom_id).unwrap();
            atom.delta = params.delta;
            atom.vdw_param = params.vdw_param;
            atom.hbond_type_id = params.hbond_type_id;
        }

        Ok(())
    }

    pub fn parameterize_rotamer(
        &self,
        rotamer: &mut Rotamer,
        residue_name: &str,
        topology: &ResidueTopology,
    ) -> Result<(), ParameterizationError> {
        assign_protein_roles_from_pool(
            &mut rotamer.atoms,
            |atom| &atom.name,
            |atom, role| atom.role = role,
            topology,
            residue_name,
        )?;

        for i in 0..rotamer.atoms.len() {
            let (delta, vdw_param) = {
                let atom = &rotamer.atoms[i];
                self.calculate_core_params(atom, residue_name)?
            };

            let hbond_type_id = self.determine_hbond_role_for_rotamer_atom(i, rotamer);

            let atom = &mut rotamer.atoms[i];
            atom.delta = delta;
            atom.vdw_param = vdw_param;
            atom.hbond_type_id = hbond_type_id;
        }

        Ok(())
    }

    fn assign_atom_roles_for_residue(
        &self,
        residue_id: ResidueId,
        system: &mut MolecularSystem,
    ) -> Result<(), ParameterizationError> {
        let (chain_type, residue_name) = {
            let residue = system.residue(residue_id).unwrap();
            let chain = system.chain(residue.chain_id).unwrap();
            (chain.chain_type, residue.name.clone())
        };

        let atom_ids: Vec<AtomId> = system.residue(residue_id).unwrap().atoms().to_vec();

        match chain_type {
            ChainType::Protein | ChainType::DNA | ChainType::RNA => {
                if let Some(topology) = self.topology_registry.get(&residue_name) {
                    let mut atoms_view: Vec<&mut Atom> = system
                        .atoms_iter_mut()
                        .filter(|(id, _)| atom_ids.contains(id))
                        .map(|(_, atom)| atom)
                        .collect();

                    assign_protein_roles_from_pool(
                        &mut atoms_view,
                        |atom| &atom.name,
                        |atom, role| atom.role = role,
                        topology,
                        &residue_name,
                    )?;
                } else {
                    warn!(
                        "Residue '{}' has no topology definition. Marking all its atoms as 'Other'.",
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

    fn calculate_core_params(
        &self,
        atom: &Atom,
        residue_name: &str,
    ) -> Result<(f64, CachedVdwParam), ParameterizationError> {
        let delta = self
            .forcefield
            .deltas
            .get(&(residue_name.to_string(), atom.name.clone()))
            .map_or(0.0, |p| p.mu + self.delta_s_factor * p.sigma);

        let vdw_param = if !atom.force_field_type.is_empty() {
            self.forcefield
                .non_bonded
                .vdw
                .get(&atom.force_field_type)
                .map(|p| p.clone().into())
                .unwrap_or(CachedVdwParam::None)
        } else {
            CachedVdwParam::None
        };

        Ok((delta, vdw_param))
    }

    fn determine_hbond_role_for_system_atom(
        &self,
        atom_id: AtomId,
        system: &MolecularSystem,
    ) -> i8 {
        let atom = system.atom(atom_id).unwrap();
        self.determine_hbond_role(&atom.name, &atom.force_field_type, || {
            system.get_bonded_neighbors(atom_id).and_then(|neighbors| {
                if neighbors.len() == 1 {
                    system
                        .atom(neighbors[0])
                        .map(|heavy_atom| heavy_atom.force_field_type.as_str())
                } else {
                    None
                }
            })
        })
    }

    fn determine_hbond_role_for_rotamer_atom(&self, atom_index: usize, rotamer: &Rotamer) -> i8 {
        let atom = &rotamer.atoms[atom_index];
        self.determine_hbond_role(&atom.name, &atom.force_field_type, || {
            let neighbors: Vec<usize> = rotamer
                .bonds
                .iter()
                .filter_map(|&(i, j)| {
                    if i == atom_index {
                        Some(j)
                    } else if j == atom_index {
                        Some(i)
                    } else {
                        None
                    }
                })
                .collect();

            if neighbors.len() == 1 {
                Some(rotamer.atoms[neighbors[0]].force_field_type.as_str())
            } else {
                None
            }
        })
    }

    fn determine_hbond_role<'b, F>(
        &self,
        atom_name: &str,
        ff_type: &str,
        get_neighbor_ff_type: F,
    ) -> i8
    where
        F: Fn() -> Option<&'b str>,
    {
        let mut hbond_type_id = -1;

        if self.forcefield.non_bonded.hbond_acceptors.contains(ff_type) {
            hbond_type_id = 1;
        }

        if atom_name.starts_with('H') {
            if let Some(heavy_atom_ff_type) = get_neighbor_ff_type() {
                if self
                    .forcefield
                    .non_bonded
                    .hbond_donors
                    .contains(heavy_atom_ff_type)
                {
                    hbond_type_id = 0;
                }
            }
        }

        hbond_type_id
    }
}

fn assign_protein_roles_from_pool<'a, T>(
    atoms: &'a mut [T],
    get_name: impl for<'b> Fn(&'b T) -> &'b str,
    set_role: impl Fn(&mut T, AtomRole),
    topology: &ResidueTopology,
    residue_name: &str,
) -> Result<(), ParameterizationError> {
    // Step 1: Create a pool mapping atom names to their indices in the slice.
    let mut atom_pool: HashMap<String, Vec<usize>> = HashMap::new();
    for (index, atom) in atoms.iter().enumerate() {
        atom_pool
            .entry(get_name(atom).to_string())
            .or_default()
            .push(index);
    }

    let mut assigned_indices = HashSet::new();

    // Step 2: First pass - identify and reserve anchor atoms. Consume them from the FRONT of the pool.
    for anchor_name in &topology.anchor_atoms {
        let indices = atom_pool.get_mut(anchor_name).ok_or_else(|| {
            ParameterizationError::InvalidAnchorAtom {
                residue_name: residue_name.to_string(),
                atom_name: anchor_name.clone(),
            }
        })?;
        if indices.is_empty() {
            return Err(ParameterizationError::InvalidAnchorAtom {
                residue_name: residue_name.to_string(),
                atom_name: anchor_name.clone(),
            });
        }
        let index_to_mark = indices.remove(0);
        set_role(&mut atoms[index_to_mark], AtomRole::Backbone);
        assigned_indices.insert(index_to_mark);
    }

    // Step 3: Second pass - identify and assign sidechain atoms. Consume them from the BACK of the pool.
    for sidechain_name in &topology.sidechain_atoms {
        if let Some(indices) = atom_pool.get_mut(sidechain_name) {
            if let Some(index_to_mark) = indices.pop() {
                if assigned_indices.contains(&index_to_mark) {
                    return Err(ParameterizationError::InvalidAnchorAtom {
                        residue_name: residue_name.to_string(),
                        atom_name: sidechain_name.clone(),
                    });
                }
                set_role(&mut atoms[index_to_mark], AtomRole::Sidechain);
                assigned_indices.insert(index_to_mark);
            }
        }
    }

    // Step 4: Final pass - any atom not yet assigned a role is part of the backbone.
    for (index, atom) in atoms.iter_mut().enumerate() {
        if !assigned_indices.contains(&index) {
            set_role(atom, AtomRole::Backbone);
        }
    }

    Ok(())
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
        models::{chain::ChainType, residue::ResidueType, system::MolecularSystem},
        rotamers::rotamer::Rotamer,
        topology::registry::TopologyRegistry,
    };
    use nalgebra::Point3;
    use std::{fs::File, io::Write};
    use tempfile::{TempDir, tempdir};

    struct TestSetup {
        forcefield: Forcefield,
        topology_registry: TopologyRegistry,
        _temp_dir: TempDir,
    }

    fn setup() -> TestSetup {
        let temp_dir = tempdir().unwrap();
        let dir = temp_dir.path();

        let ff_path = dir.join("ff.toml");
        let mut ff_file = File::create(&ff_path).unwrap();
        write!(
            ff_file,
            r#"
            [globals]
            dielectric_constant = 4.0
            potential_function = "lennard-jones-12-6"
            [vdw]
            N_R = {{ radius = 1.6, well_depth = 0.1 }}
            C_BB = {{ radius = 1.8, well_depth = 0.1 }}
            C_R = {{ radius = 1.8, well_depth = 0.1 }}
            C_SC = {{ radius = 1.9, well_depth = 0.12 }}
            H_ = {{ radius = 1.0, well_depth = 0.02 }}
            O_2 = {{ radius = 1.5, well_depth = 0.2 }}
            H_O = {{ radius = 1.0, well_depth = 0.01 }}
            [hbond]
            H_O = {{ equilibrium_distance = 2.8, well_depth = 5.0 }}
            O_2 = {{ equilibrium_distance = 2.8, well_depth = 5.0 }}
        "#
        )
        .unwrap();

        let delta_path = dir.join("delta.csv");
        let mut delta_file = File::create(&delta_path).unwrap();
        write!(
            delta_file,
            "residue_type,atom_name,mu,sigma\nALA,CB,0.1,0.05\n"
        )
        .unwrap();

        let forcefield = Forcefield::load(&ff_path, &delta_path).unwrap();

        let topo_path = dir.join("topo.toml");
        let mut topo_file = File::create(&topo_path).unwrap();
        write!(
            topo_file,
            r#"
            [ALA]
            anchor_atoms = ["N", "CA", "C"]
            sidechain_atoms = ["CB"]
            [GLY]
            anchor_atoms = ["N", "CA", "C", "HA2"]
            sidechain_atoms = ["HA1"]
        "#
        )
        .unwrap();
        let topology_registry = TopologyRegistry::load(&topo_path).unwrap();

        TestSetup {
            forcefield,
            topology_registry,
            _temp_dir: temp_dir,
        }
    }

    mod role_assignment {
        use super::*;

        #[test]
        fn assigns_roles_correctly_for_standard_protein_residue() {
            let TestSetup {
                forcefield,
                topology_registry,
                ..
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 1.0);

            let mut system = MolecularSystem::new();
            let chain_id = system.add_chain('A', ChainType::Protein);
            let ala_id = system
                .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
                .unwrap();
            for name in ["N", "CA", "C", "CB"].iter() {
                system
                    .add_atom_to_residue(ala_id, Atom::new(name, ala_id, Point3::origin()))
                    .unwrap();
            }

            parameterizer.parameterize_system(&mut system).unwrap();

            let ala_res = system.residue(ala_id).unwrap();
            assert_eq!(
                system
                    .atom(ala_res.get_first_atom_id_by_name("N").unwrap())
                    .unwrap()
                    .role,
                AtomRole::Backbone
            );
            assert_eq!(
                system
                    .atom(ala_res.get_first_atom_id_by_name("CA").unwrap())
                    .unwrap()
                    .role,
                AtomRole::Backbone
            );
            assert_eq!(
                system
                    .atom(ala_res.get_first_atom_id_by_name("C").unwrap())
                    .unwrap()
                    .role,
                AtomRole::Backbone
            );
            assert_eq!(
                system
                    .atom(ala_res.get_first_atom_id_by_name("CB").unwrap())
                    .unwrap()
                    .role,
                AtomRole::Sidechain
            );
        }

        #[test]
        fn assigns_roles_for_non_protein_chains() {
            let TestSetup {
                forcefield,
                topology_registry,
                ..
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 1.0);

            let mut system = MolecularSystem::new();
            let lig_chain = system.add_chain('L', ChainType::Ligand);
            let wat_chain = system.add_chain('W', ChainType::Water);
            let lig_id = system.add_residue(lig_chain, 1, "LIG", None).unwrap();
            let wat_id = system.add_residue(wat_chain, 2, "HOH", None).unwrap();
            let lig_atom_id = system
                .add_atom_to_residue(lig_id, Atom::new("C1", lig_id, Point3::origin()))
                .unwrap();
            let wat_atom_id = system
                .add_atom_to_residue(wat_id, Atom::new("O", wat_id, Point3::origin()))
                .unwrap();

            parameterizer.parameterize_system(&mut system).unwrap();

            assert_eq!(system.atom(lig_atom_id).unwrap().role, AtomRole::Ligand);
            assert_eq!(system.atom(wat_atom_id).unwrap().role, AtomRole::Water);
        }

        #[test]
        fn is_tolerant_of_missing_sidechain_atoms_in_system() {
            let TestSetup {
                forcefield,
                topology_registry,
                ..
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 1.0);

            let mut system = MolecularSystem::new();
            let chain_id = system.add_chain('A', ChainType::Protein);
            let ala_id = system
                .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
                .unwrap();

            for name in ["N", "CA", "C"].iter() {
                system
                    .add_atom_to_residue(ala_id, Atom::new(name, ala_id, Point3::origin()))
                    .unwrap();
            }

            let result = parameterizer.parameterize_system(&mut system);
            assert!(
                result.is_ok(),
                "Parameterization should succeed with missing sidechain atom, but failed: {:?}",
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
        fn fails_if_anchor_atom_is_missing_from_system() {
            let TestSetup {
                forcefield,
                topology_registry,
                ..
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 1.0);

            let mut system = MolecularSystem::new();
            let chain_id = system.add_chain('A', ChainType::Protein);
            let ala_id = system
                .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
                .unwrap();

            for name in ["N", "C", "CB"].iter() {
                system
                    .add_atom_to_residue(ala_id, Atom::new(name, ala_id, Point3::origin()))
                    .unwrap();
            }

            let result = parameterizer.parameterize_system(&mut system);
            assert_eq!(
                result.unwrap_err(),
                ParameterizationError::InvalidAnchorAtom {
                    residue_name: "ALA".to_string(),
                    atom_name: "CA".to_string()
                }
            );
        }
    }

    mod physicochemical_params_assignment {
        use super::*;

        #[test]
        fn assign_physicochemical_params_sets_all_fields_correctly() {
            let TestSetup {
                forcefield,
                topology_registry,
                ..
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 2.0);

            let mut cb_atom = Atom::new("CB", ResidueId::default(), Default::default());
            cb_atom.force_field_type = "C_SC".to_string();

            parameterizer
                .assign_physicochemical_params(&mut cb_atom, "ALA")
                .unwrap();

            assert!(
                (cb_atom.delta - 0.2).abs() < 1e-9,
                "Delta value for ALA CB is incorrect"
            );
            match cb_atom.vdw_param {
                CachedVdwParam::LennardJones { radius, well_depth } => {
                    assert_eq!(radius, 1.9);
                    assert_eq!(well_depth, 0.12);
                }
                _ => panic!("Incorrect VDW param type for CB"),
            }
            assert_eq!(cb_atom.hbond_type_id, -1, "CB should not be an H-bond atom");

            let mut ca_atom = Atom::new("CA", ResidueId::default(), Default::default());
            ca_atom.force_field_type = "C_BB".to_string();

            parameterizer
                .assign_physicochemical_params(&mut ca_atom, "ALA")
                .unwrap();

            assert_eq!(
                ca_atom.delta, 0.0,
                "Delta for ALA CA should be 0 as it's not in the delta file"
            );
            match ca_atom.vdw_param {
                CachedVdwParam::LennardJones { radius, well_depth } => {
                    assert_eq!(radius, 1.8);
                    assert_eq!(well_depth, 0.1);
                }
                _ => panic!("Incorrect VDW param type for CA"),
            }

            let mut h_atom = Atom::new("HG", ResidueId::default(), Default::default());
            h_atom.force_field_type = "H_O".to_string();

            parameterizer
                .assign_physicochemical_params(&mut h_atom, "SER")
                .unwrap();
            assert_eq!(
                h_atom.hbond_type_id, 0,
                "H-bond donor hydrogen should have type_id 0"
            );

            let mut o_atom = Atom::new("O", ResidueId::default(), Default::default());
            o_atom.force_field_type = "O_2".to_string();

            parameterizer
                .assign_physicochemical_params(&mut o_atom, "ALA")
                .unwrap();
            assert_eq!(
                o_atom.hbond_type_id, 1,
                "H-bond acceptor should have type_id 1"
            );
        }
    }

    mod parameterize_rotamer_tests {
        use super::*;

        #[test]
        fn parameterize_rotamer_assigns_all_properties_correctly() {
            let TestSetup {
                forcefield,
                topology_registry,
                ..
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 1.0);
            let ala_topology = topology_registry.get("ALA").unwrap();

            let mut ala_rotamer = Rotamer {
                atoms: vec![
                    Atom::new("N", ResidueId::default(), Default::default()),
                    Atom::new("CA", ResidueId::default(), Default::default()),
                    Atom::new("C", ResidueId::default(), Default::default()),
                    Atom::new("CB", ResidueId::default(), Default::default()),
                ],
                bonds: vec![],
            };
            ala_rotamer.atoms.iter_mut().for_each(|a| {
                a.force_field_type = match a.name.as_str() {
                    "N" => "N_R".to_string(),
                    "CA" => "C_BB".to_string(),
                    "C" => "C_R".to_string(),
                    "CB" => "C_SC".to_string(),
                    _ => panic!("Unknown atom"),
                }
            });

            parameterizer
                .parameterize_rotamer(&mut ala_rotamer, "ALA", ala_topology)
                .unwrap();

            let get_atom = |name: &str| ala_rotamer.atoms.iter().find(|a| a.name == name).unwrap();

            let cb_atom = get_atom("CB");
            assert_eq!(cb_atom.role, AtomRole::Sidechain);
            assert!((cb_atom.delta - (0.1 + 1.0 * 0.05)).abs() < 1e-9);
            assert!(
                matches!(cb_atom.vdw_param, CachedVdwParam::LennardJones { radius, .. } if radius == 1.9)
            );
            assert_eq!(cb_atom.hbond_type_id, -1);

            let ca_atom = get_atom("CA");
            assert_eq!(ca_atom.role, AtomRole::Backbone);
            assert_eq!(ca_atom.delta, 0.0);
            assert!(
                matches!(ca_atom.vdw_param, CachedVdwParam::LennardJones { radius, .. } if radius == 1.8)
            );
        }
    }
}
