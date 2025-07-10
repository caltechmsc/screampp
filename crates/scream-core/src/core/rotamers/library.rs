use super::placement::PlacementInfo;
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
    pub placement_info: HashMap<ResidueType, PlacementInfo>,
}

#[derive(Debug, Error)]
pub enum LibraryLoadError {
    #[error("File I/O error for '{path}': {source}")]
    Io {
        path: String,
        source: std::io::Error,
    },
    #[error("TOML parsing error for '{path}': {source}")]
    Toml {
        path: String,
        source: toml::de::Error,
    },
    #[error("Unknown residue type '{0}' found in library file")]
    UnknownResidueType(String),
    #[error(
        "Missing placement info for residue type '{0}' which is present in the rotamer library"
    )]
    MissingPlacementInfo(String),
    #[error(
        "Parameterization failed for residue '{res_type}' in rotamer from file '{path}': {source}"
    )]
    Parameterization {
        path: String,
        res_type: String,
        source: ParameterizationError,
    },
}

impl RotamerLibrary {
    pub fn load(
        rotamer_toml_path: &Path,
        placement_registry: &HashMap<String, PlacementInfo>,
        forcefield: &Forcefield,
        delta_s_factor: f64,
    ) -> Result<Self, LibraryLoadError> {
        let content =
            std::fs::read_to_string(rotamer_toml_path).map_err(|e| LibraryLoadError::Io {
                path: rotamer_toml_path.to_string_lossy().to_string(),
                source: e,
            })?;
        let raw_lib: RawRotamerFile =
            toml::from_str(&content).map_err(|e| LibraryLoadError::Toml {
                path: rotamer_toml_path.to_string_lossy().to_string(),
                source: e,
            })?;

        let parameterizer = Parameterizer::new(forcefield.clone(), delta_s_factor);
        let mut final_rotamers_map = HashMap::new();
        let mut final_placement_map = HashMap::new();

        for (res_name, raw_rotamer_list) in raw_lib {
            let res_type = ResidueType::from_str(&res_name)
                .map_err(|_| LibraryLoadError::UnknownResidueType(res_name.clone()))?;

            let placement_data = placement_registry
                .get(&res_name)
                .ok_or_else(|| LibraryLoadError::MissingPlacementInfo(res_name.clone()))?
                .clone();
            final_placement_map.insert(res_type, placement_data);

            let mut parameterized_rotamers = Vec::with_capacity(raw_rotamer_list.len());
            let placeholder_residue_id = ResidueId::default();

            for raw_rotamer in raw_rotamer_list {
                let mut atoms = Vec::with_capacity(raw_rotamer.atoms.len());
                for atom_data in raw_rotamer.atoms {
                    let mut atom = Atom::new(
                        atom_data.serial,
                        &atom_data.atom_name,
                        placeholder_residue_id,
                        Point3::from(atom_data.position),
                    );
                    atom.partial_charge = atom_data.partial_charge;

                    parameterizer
                        .parameterize_atom(&mut atom, &res_name)
                        .map_err(|e| LibraryLoadError::Parameterization {
                            path: rotamer_toml_path.to_string_lossy().to_string(),
                            res_type: res_name.clone(),
                            source: e,
                        })?;

                    atoms.push(atom);
                }
                parameterized_rotamers.push(Rotamer { atoms });
            }

            final_rotamers_map.insert(res_type, parameterized_rotamers);
        }

        Ok(Self {
            rotamers: final_rotamers_map,
            placement_info: final_placement_map,
        })
    }

    pub fn get_rotamers_for(&self, res_type: ResidueType) -> Option<&Vec<Rotamer>> {
        self.rotamers.get(&res_type)
    }

    pub fn get_placement_info_for(&self, res_type: ResidueType) -> Option<&PlacementInfo> {
        self.placement_info.get(&res_type)
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
    use crate::core::models::ids::ResidueId;
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

    fn create_dummy_forcefield_for_atom_test() -> Forcefield {
        let mut ff = create_dummy_forcefield();

        ff.deltas.insert(
            ("ALA".to_string(), "CA".to_string()),
            DeltaParam {
                res_type: "ALA".to_string(),
                atom_name: "CA".to_string(),
                mu: 1.23,
                sigma: 0.2,
            },
        );
        ff
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

        let parameterizer = Parameterizer::new(ff, 0.0);
        let result = parameterizer.parameterize_topology(&mut system);
        assert!(matches!(
            result,
            Err(ParameterizationError::AtomNotFoundInTopology { res_type, atom_name })
            if res_type == "ALA" && atom_name == "H_EXTRA"
        ));
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

    #[test]
    fn parameterize_atom_updates_atom_correctly() {
        let ff = create_dummy_forcefield_for_atom_test();
        let parameterizer = Parameterizer::new(ff, 1.0);
        let residue_id = ResidueId::default();

        let mut atom = Atom::new(1, "CA", residue_id, Point3::origin());
        atom.partial_charge = -99.9;
        atom.force_field_type = "C_from_lib".to_string();

        let result = parameterizer.parameterize_atom(&mut atom, "ALA");

        assert!(result.is_ok(), "Parameterization should succeed");

        assert_eq!(
            atom.force_field_type, "C_SP3",
            "Force field type should be overwritten by topology"
        );

        assert_eq!(
            atom.partial_charge, -99.9,
            "Partial charge should NOT be modified by parameterize_atom"
        );

        assert!(
            (atom.delta - 1.43).abs() < 1e-9,
            "Delta value should be mu + s * sigma"
        );

        assert_eq!(
            atom.vdw_radius, 2.0,
            "VDW radius should be set from forcefield"
        );
        assert_eq!(
            atom.vdw_well_depth, 0.1,
            "VDW well depth should be set from forcefield"
        );

        assert_eq!(
            atom.hbond_type_id, -1,
            "HBond type ID for non-HBond atom should be -1"
        );
    }

    #[test]
    fn parameterize_atom_handles_hbond_donor_hydrogen() {
        let residue_id = ResidueId::default();

        let mut atom = Atom::new(1, "HN", residue_id, Point3::origin());
        atom.force_field_type = "H___A".to_string();
        atom.partial_charge = 0.35;

        let mut ff_for_hbond = create_dummy_forcefield_for_atom_test();
        ff_for_hbond
            .topology
            .get_mut("ALA")
            .unwrap()
            .atoms
            .push(TopologyAtomParam {
                name: "HN".to_string(),
                ff_type: "H___A".to_string(),
            });
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
