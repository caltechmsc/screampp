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

/// Defines the force field type identifier for hydrogen atoms that can participate as donors in hydrogen bonding.
///
/// This constant represents the specific force field type used to identify hydrogen atoms
/// that are capable of forming hydrogen bonds as donors in the Dreiding force field.
const DREIDING_HBOND_DONOR_HYDROGEN: &str = "H___A";

/// Represents errors that can occur during the parameterization process.
///
/// This enum encapsulates various error conditions that may arise when assigning
/// force field parameters to molecular systems or rotamers, providing detailed
/// context for debugging and error handling in molecular simulations.
#[derive(Debug, Error, PartialEq, Eq)]
pub enum ParameterizationError {
    /// Indicates that a required van der Waals parameter is missing for a specific force field type.
    ///
    /// This error occurs when an atom's force field type is defined but no corresponding
    /// VDW parameters (Lennard-Jones or Buckingham) are found in the force field data.
    #[error(
        "Missing VDW parameter for force field type: '{ff_type}' in atom '{atom_name}' of residue {residue_name}"
    )]
    MissingVdwParams {
        /// The force field type that is missing parameters.
        ff_type: String,
        /// The name of the atom with missing parameters.
        atom_name: String,
        /// The name of the residue containing the atom.
        residue_name: String,
    },
    /// Indicates that a required anchor atom is missing or misclassified in the residue topology.
    ///
    /// This error occurs when the topology definition specifies an anchor atom that is either
    /// not present in the system or has been incorrectly classified as a sidechain atom.
    #[error(
        "Missing or misclassified anchor atom in residue '{residue_name}': Cannot find required anchor atom '{atom_name}', or it was incorrectly defined as a sidechain atom in the topology."
    )]
    InvalidAnchorAtom {
        /// The name of the residue with the invalid anchor atom.
        residue_name: String,
        /// The name of the missing or misclassified anchor atom.
        atom_name: String,
    },
}

/// Handles the assignment of force field parameters to molecular systems and rotamers.
///
/// This struct provides methods to parameterize atoms with their physicochemical properties,
/// including van der Waals parameters, delta values for flat-bottom potentials, and hydrogen
/// bonding classifications. It uses topology information to assign atom roles and force field
/// data to compute interaction parameters.
pub struct Parameterizer<'a> {
    /// Reference to the force field containing parameter definitions.
    forcefield: &'a Forcefield,
    /// Reference to the topology registry for residue definitions.
    topology_registry: &'a TopologyRegistry,
    /// Scaling factor for the sigma component of delta parameters.
    delta_s_factor: f64,
}

/// Internal structure for storing calculated atom parameters before applying them.
///
/// This struct holds the computed parameters for an atom during the parameterization
/// process, allowing for a two-phase approach where parameters are calculated first
/// and then applied to avoid borrowing conflicts.
struct CalculatedAtomParams {
    /// The delta value for flat-bottom potential modifications.
    delta: f64,
    /// The cached van der Waals parameter for the atom.
    vdw_param: CachedVdwParam,
    /// The hydrogen bonding classification ID (-1: none, 0: donor, 1: acceptor).
    hbond_type_id: i8,
}

impl<'a> Parameterizer<'a> {
    /// Creates a new parameterizer with the given force field and topology registry.
    ///
    /// # Arguments
    ///
    /// * `forcefield` - The force field containing parameter definitions
    /// * `topology_registry` - The registry of residue topologies
    /// * `delta_s_factor` - Scaling factor for delta sigma values
    ///
    /// # Return
    ///
    /// Returns a new `Parameterizer` instance configured with the provided references.
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

    /// Parameterizes all atoms in a molecular system with force field properties.
    ///
    /// This method performs a complete parameterization of the system in two passes:
    /// first assigning atom roles based on topology, then computing physicochemical
    /// parameters for each atom. The two-pass approach ensures that mutable borrowing
    /// conflicts are avoided during parameter calculation.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system to parameterize (modified in place)
    ///
    /// # Return
    ///
    /// Returns `Ok(())` if parameterization succeeds, or an error if required parameters
    /// or topology information is missing.
    ///
    /// # Errors
    ///
    /// Returns `ParameterizationError::InvalidAnchorAtom` if required anchor atoms
    /// are missing from any residue.
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

    /// Parameterizes a rotamer with force field properties.
    ///
    /// This method assigns atom roles based on the provided topology and computes
    /// physicochemical parameters for each atom in the rotamer. It is designed for
    /// use with individual rotamer conformations during side-chain placement.
    ///
    /// # Arguments
    ///
    /// * `rotamer` - The rotamer to parameterize (modified in place)
    /// * `residue_name` - The name of the residue this rotamer represents
    /// * `topology` - The topology definition for the residue
    ///
    /// # Return
    ///
    /// Returns `Ok(())` if parameterization succeeds, or an error if required parameters
    /// or topology information is missing.
    ///
    /// # Errors
    ///
    /// Returns `ParameterizationError::InvalidAnchorAtom` if required anchor atoms
    /// are missing from the rotamer.
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

    /// Assigns atom roles for a specific residue in the molecular system.
    ///
    /// This method determines the role (backbone, sidechain, etc.) of each atom in a residue
    /// based on the chain type and available topology information. For protein residues,
    /// it uses topology definitions to classify atoms; for other chain types, it assigns
    /// roles based on the chain category.
    ///
    /// # Arguments
    ///
    /// * `residue_id` - The ID of the residue to parameterize
    /// * `system` - The molecular system containing the residue
    ///
    /// # Return
    ///
    /// Returns `Ok(())` if role assignment succeeds, or an error if required topology
    /// information is missing.
    ///
    /// # Errors
    ///
    /// Returns `ParameterizationError::InvalidAnchorAtom` if required anchor atoms
    /// are missing from the residue.
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

    /// Calculates the core physicochemical parameters for an atom.
    ///
    /// This method computes the delta value (for flat-bottom potentials) and van der Waals
    /// parameters for a single atom based on its force field type and residue-specific
    /// delta parameters.
    ///
    /// # Arguments
    ///
    /// * `atom` - The atom to parameterize
    /// * `residue_name` - The name of the residue containing the atom
    ///
    /// # Return
    ///
    /// Returns a tuple of `(delta, vdw_param)` if calculation succeeds.
    ///
    /// # Errors
    ///
    /// This method does not return errors directly, but the VDW parameter lookup
    /// may result in `CachedVdwParam::None` if the force field type is not found.
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

    /// Determines the hydrogen bonding role for an atom in the molecular system.
    ///
    /// This method analyzes the atom's force field type and its bonding environment
    /// to classify it as a hydrogen bond donor, acceptor, or neither.
    ///
    /// # Arguments
    ///
    /// * `atom_id` - The ID of the atom to classify
    /// * `system` - The molecular system containing the atom
    ///
    /// # Return
    ///
    /// Returns an integer indicating the hydrogen bonding role:
    /// - `-1`: Not involved in hydrogen bonding
    /// - `0`: Hydrogen bond donor
    /// - `1`: Hydrogen bond acceptor
    fn determine_hbond_role_for_system_atom(
        &self,
        atom_id: AtomId,
        system: &MolecularSystem,
    ) -> i8 {
        let atom = system.atom(atom_id).unwrap();
        self.determine_hbond_role(&atom.force_field_type, || {
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

    /// Determines the hydrogen bonding role for an atom in a rotamer.
    ///
    /// This method analyzes the atom's force field type and its bonding environment
    /// within the rotamer to classify it as a hydrogen bond donor, acceptor, or neither.
    ///
    /// # Arguments
    ///
    /// * `atom_index` - The index of the atom in the rotamer's atom list
    /// * `rotamer` - The rotamer containing the atom
    ///
    /// # Return
    ///
    /// Returns an integer indicating the hydrogen bonding role:
    /// - `-1`: Not involved in hydrogen bonding
    /// - `0`: Hydrogen bond donor
    /// - `1`: Hydrogen bond acceptor
    fn determine_hbond_role_for_rotamer_atom(&self, atom_index: usize, rotamer: &Rotamer) -> i8 {
        let atom = &rotamer.atoms[atom_index];
        self.determine_hbond_role(&atom.force_field_type, || {
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

    /// Determines the hydrogen bonding role based on force field type and neighbor information.
    ///
    /// This generic method implements the logic for classifying atoms as hydrogen bond
    /// donors or acceptors based on their force field types and the types of their bonded
    /// neighbors. It handles the special case of donor hydrogens that require checking
    /// the heavy atom they are attached to.
    ///
    /// # Arguments
    ///
    /// * `ff_type` - The force field type of the atom
    /// * `get_neighbor_ff_type` - A closure that returns the force field type of the bonded neighbor
    ///
    /// # Return
    ///
    /// Returns an integer indicating the hydrogen bonding role:
    /// - `-1`: Not involved in hydrogen bonding
    /// - `0`: Hydrogen bond donor
    /// - `1`: Hydrogen bond acceptor
    fn determine_hbond_role<'b, F>(&self, ff_type: &str, get_neighbor_ff_type: F) -> i8
    where
        F: Fn() -> Option<&'b str>,
    {
        let is_donor_hydrogen = ff_type == DREIDING_HBOND_DONOR_HYDROGEN;

        let mut hbond_type_id = -1;

        if is_donor_hydrogen {
            if let Some(heavy_atom_ff_type) = get_neighbor_ff_type() {
                if self
                    .forcefield
                    .non_bonded
                    .hbond_donors
                    .contains(heavy_atom_ff_type)
                {
                    return 0;
                }
            }
        }

        if !is_donor_hydrogen && self.forcefield.non_bonded.hbond_acceptors.contains(ff_type) {
            hbond_type_id = 1;
        }

        hbond_type_id
    }
}

/// Assigns atom roles to a collection of atoms based on residue topology.
///
/// This function implements a multi-step algorithm to classify atoms as backbone or
/// sidechain based on topology definitions. It uses a pool-based approach to handle
/// cases where multiple atoms may have the same name, ensuring proper assignment
/// by consuming atoms from different ends of the pool for anchor vs sidechain atoms.
///
/// # Arguments
///
/// * `atoms` - Mutable slice of atoms to assign roles to
/// * `get_name` - Function to extract the atom name from an atom
/// * `set_role` - Function to set the role on an atom
/// * `topology` - The topology definition for the residue
/// * `residue_name` - The name of the residue for error reporting
///
/// # Return
///
/// Returns `Ok(())` if role assignment succeeds, or an error if required anchor atoms
/// are missing.
///
/// # Errors
///
/// Returns `ParameterizationError::InvalidAnchorAtom` if any required anchor atom
/// is missing from the atom collection.
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

/// Converts a force field VDW parameter into a cached parameter for efficient computation.
///
/// This implementation provides a direct conversion from the configuration-based
/// `VdwParam` enum to the runtime-optimized `CachedVdwParam` enum used during
/// energy calculations.
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
        forcefield::params::EnergyWeights,
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
            [hbond.O_2-H_O]
            equilibrium_distance = 2.8
            well_depth = 5.0
            [hbond.H_O-O_2]
            equilibrium_distance = 2.8
            well_depth = 5.0
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

        let forcefield =
            Forcefield::load(&ff_path, &delta_path, &EnergyWeights::default()).unwrap();

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
        fn parameterize_rotamer_assigns_physicochemical_and_hbond_fields() {
            let TestSetup {
                forcefield,
                topology_registry,
                ..
            } = setup();
            let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 2.0);
            let ala_topology = topology_registry.get("ALA").unwrap();

            let mut rotamer = Rotamer {
                atoms: vec![
                    Atom::new("N", ResidueId::default(), Default::default()),
                    Atom::new("CA", ResidueId::default(), Default::default()),
                    Atom::new("C", ResidueId::default(), Default::default()),
                    Atom::new("CB", ResidueId::default(), Default::default()),
                    Atom::new("O", ResidueId::default(), Default::default()),
                    Atom::new("H", ResidueId::default(), Default::default()),
                ],
                bonds: vec![(4, 5)],
            };

            let mappings = [
                ("N", "N_R"),
                ("CA", "C_BB"),
                ("C", "C_R"),
                ("CB", "C_SC"),
                ("O", "O_2"),
                ("H", "H___A"),
            ];
            for atom in rotamer.atoms.iter_mut() {
                for (name, ff) in &mappings {
                    if atom.name == *name {
                        atom.force_field_type = (*ff).to_string();
                        break;
                    }
                }
            }

            parameterizer
                .parameterize_rotamer(&mut rotamer, "ALA", ala_topology)
                .unwrap();

            let get = |n: &str| rotamer.atoms.iter().find(|a| a.name == n).unwrap();

            let cb = get("CB");
            assert_eq!(cb.role, AtomRole::Sidechain);
            assert!((cb.delta - 0.2).abs() < 1e-9, "CB delta incorrect");
            assert!(
                matches!(cb.vdw_param, CachedVdwParam::LennardJones { radius, .. } if (radius-1.9).abs()<1e-12)
            );
            assert_eq!(cb.hbond_type_id, -1, "CB should not be H-bond classified");

            let ca = get("CA");
            assert_eq!(ca.role, AtomRole::Backbone);
            assert_eq!(ca.delta, 0.0, "CA delta should be zero");
            assert!(
                matches!(ca.vdw_param, CachedVdwParam::LennardJones { radius, .. } if (radius-1.8).abs()<1e-12)
            );

            let h = get("H");
            assert_eq!(h.hbond_type_id, 0, "Hydrogen should be donor (0)");

            let o = get("O");
            assert_eq!(o.hbond_type_id, 1, "Oxygen should be acceptor (1)");
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
