use super::rotamer::{Rotamer, RotamerData};
use crate::core::{
    forcefield::{
        parameterization::{ParameterizationError, Parameterizer},
        params::Forcefield,
    },
    models::{
        atom::Atom,
        ids::{AtomId, ResidueId},
        residue::ResidueType,
        system::MolecularSystem,
    },
    topology::registry::TopologyRegistry,
};
use nalgebra::Point3;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::str::FromStr;
use thiserror::Error;

/// Type alias for the raw rotamer data loaded from TOML files.
///
/// This represents the deserialized structure of rotamer library files,
/// mapping residue names to their corresponding rotamer data arrays.
/// Used internally during the loading process before conversion to
/// the runtime `Rotamer` format.
type RawRotamerFile = HashMap<String, Vec<RotamerData>>;

/// Manages a collection of pre-computed rotamer conformations for protein side-chain modeling.
///
/// This library serves as the central repository for rotamer data in SCREAM++,
/// providing access to discrete conformational states of amino acid side chains.
/// Rotamers are pre-computed low-energy conformations that represent the most
/// probable side-chain orientations, enabling efficient protein structure prediction
/// and refinement algorithms.
#[derive(Debug, Default, Clone)]
pub struct RotamerLibrary {
    /// The core storage mapping residue types to their available rotamers.
    ///
    /// Each residue type (e.g., Alanine, Glycine) is associated with a vector
    /// of possible rotamer conformations. These rotamers are fully parameterized
    /// and ready for use in molecular mechanics calculations.
    pub rotamers: HashMap<ResidueType, Vec<Rotamer>>,
}

/// Represents errors that can occur during rotamer library loading and processing.
///
/// This enum encompasses all possible failure modes when loading rotamer data
/// from external files, validating topology consistency, and parameterizing
/// rotamers for molecular mechanics calculations.
#[derive(Debug, Error)]
pub enum LibraryLoadError {
    /// Indicates that the rotamer library file could not be read from disk.
    ///
    /// This error occurs when there are permission issues, the file doesn't exist,
    /// or other I/O-related problems prevent reading the rotamer configuration.
    #[error("File I/O error for '{path}': {source}")]
    Io {
        /// The path to the file that could not be read.
        path: String,
        /// The underlying I/O error that occurred.
        source: std::io::Error,
    },
    /// Indicates that the rotamer library file content is not valid TOML.
    ///
    /// This error occurs when the file exists but cannot be parsed as TOML,
    /// or when the parsed data doesn't match the expected `RotamerData` structure.
    #[error("TOML parsing error for '{path}': {source}")]
    Toml {
        /// The path to the file that could not be parsed.
        path: String,
        /// The underlying TOML parsing error that occurred.
        source: toml::de::Error,
    },
    /// Indicates that an unknown residue type was encountered in the rotamer file.
    ///
    /// This error occurs when the rotamer library contains a residue name that
    /// cannot be mapped to a known `ResidueType` enum variant, indicating
    /// either a typo or missing residue type definition.
    #[error("Unknown residue type '{0}' found in library file")]
    UnknownResidueType(String),
    /// Indicates that a residue type in the rotamer library lacks a topology definition.
    ///
    /// This error occurs when rotamer data exists for a residue type, but no
    /// corresponding topology information is available in the topology registry.
    /// Topology definitions are required for proper atom classification and parameterization.
    #[error(
        "Missing topology definition for residue type '{0}', which is present in the rotamer library"
    )]
    MissingTopology(String),
    /// Indicates that parameterization of a rotamer failed during loading.
    ///
    /// This error occurs when the force field parameterization process encounters
    /// issues, such as missing parameters or invalid atom types, preventing
    /// the rotamer from being used in energy calculations.
    #[error(
        "Parameterization failed for residue '{residue_type}' in rotamer from file '{path}': {source}"
    )]
    Parameterization {
        /// The path to the rotamer file being processed.
        path: String,
        /// The residue type that failed parameterization.
        residue_type: String,
        /// The underlying parameterization error that occurred.
        source: ParameterizationError,
    },
    /// Indicates that a bond definition references a non-existent atom serial.
    ///
    /// This error occurs when a bond in the rotamer data refers to an atom
    /// serial number that doesn't exist in the atom list, indicating corrupted
    /// or malformed rotamer data.
    #[error(
        "Invalid bond definition in rotamer library for residue '{residue_type}': bond references non-existent atom serial '{serial}'"
    )]
    InvalidBondSerial { residue_type: String, serial: usize },
    /// Indicates that duplicate atom serial numbers were found in a rotamer definition.
    ///
    /// This error occurs when multiple atoms in the same rotamer have the same
    /// serial number, which would cause conflicts in bond definitions and
    /// atom indexing.
    #[error(
        "Duplicate atom serial '{serial}' found in rotamer definition for residue '{residue_type}'"
    )]
    DuplicateAtomSerial { residue_type: String, serial: usize },
}

impl RotamerLibrary {
    /// Loads and parameterizes a rotamer library from a TOML configuration file.
    ///
    /// This method performs a multi-phase process to load rotamer data, validate
    /// consistency with topology definitions, and parameterize rotamers for use
    /// in molecular mechanics calculations. The process includes error checking
    /// and cross-validation to ensure data integrity.
    ///
    /// # Arguments
    ///
    /// * `rotamer_toml_path` - Path to the TOML file containing rotamer definitions
    /// * `topology_registry` - Registry containing residue topology information
    /// * `forcefield` - Force field parameters for rotamer parameterization
    /// * `delta_s_factor` - Scaling factor for delta parameters in parameterization
    ///
    /// # Return
    ///
    /// Returns a fully loaded and parameterized `RotamerLibrary` ready for use,
    /// or an error if loading, validation, or parameterization fails.
    ///
    /// # Errors
    ///
    /// Returns `LibraryLoadError::Io` if the file cannot be read.
    /// Returns `LibraryLoadError::Toml` if the file is not valid TOML.
    /// Returns `LibraryLoadError::UnknownResidueType` if unknown residues are found.
    /// Returns `LibraryLoadError::MissingTopology` if topology definitions are missing.
    /// Returns `LibraryLoadError::Parameterization` if parameterization fails.
    /// Returns `LibraryLoadError::InvalidBondSerial` if bond definitions are invalid.
    /// Returns `LibraryLoadError::DuplicateAtomSerial` if duplicate serials exist.
    pub fn load(
        rotamer_toml_path: &Path,
        topology_registry: &TopologyRegistry,
        forcefield: &Forcefield,
        delta_s_factor: f64,
    ) -> Result<Self, LibraryLoadError> {
        // Phase 1: Load raw rotamer data from TOML file
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

        // Phase 2: Create parameterizer for pre-parameterizing rotamers
        let parameterizer = Parameterizer::new(forcefield, topology_registry, delta_s_factor);
        let mut final_rotamers_map = HashMap::new();

        // Phase 3: Process and parameterize each rotamer
        for (res_name, raw_rotamer_list) in raw_lib {
            let residue_type = ResidueType::from_str(&res_name)
                .map_err(|_| LibraryLoadError::UnknownResidueType(res_name.clone()))?;

            let topology = topology_registry
                .get(&res_name)
                .ok_or_else(|| LibraryLoadError::MissingTopology(res_name.clone()))?;

            let mut processed_rotamers = Vec::with_capacity(raw_rotamer_list.len());
            for raw_rotamer_data in raw_rotamer_list {
                let rotamer = Self::process_raw_rotamer(
                    &raw_rotamer_data,
                    &parameterizer,
                    &res_name,
                    topology,
                    rotamer_toml_path,
                )?;
                processed_rotamers.push(rotamer);
            }
            final_rotamers_map.insert(residue_type, processed_rotamers);
        }

        Ok(Self {
            rotamers: final_rotamers_map,
        })
    }

    /// Processes and parameterizes a single raw rotamer from the library file.
    ///
    /// This method converts raw rotamer data into a fully parameterized `Rotamer`
    /// structure suitable for molecular mechanics calculations. It handles atom
    /// creation, bond validation, and force field parameterization.
    ///
    /// # Arguments
    ///
    /// * `raw_rotamer_data` - The raw rotamer data from the TOML file
    /// * `parameterizer` - The parameterizer configured for this library
    /// * `res_name` - Name of the residue type being processed
    /// * `topology` - Topology definition for the residue
    /// * `path_for_error` - File path for error reporting
    ///
    /// # Return
    ///
    /// Returns a fully parameterized `Rotamer` ready for use, or an error
    /// if processing or parameterization fails.
    ///
    /// # Errors
    ///
    /// Returns `LibraryLoadError::InvalidBondSerial` if bond references invalid atoms.
    /// Returns `LibraryLoadError::DuplicateAtomSerial` if duplicate serials exist.
    /// Returns `LibraryLoadError::Parameterization` if parameterization fails.
    fn process_raw_rotamer(
        raw_rotamer_data: &RotamerData,
        parameterizer: &Parameterizer,
        res_name: &str,
        topology: &crate::core::topology::registry::ResidueTopology,
        path_for_error: &Path,
    ) -> Result<Rotamer, LibraryLoadError> {
        let mut atoms = Vec::with_capacity(raw_rotamer_data.atoms.len());
        let mut serial_to_index_map = HashMap::with_capacity(raw_rotamer_data.atoms.len());
        let placeholder_residue_id = ResidueId::default();

        for (index, atom_data) in raw_rotamer_data.atoms.iter().enumerate() {
            if serial_to_index_map
                .insert(atom_data.serial, index)
                .is_some()
            {
                return Err(LibraryLoadError::DuplicateAtomSerial {
                    residue_type: res_name.to_string(),
                    serial: atom_data.serial,
                });
            }

            let mut atom = Atom::new(
                &atom_data.atom_name,
                placeholder_residue_id,
                Point3::from(atom_data.position),
            );
            atom.partial_charge = atom_data.partial_charge;
            atom.force_field_type = atom_data.force_field_type.clone();

            atoms.push(atom);
        }

        let mut bonds = Vec::with_capacity(raw_rotamer_data.bonds.len());
        for bond_serials in &raw_rotamer_data.bonds {
            let index1 = *serial_to_index_map.get(&bond_serials[0]).ok_or_else(|| {
                LibraryLoadError::InvalidBondSerial {
                    residue_type: res_name.to_string(),
                    serial: bond_serials[0],
                }
            })?;
            let index2 = *serial_to_index_map.get(&bond_serials[1]).ok_or_else(|| {
                LibraryLoadError::InvalidBondSerial {
                    residue_type: res_name.to_string(),
                    serial: bond_serials[1],
                }
            })?;
            bonds.push((index1, index2));
        }

        let mut rotamer = Rotamer { atoms, bonds };

        parameterizer
            .parameterize_rotamer(&mut rotamer, res_name, topology)
            .map_err(|e| LibraryLoadError::Parameterization {
                path: path_for_error.to_string_lossy().to_string(),
                residue_type: res_name.to_string(),
                source: e,
            })?;

        Ok(rotamer)
    }

    /// Retrieves all available rotamers for a specific residue type.
    ///
    /// This method provides access to the pre-computed rotamer conformations
    /// for a given amino acid type, enabling side-chain placement algorithms
    /// to evaluate different conformational possibilities.
    ///
    /// # Arguments
    ///
    /// * `residue_type` - The type of amino acid residue to get rotamers for
    ///
    /// # Return
    ///
    /// Returns `Some(&Vec<Rotamer>)` containing all rotamers for the residue type,
    /// or `None` if no rotamers are available for the requested type.
    pub fn get_rotamers_for(&self, residue_type: ResidueType) -> Option<&Vec<Rotamer>> {
        self.rotamers.get(&residue_type)
    }

    /// Extracts and includes rotamer conformations from existing molecular system residues.
    ///
    /// This method allows incorporating experimentally determined or user-provided
    /// conformations from a molecular system into the rotamer library. It extracts
    /// side-chain conformations from specified residues and adds them as additional
    /// rotamer options for the corresponding residue types.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system containing the residues to extract from
    /// * `active_residues` - Set of residue IDs to extract conformations from
    /// * `topology_registry` - Registry for residue topology information
    /// * `parameterizer` - Parameterizer for force field parameterization
    pub fn include_system_conformations(
        &mut self,
        system: &MolecularSystem,
        active_residues: &HashSet<ResidueId>,
        topology_registry: &TopologyRegistry,
        parameterizer: &Parameterizer,
    ) {
        for &residue_id in active_residues {
            if let Some(mut rotamer) =
                self.extract_rotamer_from_system(system, residue_id, topology_registry)
            {
                let residue = system.residue(residue_id).unwrap();
                let residue_type = residue.residue_type.unwrap();
                let topology = topology_registry.get(&residue.name).unwrap();

                if let Err(e) =
                    parameterizer.parameterize_rotamer(&mut rotamer, &residue.name, topology)
                {
                    tracing::error!(
                        "Failed to parameterize system-extracted rotamer for residue {} {}: {}",
                        residue.residue_number,
                        residue.name,
                        e
                    );
                    continue;
                }

                self.rotamers.entry(residue_type).or_default().push(rotamer);
            }
        }
    }

    /// Extracts a complete rotamer conformation from a molecular system residue.
    ///
    /// This method reconstructs a rotamer from the atoms and bonds of a specific
    /// residue in a molecular system, following the topology definition to ensure
    /// all required atoms are present and properly classified.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system containing the residue
    /// * `residue_id` - ID of the residue to extract the rotamer from
    /// * `topology_registry` - Registry containing topology definitions
    ///
    /// # Return
    ///
    /// Returns `Some(Rotamer)` if extraction succeeds, or `None` if the residue
    /// is missing required atoms or topology information.
    fn extract_rotamer_from_system(
        &self,
        system: &MolecularSystem,
        residue_id: ResidueId,
        topology_registry: &TopologyRegistry,
    ) -> Option<Rotamer> {
        let residue = system.residue(residue_id)?;
        let topology = topology_registry.get(&residue.name)?;

        // Step 1: Build the consuming pool of atoms from the source residue
        let mut atom_pool: HashMap<String, Vec<AtomId>> = HashMap::new();
        for &atom_id in residue.atoms() {
            if let Some(atom) = system.atom(atom_id) {
                atom_pool
                    .entry(atom.name.clone())
                    .or_default()
                    .push(atom_id);
            }
        }

        let mut extracted_atoms = Vec::new();
        let mut old_id_to_new_index = HashMap::new();
        let mut consumed_atom_ids = HashSet::new();

        // Step 2: Extract ANCHOR atoms (Mandatory)
        for anchor_name in &topology.anchor_atoms {
            match atom_pool.get_mut(anchor_name) {
                Some(ids) if !ids.is_empty() => {
                    let atom_id = ids.remove(0);
                    if consumed_atom_ids.insert(atom_id) {
                        let atom = system.atom(atom_id)?;
                        old_id_to_new_index.insert(atom_id, extracted_atoms.len());
                        extracted_atoms.push(atom.clone());
                    }
                }
                _ => {
                    tracing::warn!(
                        "Skipping conformation extraction for residue {} {} (ID: {:?}): Missing mandatory ANCHOR atom '{}'.",
                        residue.residue_number,
                        residue.name,
                        residue_id,
                        anchor_name
                    );
                    return None;
                }
            }
        }

        // Step 3: Extract SIDECHAIN atoms (Mandatory)
        for sidechain_name in &topology.sidechain_atoms {
            match atom_pool.get_mut(sidechain_name) {
                Some(ids) if !ids.is_empty() => {
                    let atom_id = ids.pop().unwrap();
                    if consumed_atom_ids.insert(atom_id) {
                        let atom = system.atom(atom_id)?;
                        old_id_to_new_index.insert(atom_id, extracted_atoms.len());
                        extracted_atoms.push(atom.clone());
                    }
                }
                _ => {
                    tracing::warn!(
                        "Skipping conformation extraction for residue {} {} (ID: {:?}): Incomplete sidechain, missing atom '{}'.",
                        residue.residue_number,
                        residue.name,
                        residue_id,
                        sidechain_name
                    );
                    return None;
                }
            }
        }

        // Step 4: Reconstruct bonds for the new, complete rotamer
        let mut new_rotamer_bonds = Vec::new();
        for (&old_id_a, &new_idx_a) in &old_id_to_new_index {
            if let Some(neighbors) = system.get_bonded_neighbors(old_id_a) {
                for &old_id_b in neighbors {
                    if let Some(&new_idx_b) = old_id_to_new_index.get(&old_id_b) {
                        if new_idx_a < new_idx_b {
                            new_rotamer_bonds.push((new_idx_a, new_idx_b));
                        }
                    }
                }
            }
        }

        Some(Rotamer {
            atoms: extracted_atoms,
            bonds: new_rotamer_bonds,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::parameterization::Parameterizer;
    use crate::core::forcefield::params::{DeltaParam, GlobalParams, NonBondedParams, VdwParam};
    use crate::core::{
        models::{
            atom::{Atom, AtomRole},
            chain::ChainType,
            residue::ResidueType,
        },
        topology::registry::TopologyRegistry,
    };
    use std::fs::File;
    use std::io::Write;
    use tempfile::TempDir;

    struct TestSetup {
        forcefield: Forcefield,
        topology_registry: TopologyRegistry,
        temp_dir: TempDir,
    }

    fn setup() -> TestSetup {
        let temp_dir = TempDir::new().unwrap();

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
        let forcefield = Forcefield {
            deltas,
            non_bonded: NonBondedParams {
                globals: GlobalParams {
                    dielectric_constant: 1.0,
                    potential_function: "".to_string(),
                },
                vdw,
                hbond: HashMap::new(),
                hbond_donors: HashSet::new(),
                hbond_acceptors: HashSet::new(),
            },
            weight_map: HashMap::new(),
        };

        let topology_toml_path = temp_dir.path().join("topology.toml");
        let mut topology_file = File::create(&topology_toml_path).unwrap();
        write!(
            topology_file,
            r#"
[ALA]
anchor_atoms = ["N", "CA", "C"]
sidechain_atoms = ["CB", "HB1", "HB2", "HB3"]
[GLY]
anchor_atoms = ["N", "CA", "C", "HA2"]
sidechain_atoms = ["HA1"]
[TRP]
anchor_atoms = []
sidechain_atoms = []
"#
        )
        .unwrap();
        let topology_registry = TopologyRegistry::load(&topology_toml_path).unwrap();

        TestSetup {
            forcefield,
            topology_registry,
            temp_dir,
        }
    }

    fn write_rotamer_file(dir: &Path, name: &str, content: &str) -> std::path::PathBuf {
        let path = dir.join(name);
        let mut file = File::create(&path).unwrap();
        write!(file, "{}", content).unwrap();
        path
    }

    mod load_tests {
        use super::*;

        #[test]
        fn load_and_parameterize_rotamers_successfully() {
            let setup = setup();
            let rotamer_content = r#"
[[ALA]]
atoms = [
    { serial = 1, atom_name = "N", position = [0.0, 1.0, 0.0], partial_charge = -0.3, force_field_type = "C_BB" },
    { serial = 2, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.1, force_field_type = "C_BB" },
    { serial = 3, atom_name = "C", position = [1.0, 0.0, 0.0], partial_charge = 0.5, force_field_type = "C_BB" },
    { serial = 4, atom_name = "CB", position = [-1.0, -0.5, 0.0], partial_charge = -0.3, force_field_type = "C_SC" }
]
bonds = [ [2, 1], [2, 3], [2, 4] ]
"#;
            let rotamer_path =
                write_rotamer_file(&setup.temp_dir.path(), "rot.toml", rotamer_content);

            let library = RotamerLibrary::load(
                &rotamer_path,
                &setup.topology_registry,
                &setup.forcefield,
                1.0,
            )
            .unwrap();

            assert!(library.get_rotamers_for(ResidueType::Alanine).is_some());
            let rotamer = &library.get_rotamers_for(ResidueType::Alanine).unwrap()[0];
            assert_eq!(
                rotamer.atoms.iter().find(|a| a.name == "CB").unwrap().role,
                AtomRole::Sidechain
            );
            assert_eq!(
                rotamer.atoms.iter().find(|a| a.name == "CA").unwrap().role,
                AtomRole::Backbone
            );
        }

        #[test]
        fn load_fails_if_topology_is_missing() {
            let setup = setup();
            let rotamer_content = r#"
[[LYS]]
atoms = []
bonds = []
"#;
            let rotamer_path =
                write_rotamer_file(&setup.temp_dir.path(), "lys_rotamer.toml", rotamer_content);
            let result = RotamerLibrary::load(
                &rotamer_path,
                &setup.topology_registry,
                &setup.forcefield,
                0.0,
            );

            assert!(
                matches!(result, Err(LibraryLoadError::MissingTopology(name)) if name == "LYS"),
                "Expected MissingTopology error for LYS, but got something else or Ok."
            );
        }

        #[test]
        fn load_fails_on_duplicate_atom_serial() {
            let setup = setup();
            let rotamer_content = r#"
[[ALA]]
atoms = [
    { serial = 1, atom_name = "N", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_BB" },
    { serial = 1, atom_name = "CA", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_BB" }
]
bonds = []"#;
            let rotamer_path =
                write_rotamer_file(&setup.temp_dir.path(), "rot.toml", rotamer_content);
            let result = RotamerLibrary::load(
                &rotamer_path,
                &setup.topology_registry,
                &setup.forcefield,
                1.0,
            );

            assert!(matches!(
                result,
                Err(LibraryLoadError::DuplicateAtomSerial { .. })
            ));
        }

        #[test]
        fn load_fails_on_invalid_bond_serial() {
            let setup = setup();
            let rotamer_content = r#"
[[ALA]]
atoms = [{ serial = 1, atom_name = "N", position = [0.0, 0.0, 0.0], partial_charge = 0.0, force_field_type = "C_BB" }]
bonds = [[1, 99]]"#;
            let rotamer_path =
                write_rotamer_file(&setup.temp_dir.path(), "rot.toml", rotamer_content);
            let result = RotamerLibrary::load(
                &rotamer_path,
                &setup.topology_registry,
                &setup.forcefield,
                1.0,
            );

            assert!(
                matches!(result, Err(LibraryLoadError::InvalidBondSerial { serial, .. }) if serial == 99)
            );
        }
    }

    mod include_system_conformations_tests {
        use super::*;

        #[test]
        fn extracts_standard_residue_ala() {
            let setup = setup();
            let mut library = RotamerLibrary::default();
            let mut system = MolecularSystem::new();
            let chain_id = system.add_chain('A', ChainType::Protein);
            let ala_id = system
                .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
                .unwrap();
            for name in &["N", "CA", "C", "CB", "HB1", "HB2", "HB3"] {
                system
                    .add_atom_to_residue(ala_id, Atom::new(name, ala_id, Default::default()))
                    .unwrap();
            }
            let active = [ala_id].iter().cloned().collect();

            let parameterizer =
                Parameterizer::new(&setup.forcefield, &setup.topology_registry, 1.0);
            library.include_system_conformations(
                &system,
                &active,
                &setup.topology_registry,
                &parameterizer,
            );

            let rotamers = library.get_rotamers_for(ResidueType::Alanine).unwrap();
            assert_eq!(rotamers.len(), 1);
            assert_eq!(rotamers[0].atoms.len(), 7);
        }

        #[test]
        fn extracts_glycine_with_complex_topology() {
            let setup = setup();
            let mut library = RotamerLibrary::default();
            let mut system = MolecularSystem::new();
            let chain_id = system.add_chain('A', ChainType::Protein);
            let gly_id = system
                .add_residue(chain_id, 1, "GLY", Some(ResidueType::Glycine))
                .unwrap();
            for name in &["N", "CA", "C", "HA1", "HA2"] {
                system
                    .add_atom_to_residue(gly_id, Atom::new(name, gly_id, Default::default()))
                    .unwrap();
            }
            let active = [gly_id].iter().cloned().collect();

            let parameterizer =
                Parameterizer::new(&setup.forcefield, &setup.topology_registry, 1.0);
            library.include_system_conformations(
                &system,
                &active,
                &setup.topology_registry,
                &parameterizer,
            );

            let rotamers = library.get_rotamers_for(ResidueType::Glycine).unwrap();
            assert_eq!(rotamers.len(), 1);
            assert_eq!(rotamers[0].atoms.len(), 5);
        }

        #[test]
        fn include_is_safe_if_called_multiple_times() {
            let setup = setup();
            let mut library = RotamerLibrary::default();
            let mut system = MolecularSystem::new();
            let chain_id = system.add_chain('A', ChainType::Protein);
            let ala_id = system
                .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
                .unwrap();
            for name in &["N", "CA", "C", "CB", "HB1", "HB2", "HB3"] {
                system
                    .add_atom_to_residue(ala_id, Atom::new(name, ala_id, Default::default()))
                    .unwrap();
            }
            let active = [ala_id].iter().cloned().collect();

            let parameterizer =
                Parameterizer::new(&setup.forcefield, &setup.topology_registry, 1.0);
            library.include_system_conformations(
                &system,
                &active,
                &setup.topology_registry,
                &parameterizer,
            );
            assert_eq!(
                library
                    .get_rotamers_for(ResidueType::Alanine)
                    .unwrap()
                    .len(),
                1
            );

            library.include_system_conformations(
                &system,
                &active,
                &setup.topology_registry,
                &parameterizer,
            );
            assert_eq!(
                library
                    .get_rotamers_for(ResidueType::Alanine)
                    .unwrap()
                    .len(),
                2
            );
        }

        #[test]
        fn include_skips_residue_if_topology_is_missing() {
            let setup = setup();
            let mut library = RotamerLibrary::default();
            let mut system = MolecularSystem::new();
            let chain_id = system.add_chain('A', ChainType::Protein);
            let pro_id = system
                .add_residue(chain_id, 1, "PRO", Some(ResidueType::Proline))
                .unwrap();
            let active = [pro_id].iter().cloned().collect();

            let parameterizer =
                Parameterizer::new(&setup.forcefield, &setup.topology_registry, 1.0);
            library.include_system_conformations(
                &system,
                &active,
                &setup.topology_registry,
                &parameterizer,
            );

            assert!(library.get_rotamers_for(ResidueType::Proline).is_none());
        }
    }
}
