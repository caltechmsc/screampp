use crate::core::models::atom::Atom;
use serde::Deserialize;

/// Represents atom data for a single atom in a rotamer, suitable for deserialization from external files.
///
/// This struct contains the essential information needed to reconstruct an atom
/// within a protein side-chain rotamer, including its position, charge, and
/// force field classification. It serves as an intermediate representation
/// between serialized rotamer data and the full `Atom` structure used in molecular systems.
#[derive(Debug, Clone, Deserialize)]
pub struct RotamerAtomData {
    /// The serial number of the atom within the rotamer.
    ///
    /// This corresponds to the atom's position in the rotamer's atom list
    /// and is used for establishing connectivity relationships.
    pub serial: usize,
    /// The name of the atom (e.g., "CA", "CB", "CG").
    ///
    /// This follows standard PDB atom naming conventions and helps
    /// identify the atom's role in the amino acid side chain.
    pub atom_name: String,
    /// The partial charge of the atom in atomic units.
    ///
    /// This value is used in electrostatic energy calculations and
    /// is typically derived from quantum mechanical calculations or
    /// empirical force field parameters.
    pub partial_charge: f64,
    /// The 3D coordinates of the atom in Cartesian space.
    ///
    /// Coordinates are stored as [x, y, z] in Angstroms and represent
    /// the atom's position relative to the rotamer's local coordinate system.
    pub position: [f64; 3],
    /// The force field atom type identifier.
    ///
    /// This string identifies the atom's type in the molecular mechanics
    /// force field (e.g., "C_3", "N_R", "O_2") and determines which
    /// parameters to use for energy calculations.
    pub force_field_type: String,
}

/// Represents the complete data for a protein side-chain rotamer, suitable for deserialization.
///
/// This struct encapsulates all the information needed to define a specific
/// conformational state of an amino acid side chain, including both atomic
/// properties and connectivity. Rotamers are discrete conformational states
/// that amino acid side chains can adopt, and this structure provides the
/// data needed to reconstruct them from external storage formats.
#[derive(Debug, Clone, Deserialize)]
pub struct RotamerData {
    /// The atoms that make up this rotamer.
    ///
    /// Each atom contains its name, position, charge, and force field type.
    /// The atoms are ordered by their serial numbers for consistent indexing.
    pub atoms: Vec<RotamerAtomData>,
    /// The bonds connecting atoms within the rotamer.
    ///
    /// Each bond is represented as a pair of atom indices [atom1, atom2],
    /// where the indices correspond to positions in the `atoms` vector.
    /// Only bonds within the side chain are included (backbone connectivity
    /// is handled separately).
    pub bonds: Vec<[usize; 2]>,
}

/// Represents a protein side-chain rotamer in its runtime form.
///
/// This struct contains the complete molecular representation of a specific
/// conformational state of an amino acid side chain, ready for use in
/// molecular modeling and energy calculations. Unlike `RotamerData`, this
/// structure uses the full `Atom` type from the molecular system, making
/// it suitable for integration with the broader molecular modeling framework.
#[derive(Debug, Clone)]
pub struct Rotamer {
    /// The atoms that constitute this rotamer.
    ///
    /// Each atom is a complete `Atom` object with all necessary properties
    /// for molecular modeling, including coordinates, charges, and force
    /// field parameters. The atoms are ordered consistently for reliable
    /// indexing and energy calculations.
    pub atoms: Vec<Atom>,
    /// The bonds connecting atoms within the rotamer.
    ///
    /// Each bond is represented as a pair of indices (atom1_idx, atom2_idx)
    /// that reference positions in the `atoms` vector. These define the
    /// connectivity of the side chain atoms and are essential for molecular
    /// mechanics calculations and structural analysis.
    pub bonds: Vec<(usize, usize)>,
}
