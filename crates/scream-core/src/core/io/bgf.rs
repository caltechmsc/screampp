use super::sorting::sorter::sort_system_atoms;
use super::traits::MolecularFile;
use crate::core::models::atom::Atom;
use crate::core::models::chain::{Chain, ChainType};
use crate::core::models::ids::{AtomId, ChainId, ResidueId};
use crate::core::models::residue::{Residue, ResidueType};
use crate::core::models::system::MolecularSystem;
use crate::core::models::topology::BondOrder;
use nalgebra::Point3;
use std::collections::{BTreeSet, HashMap};
use std::io::{self, BufRead, Write};
use thiserror::Error;

/// Metadata associated with a BGF file, containing header information and other non-structural data.
///
/// This struct holds information that is not part of the molecular system's topology but is
/// preserved from the original BGF file, such as header lines that may contain remarks or
/// force field information.
#[derive(Debug, Default, Clone)]
pub struct BgfMetadata {
    /// Lines from the BGF file header that appear before the atom records.
    ///
    /// These typically include remarks, force field specifications, and other metadata
    /// that should be preserved when writing the file back out.
    pub header_lines: Vec<String>,
}

/// A handler for reading and writing BGF (Biograf) molecular structure files.
///
/// This struct implements the `MolecularFile` trait to provide support for the BGF format,
/// which is commonly used in molecular simulations and modeling software. It handles
/// parsing atom records, connectivity information, and metadata while building a
/// `MolecularSystem` representation.
#[derive(Debug, Default)]
pub struct BgfFile;

/// Errors that can occur during BGF file processing.
///
/// This enum covers I/O errors, parsing failures, and logical inconsistencies
/// encountered when reading or writing BGF files.
#[derive(Debug, Error)]
pub enum BgfError {
    /// An I/O error occurred during file operations.
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    /// A parsing error occurred on a specific line.
    #[error("Parse error on line {line_num}: {message}")]
    Parse { line_num: usize, message: String },
    /// A logical error occurred during file processing.
    #[error("Logic error during file processing: {0}")]
    Logic(String),
}

impl MolecularFile for BgfFile {
    type Metadata = BgfMetadata;
    type Error = BgfError;

    /// Reads a molecular system from a BGF file.
    ///
    /// This method parses the BGF file format, extracting atoms, residues, chains,
    /// and connectivity information to construct a `MolecularSystem`. It handles
    /// both ATOM and HETATM records, automatically categorizing chains based on
    /// residue types, and processes CONECT records to establish bonds.
    ///
    /// # Arguments
    ///
    /// * `reader` - A buffered reader providing the BGF file content.
    ///
    /// # Return
    ///
    /// Returns a tuple containing the constructed `MolecularSystem` and any
    /// associated metadata extracted from the file.
    ///
    /// # Errors
    ///
    /// Returns a `BgfError` if parsing fails due to malformed input, I/O issues,
    /// or logical inconsistencies in the file structure.
    fn read_from(
        reader: &mut impl BufRead,
    ) -> Result<(MolecularSystem, Self::Metadata), Self::Error> {
        let mut system = MolecularSystem::new();
        let mut metadata = BgfMetadata::default();

        let mut serial_to_atom_id = HashMap::<usize, AtomId>::new();

        // --- Phase 1: Pre-scan and categorize all lines ---
        let lines: Vec<_> = reader.lines().collect::<Result<_, _>>()?;
        let mut atom_lines = Vec::new();
        let mut conect_lines = Vec::new();
        let mut atom_section_started = false;

        for (i, line) in lines.iter().enumerate() {
            let line_num = i + 1;
            if line.trim().is_empty() || line.trim() == "END" {
                continue;
            }

            let record = line.get(0..6).unwrap_or("").trim();
            match record {
                "ATOM" | "HETATM" => {
                    atom_section_started = true;
                    atom_lines.push((line_num, line));
                }
                "CONECT" => conect_lines.push((line_num, line)),
                "ORDER" => {}
                "FORMAT" => {}
                _ => {
                    if !atom_section_started {
                        metadata.header_lines.push(line.clone());
                    }
                }
            }
        }

        // --- Phase 2: Process ATOM/HETATM records to build the system topology ---
        let mut current_chain_id: Option<ChainId> = None;
        let mut current_residue_id: Option<ResidueId> = None;
        let mut last_chain_char = '\0';
        let mut last_res_seq = isize::MIN;

        for (line_num, line) in atom_lines {
            let atom_info = parse_atom_line(line).map_err(|msg| BgfError::Parse {
                line_num,
                message: msg,
            })?;

            if last_chain_char != atom_info.chain_char {
                let chain_type = match atom_info.record_type.as_str() {
                    "ATOM" => ChainType::Protein,
                    "HETATM" => match atom_info.res_name.as_str() {
                        "HOH" | "TIP3" => ChainType::Water,
                        _ => ChainType::Ligand,
                    },
                    _ => ChainType::Other,
                };
                current_chain_id = Some(system.add_chain(atom_info.chain_char, chain_type));
                last_chain_char = atom_info.chain_char;
                last_res_seq = isize::MIN;
            }

            if last_res_seq != atom_info.res_seq {
                current_residue_id = system.add_residue(
                    current_chain_id.unwrap(),
                    atom_info.res_seq,
                    &atom_info.res_name,
                    ResidueType::from_str_optional(&atom_info.res_name),
                );
                last_res_seq = atom_info.res_seq;
            }

            let active_residue_id = current_residue_id
                .ok_or_else(|| BgfError::Logic("Failed to create or find residue".to_string()))?;

            let mut atom = Atom::new(&atom_info.name, active_residue_id, atom_info.pos);

            atom.partial_charge = atom_info.charge;
            atom.force_field_type = atom_info.ff_type;

            let atom_id = system
                .add_atom_to_residue(active_residue_id, atom)
                .ok_or_else(|| BgfError::Logic("Failed to add atom".to_string()))?;

            serial_to_atom_id.insert(atom_info.serial, atom_id);
        }

        // --- Phase 3: Process CONECT records using the temporary serial map ---
        for (line_num, line) in conect_lines {
            let (base_serial, connected_serials) =
                parse_conect_line(line).map_err(|msg| BgfError::Parse {
                    line_num,
                    message: msg,
                })?;

            let base_id = serial_to_atom_id
                .get(&base_serial)
                .ok_or_else(|| BgfError::Parse {
                    line_num,
                    message: format!("Atom serial {} not found for CONECT record", base_serial),
                })?;

            for &connected_serial in &connected_serials {
                let connected_id =
                    serial_to_atom_id
                        .get(&connected_serial)
                        .ok_or_else(|| BgfError::Parse {
                            line_num,
                            message: format!(
                                "Atom serial {} not found for CONECT record",
                                connected_serial
                            ),
                        })?;
                system.add_bond(*base_id, *connected_id, BondOrder::Single);
            }
        }

        Ok((system, metadata))
    }

    /// Writes a molecular system to a BGF file with associated metadata.
    ///
    /// This method serializes the `MolecularSystem` into the BGF format, including
    /// atom records, connectivity information, and any provided metadata. Atoms are
    /// sorted into canonical order for consistent output, and serial numbers are
    /// reassigned sequentially.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system to write.
    /// * `metadata` - Metadata to include in the output file.
    /// * `writer` - A writer to output the BGF file content.
    ///
    /// # Return
    ///
    /// Returns `Ok(())` on successful writing.
    ///
    /// # Errors
    ///
    /// Returns a `BgfError` if writing fails due to I/O issues.
    fn write_to(
        system: &MolecularSystem,
        metadata: &Self::Metadata,
        writer: &mut impl Write,
    ) -> Result<(), Self::Error> {
        // --- Phase 1: Sort all atoms into a canonical order ---
        let sorted_atoms = sort_system_atoms(system);

        // --- Phase 2: Renumber atoms and create a map for CONECT records ---
        let mut atom_id_to_new_serial = HashMap::<AtomId, usize>::new();
        for (index, canonical_atom) in sorted_atoms.iter().enumerate() {
            atom_id_to_new_serial.insert(canonical_atom.id, index + 1);
        }

        // --- Phase 3: Write header and formatted ATOM records ---
        for line in &metadata.header_lines {
            writeln!(writer, "{}", line)?;
        }

        writeln!(
            writer,
            "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)"
        )?;

        for canonical_atom in &sorted_atoms {
            let atom = canonical_atom.source;
            let new_serial = atom_id_to_new_serial[&canonical_atom.id];

            let residue = system.residue(atom.residue_id).unwrap();
            let chain = system.chain(residue.chain_id).unwrap();

            let record_type = match chain.chain_type {
                ChainType::Protein | ChainType::DNA | ChainType::RNA => "ATOM  ",
                _ => "HETATM",
            };

            let atoms_connected = system
                .get_bonded_neighbors(canonical_atom.id)
                .map_or(0, |n| n.len() as u8);
            let lone_pairs = 0; // TODO: Replace with actual logic (placeholder)

            let atom_line = format_atom_line(
                record_type,
                new_serial,
                atom,
                residue,
                chain,
                atoms_connected,
                lone_pairs,
            );
            writeln!(writer, "{}", atom_line)?;
        }

        // --- Phase 4: Write CONECT records using the new serial numbers ---
        let mut bond_pairs = BTreeSet::new();
        for bond in system.bonds() {
            if let (Some(&serial1), Some(&serial2)) = (
                atom_id_to_new_serial.get(&bond.atom1_id),
                atom_id_to_new_serial.get(&bond.atom2_id),
            ) {
                let bond_pair = if serial1 < serial2 {
                    (serial1, serial2)
                } else {
                    (serial2, serial1)
                };
                bond_pairs.insert(bond_pair);
            }
        }

        writeln!(writer, "FORMAT CONECT (a6,12i6)")?;

        for (s1, s2) in bond_pairs {
            writeln!(writer, "CONECT {:>5} {:>6}", s1, s2)?;
        }

        writeln!(writer, "END")?;
        Ok(())
    }

    /// Writes a molecular system to a BGF file with default metadata.
    ///
    /// This is a convenience method that writes the system using minimal default
    /// metadata, including standard BGF header lines for BIOGRF version and
    /// force field information.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system to write.
    /// * `writer` - A writer to output the BGF file content.
    ///
    /// # Return
    ///
    /// Returns `Ok(())` on successful writing.
    ///
    /// # Errors
    ///
    /// Returns a `BgfError` if writing fails due to I/O issues.
    fn write_system_to(
        system: &MolecularSystem,
        writer: &mut impl Write,
    ) -> Result<(), Self::Error> {
        let minimal_metadata = BgfMetadata {
            header_lines: vec!["BIOGRF  332".to_string(), "FORCEFIELD DREIDING".to_string()],
            ..Default::default()
        };
        Self::write_to(system, &minimal_metadata, writer)
    }
}

/// Parses a single ATOM or HETATM line from a BGF file.
///
/// This function extracts all relevant information from a BGF atom record,
/// including position, charge, and force field type, returning a structured
/// representation suitable for building the molecular system.
///
/// # Arguments
///
/// * `line` - The BGF atom record line to parse.
///
/// # Return
///
/// Returns a `ParsedAtomInfo` struct containing the parsed data.
///
/// # Errors
///
/// Returns a `String` describing the parsing error if the line is malformed
/// or contains invalid data.
fn parse_atom_line(line: &str) -> Result<ParsedAtomInfo, String> {
    let get_slice = |start, end| {
        line.get(start..end)
            .ok_or_else(|| format!("Line is too short for slice {}-{}", start, end))
    };

    let record_type = get_slice(0, 6)?.trim().to_string();
    let serial = get_slice(7, 12)?
        .trim()
        .parse()
        .map_err(|e| format!("Invalid serial: {}", e))?;
    let name = get_slice(13, 18)?.trim().to_string();
    let res_name = get_slice(19, 22)?.trim().to_string();
    let chain_char = get_slice(23, 24)?.chars().next().unwrap_or(' ');
    let res_seq = get_slice(25, 30)?
        .trim()
        .parse()
        .map_err(|e| format!("Invalid residue number: {}", e))?;
    let x = get_slice(30, 40)?
        .trim()
        .parse()
        .map_err(|e| format!("Invalid X coordinate: {}", e))?;
    let y = get_slice(40, 50)?
        .trim()
        .parse()
        .map_err(|e| format!("Invalid Y coordinate: {}", e))?;
    let z = get_slice(50, 60)?
        .trim()
        .parse()
        .map_err(|e| format!("Invalid Z coordinate: {}", e))?;
    let ff_type = get_slice(61, 66)?.trim().to_string();
    let charge = get_slice(72, 80)?.trim().parse().unwrap_or(0.0);

    Ok(ParsedAtomInfo {
        record_type,
        serial,
        name,
        res_name,
        chain_char,
        res_seq,
        pos: Point3::new(x, y, z),
        charge,
        ff_type,
    })
}

/// Parses a CONECT line from a BGF file to extract connectivity information.
///
/// This function processes CONECT records, which specify bonds between atoms
/// using their serial numbers, returning the base atom and its connected atoms.
///
/// # Arguments
///
/// * `line` - The BGF CONECT record line to parse.
///
/// # Return
///
/// Returns a tuple containing the base atom serial and a vector of connected
/// atom serials.
///
/// # Errors
///
/// Returns a `String` describing the parsing error if the line is malformed
/// or contains invalid serial numbers.
fn parse_conect_line(line: &str) -> Result<(usize, Vec<usize>), String> {
    let mut parts = line.split_whitespace().skip(1);
    let base_serial_str = parts.next().ok_or("Missing base serial in CONECT")?;
    let base_serial = base_serial_str
        .parse()
        .map_err(|e| format!("Invalid base serial '{}': {}", base_serial_str, e))?;

    let connected_serials = parts
        .map(|s| s.parse::<usize>())
        .collect::<Result<Vec<usize>, _>>()
        .map_err(|e| format!("Invalid connected serial: {}", e))?;

    Ok((base_serial, connected_serials))
}

/// Formats an atom record line for output in BGF format.
///
/// This function constructs a properly formatted BGF atom line from the
/// provided atom, residue, and chain information, including connectivity
/// and charge data.
///
/// # Arguments
///
/// * `record_type` - The record type ("ATOM" or "HETATM").
/// * `serial` - The sequential serial number for the atom.
/// * `atom` - The atom to format.
/// * `residue` - The residue containing the atom.
/// * `chain` - The chain containing the residue.
/// * `atoms_connected` - Number of atoms connected to this atom.
/// * `lone_pairs` - Number of lone pairs (currently placeholder).
///
/// # Return
///
/// Returns the formatted BGF atom line as a string.
fn format_atom_line(
    record_type: &str,
    serial: usize,
    atom: &Atom,
    residue: &Residue,
    chain: &Chain,
    atoms_connected: u8,
    lone_pairs: u8,
) -> String {
    format!(
        "{:<6} {:>5} {:<5} {:>3} {:1} {:>5}{:>10.5}{:>10.5}{:>10.5} {:<5}{:>3}{:>2} {:>8.5}",
        record_type,
        serial,
        atom.name,
        residue.name,
        chain.id,
        residue.residue_number,
        atom.position.x,
        atom.position.y,
        atom.position.z,
        atom.force_field_type,
        atoms_connected,
        lone_pairs,
        atom.partial_charge
    )
}

/// Internal structure holding parsed information from a BGF atom line.
///
/// This struct is used internally during parsing to temporarily store
/// atom information before creating the actual `Atom` instance.
#[derive(Debug, Clone)]
struct ParsedAtomInfo {
    /// The record type ("ATOM" or "HETATM").
    record_type: String,
    /// The original serial number from the file.
    serial: usize,
    /// The atom name.
    name: String,
    /// The residue name.
    res_name: String,
    /// The chain identifier character.
    chain_char: char,
    /// The residue sequence number.
    res_seq: isize,
    /// The 3D position of the atom.
    pos: Point3<f64>,
    /// The partial charge of the atom.
    charge: f64,
    /// The force field type of the atom.
    ff_type: String,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::chain::ChainType;
    use crate::core::models::residue::ResidueType;
    use std::io::Cursor;

    const CANONICAL_BGF_DATA: &str = r#"BIOGRF 332
REMARK A standard, clean BGF file
FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)
ATOM       1 N     GLY A     1  -0.41600  -0.53500   0.00000 N_R    1 3 -0.35000
ATOM       2 CA    GLY A     1   0.00000   0.00000   1.00000 C_3    1 4  0.07000
ATOM       3 C     GLY A     1   1.00000   0.00000   0.00000 C_R    1 3  0.51000
ATOM       4 O     GLY A     1   1.50000   1.00000   0.00000 O_2    1 1 -0.51000
ATOM       5 N     ALA B     2   2.00000   0.00000   0.00000 N_R    1 3 -0.35000
ATOM       6 CA    ALA B     2   3.00000   0.00000   1.00000 C_3    1 4  0.07000
ATOM       7 CB    ALA B     2   3.50000   1.50000   1.00000 C_3    1 4 -0.18000
FORMAT CONECT (a6,12i6)
CONECT     1     2
CONECT     2     3
CONECT     3     4
CONECT     3     5
CONECT     5     6
CONECT     6     7
END
"#;

    const DISORDERED_BGF_DATA: &str = r#"BIOGRF 332
REMARK A messy BGF file with non-sequential serials and disordered records.
HETATM   999 O     HOH W     3      10.0      10.0      10.0 O_3    0 2   -0.834
ATOM      10 CA    ALA A     2       1.0       1.0       1.0 C_3    0 0    0.070
ATOM       5 N     ALA A     2       0.0       0.0       0.0 N_R    0 0   -0.350
ATOM     150 C     GLY A     1       3.0       3.0       3.0 C_R    0 0    0.510
ATOM       1 N     GLY A     1       2.0       2.0       2.0 N_R    0 0   -0.350
CONECT    10     5
CONECT     1   150
END
"#;

    #[test]
    fn read_from_parses_canonical_file_correctly() {
        let mut reader = Cursor::new(CANONICAL_BGF_DATA);
        let result = BgfFile::read_from(&mut reader);
        assert!(
            result.is_ok(),
            "Parsing a canonical BGF file should succeed"
        );

        let (system, metadata) = result.unwrap();

        assert_eq!(metadata.header_lines.len(), 2);
        assert_eq!(
            metadata.header_lines[1],
            "REMARK A standard, clean BGF file"
        );

        assert_eq!(system.atoms_iter().count(), 7);
        assert_eq!(system.chains_iter().count(), 2);
        assert_eq!(system.residues_iter().count(), 2);
        assert_eq!(system.bonds().len(), 6);

        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let gly_res_id = system.find_residue_by_id(chain_a_id, 1).unwrap();
        let gly = system.residue(gly_res_id).unwrap();

        assert_eq!(gly.name, "GLY");
        assert_eq!(gly.residue_type, Some(ResidueType::Glycine));
        assert_eq!(gly.atoms().len(), 4);

        let n_atom_id = gly.get_first_atom_id_by_name("N").unwrap();
        let ca_atom_id = gly.get_first_atom_id_by_name("CA").unwrap();

        assert!(
            system
                .bonds()
                .iter()
                .any(|b| b.contains(n_atom_id) && b.contains(ca_atom_id))
        );
    }

    #[test]
    fn read_from_handles_disordered_file_and_discards_serials() {
        let mut reader = Cursor::new(DISORDERED_BGF_DATA);
        let result = BgfFile::read_from(&mut reader);
        assert!(
            result.is_ok(),
            "Parsing a disordered BGF file should succeed"
        );

        let (system, _) = result.unwrap();

        assert_eq!(system.atoms_iter().count(), 5);
        assert_eq!(system.chains_iter().count(), 2);
        assert_eq!(system.residues_iter().count(), 3);

        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let ala_res_id = system.find_residue_by_id(chain_a_id, 2).unwrap();
        let ala = system.residue(ala_res_id).unwrap();

        let ala_n_id = ala.get_first_atom_id_by_name("N").unwrap();
        let ala_ca_id = ala.get_first_atom_id_by_name("CA").unwrap();

        assert_eq!(system.bonds().len(), 2);
        assert!(
            system
                .bonds()
                .iter()
                .any(|b| b.contains(ala_n_id) && b.contains(ala_ca_id))
        );

        let any_atom = system.atom(ala_n_id).unwrap();

        assert!(any_atom.name == "N");
    }

    #[test]
    fn read_from_fails_on_conect_with_unknown_serial() {
        let bad_data = "ATOM       1 N     GLY A     1       0.0       0.0       0.0 N_R    1 3 -0.35000\nCONECT     1    99\nEND\n";
        let mut reader = Cursor::new(bad_data);
        let result = BgfFile::read_from(&mut reader);

        assert!(matches!(result, Err(BgfError::Parse { line_num: 2, .. })));
    }

    #[test]
    fn write_to_produces_canonically_sorted_and_numbered_output() {
        let mut reader = Cursor::new(DISORDERED_BGF_DATA);
        let (system, metadata) = BgfFile::read_from(&mut reader).unwrap();

        let mut buffer = Vec::new();
        let result = BgfFile::write_to(&system, &metadata, &mut buffer);
        assert!(result.is_ok());

        let output = String::from_utf8(buffer).unwrap();
        let output_lines: Vec<_> = output.lines().collect();

        assert_eq!(output_lines[0], "BIOGRF 332");
        assert_eq!(
            output_lines[1],
            "REMARK A messy BGF file with non-sequential serials and disordered records."
        );

        let atom_lines: Vec<_> = output_lines
            .iter()
            .filter(|l| l.starts_with("ATOM") || l.starts_with("HETATM"))
            .collect();
        assert_eq!(atom_lines.len(), 5);

        for (i, line) in atom_lines.iter().enumerate() {
            let serial_str = line.get(7..12).unwrap().trim();
            assert_eq!(
                serial_str.parse::<usize>().unwrap(),
                i + 1,
                "Serials should be sequential"
            );
        }

        assert!(atom_lines[0].contains(" N     GLY A     1"));
        assert!(atom_lines[1].contains(" C     GLY A     1"));
        assert!(atom_lines[2].contains(" N     ALA A     2"));
        assert!(atom_lines[3].contains(" CA    ALA A     2"));
        assert!(atom_lines[4].contains(" O     HOH W     3"));

        let conect_lines: Vec<_> = output_lines
            .iter()
            .filter(|l| l.starts_with("CONECT"))
            .collect();
        assert_eq!(conect_lines.len(), 2);

        let bond1_found = conect_lines
            .iter()
            .any(|l| l.contains("    1") && l.contains("    2"));
        assert!(
            bond1_found,
            "CONECT record for bond between new serials 1 and 2 not found"
        );

        let bond2_found = conect_lines
            .iter()
            .any(|l| l.contains("    3") && l.contains("    4"));
        assert!(
            bond2_found,
            "CONECT record for bond between new serials 3 and 4 not found"
        );
    }

    #[test]
    fn write_system_to_produces_correct_default_header() {
        let mut system = MolecularSystem::new();
        system.add_chain('A', ChainType::Protein);

        let mut buffer = Vec::new();
        BgfFile::write_system_to(&system, &mut buffer).unwrap();
        let output = String::from_utf8(buffer).unwrap();

        assert!(output.contains("BIOGRF  332"));
        assert!(output.contains("FORCEFIELD DREIDING"));
    }

    #[test]
    fn write_to_includes_format_lines() {
        let mut system = MolecularSystem::new();
        system.add_chain('A', ChainType::Protein);

        let mut buffer = Vec::new();
        BgfFile::write_system_to(&system, &mut buffer).unwrap();
        let output = String::from_utf8(buffer).unwrap();

        assert!(output.contains(
            "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)"
        ));
        assert!(output.contains("FORMAT CONECT (a6,12i6)"));
    }
}
