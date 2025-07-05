use super::traits::MolecularFile;
use crate::core::models::atom::Atom;
use crate::core::models::chain::ChainType;
use crate::core::models::ids::{ChainId, ResidueId};
use crate::core::models::residue::ResidueType;
use crate::core::models::system::MolecularSystem;
use crate::core::models::topology::BondOrder;
use nalgebra::Point3;
use std::collections::HashMap;
use std::io::{self, BufRead, Write};
use thiserror::Error;

#[derive(Debug, Default, Clone)]
pub struct BgfMetadata {
    pub header_lines: Vec<String>, // Lines before the ATOM records
    pub order_lines: Vec<String>,  // Lines containing bond orders
    pub footer_lines: Vec<String>, // Lines after the ATOM records
}

#[derive(Debug, Default)]
pub struct BgfFile;

#[derive(Debug, Error)]
pub enum BgfError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    #[error("Parse error on line {line_num}: {message}")]
    Parse { line_num: usize, message: String },
    #[error("Logic error: {0}")]
    Logic(String),
}

impl MolecularFile for BgfFile {
    type Metadata = BgfMetadata;
    type Error = BgfError;

    fn read_from(
        reader: &mut impl BufRead,
    ) -> Result<(MolecularSystem, Self::Metadata), Self::Error> {
        let mut system = MolecularSystem::new();
        let mut metadata = BgfMetadata::default();

        let mut atom_lines = Vec::new();
        let mut conect_lines = Vec::new();
        let mut atom_section_started = false;

        for (i, line_result) in reader.lines().enumerate() {
            let line = line_result?;
            if line.trim().is_empty() {
                continue;
            }
            if line.trim() == "END" {
                break;
            }

            let record = line.get(0..6).unwrap_or("").trim();
            match record {
                "ATOM" | "HETATM" => {
                    atom_section_started = true;
                    atom_lines.push((i + 1, line));
                }
                "CONECT" => conect_lines.push((i + 1, line)),
                "ORDER" => metadata.order_lines.push(line),
                "FORMAT" => {}
                _ => {
                    if !atom_section_started {
                        metadata.header_lines.push(line);
                    } else {
                        metadata.footer_lines.push(line);
                    }
                }
            }
        }

        let mut current_chain_id: Option<ChainId> = None;
        let mut current_residue_id: Option<ResidueId> = None;

        for (line_num, line) in atom_lines {
            let atom_info = parse_atom_line(&line).map_err(|msg| BgfError::Parse {
                line_num,
                message: msg,
            })?;

            let chain_char = atom_info.chain_char;
            let mut new_chain = false;
            if let Some(id) = current_chain_id {
                if system.chain(id).unwrap().id != chain_char {
                    new_chain = true;
                }
            } else {
                new_chain = true;
            }

            if new_chain {
                let chain_type = match atom_info.record_type.as_str() {
                    "ATOM" => ChainType::Protein,
                    "HETATM" => match atom_info.res_name.as_str() {
                        "HOH" | "TIP3" => ChainType::Water,
                        _ => ChainType::Ligand,
                    },
                    _ => ChainType::Other,
                };
                current_chain_id = Some(system.add_chain(chain_char, chain_type));
                current_residue_id = None;
            }
            let active_chain_id = current_chain_id.unwrap();

            let res_seq = atom_info.res_seq;
            let mut new_residue = false;
            if let Some(id) = current_residue_id {
                if system.residue(id).unwrap().id != res_seq {
                    new_residue = true;
                }
            } else {
                new_residue = true;
            }

            if new_residue {
                current_residue_id = system.add_residue(
                    active_chain_id,
                    res_seq,
                    &atom_info.res_name,
                    ResidueType::from_str_optional(&atom_info.res_name),
                );
            }
            let active_residue_id = current_residue_id
                .ok_or_else(|| BgfError::Logic("Failed to create or find residue".to_string()))?;

            let mut atom = Atom::new(
                atom_info.serial,
                &atom_info.name,
                active_residue_id,
                atom_info.pos,
            );
            atom.partial_charge = atom_info.charge;
            atom.force_field_type = atom_info.ff_type;

            system
                .add_atom_to_residue(active_residue_id, atom)
                .ok_or_else(|| BgfError::Logic("Failed to add atom".to_string()))?;
        }

        for (line_num, line) in conect_lines {
            let (base_serial, connected_serials) =
                parse_conect_line(&line).map_err(|msg| BgfError::Parse {
                    line_num,
                    message: msg,
                })?;

            let base_id =
                system
                    .find_atom_by_serial(base_serial)
                    .ok_or_else(|| BgfError::Parse {
                        line_num,
                        message: format!("Atom serial {} not found for CONECT record", base_serial),
                    })?;

            for &connected_serial in &connected_serials {
                let connected_id =
                    system
                        .find_atom_by_serial(connected_serial)
                        .ok_or_else(|| BgfError::Parse {
                            line_num,
                            message: format!(
                                "Atom serial {} not found for CONECT record",
                                connected_serial
                            ),
                        })?;
                system.add_bond(base_id, connected_id, BondOrder::Single);
            }
        }

        Ok((system, metadata))
    }

    fn write_to(
        system: &MolecularSystem,
        metadata: &Self::Metadata,
        writer: &mut impl Write,
    ) -> Result<(), Self::Error> {
        for line in &metadata.header_lines {
            writeln!(writer, "{}", line)?;
        }

        for (_chain_id, chain) in system.chains_iter() {
            let record_type = match chain.chain_type {
                ChainType::Protein | ChainType::DNA | ChainType::RNA => "ATOM  ",
                _ => "HETATM",
            };
            for &residue_id in chain.residues() {
                let residue = system.residue(residue_id).unwrap();
                for &atom_id in residue.atoms() {
                    let atom = system.atom(atom_id).unwrap();
                    let atoms_connected = 0; // TODO: Replace with actual logic (placeholder)
                    let lone_pair = 0; // TODO: Replace with actual logic (placeholder)
                    let atom_line = format_atom_line(
                        record_type,
                        atom,
                        residue,
                        chain,
                        atoms_connected,
                        lone_pair,
                    );
                    writeln!(writer, "{}", atom_line)?;
                }
            }
        }

        for line in &metadata.footer_lines {
            writeln!(writer, "{}", line)?;
        }

        let mut bond_serials_written: HashMap<usize, Vec<usize>> = HashMap::new();
        for bond in system.bonds() {
            let atom1 = system.atom(bond.atom1_id).unwrap();
            let atom2 = system.atom(bond.atom2_id).unwrap();

            if atom1.serial < atom2.serial {
                bond_serials_written
                    .entry(atom1.serial)
                    .or_default()
                    .push(atom2.serial);
            } else if atom2.serial < atom1.serial {
                bond_serials_written
                    .entry(atom2.serial)
                    .or_default()
                    .push(atom1.serial);
            }
        }

        let mut sorted_base_serials: Vec<_> = bond_serials_written.keys().copied().collect();
        sorted_base_serials.sort_unstable();

        for base_serial in sorted_base_serials {
            let neighbors = bond_serials_written.get_mut(&base_serial).unwrap();
            neighbors.sort_unstable();
            for chunk in neighbors.chunks(14) {
                let mut connect_line = format!("CONECT {:>5}", base_serial);
                for &neighbor_serial in chunk {
                    connect_line.push_str(&format!("{:>6}", neighbor_serial));
                }
                writeln!(writer, "{}", connect_line)?;
            }
        }

        for line in &metadata.order_lines {
            writeln!(writer, "{}", line)?;
        }

        writeln!(writer, "END")?;
        Ok(())
    }

    fn write_system_to(
        system: &MolecularSystem,
        writer: &mut impl Write,
    ) -> Result<(), Self::Error> {
        let minimal_metadata = BgfMetadata {
            header_lines: vec![
                "BIOGRF  332".to_string(),
                "FORCEFIELD DREIDING".to_string(),
                "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)"
                    .to_string(),
            ],
            ..Default::default()
        };
        Self::write_to(system, &minimal_metadata, writer)
    }
}

#[derive(Debug, Clone)]
struct ParsedAtomInfo {
    record_type: String,
    serial: usize,
    name: String,
    res_name: String,
    chain_char: char,
    res_seq: isize,
    pos: Point3<f64>,
    charge: f64,
    ff_type: String,
}

fn parse_atom_line(line: &str) -> Result<ParsedAtomInfo, String> {
    let get_slice = |start: usize, end: usize| -> Result<&str, String> {
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
        .map_err(|e| format!("Invalid residue ID: {}", e))?;
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

fn parse_conect_line(line: &str) -> Result<(usize, Vec<usize>), String> {
    let mut parts = line.split_whitespace().skip(1).map(str::parse);
    let base_serial = parts
        .next()
        .ok_or("Missing base serial in CONECT")?
        .map_err(|e| format!("Invalid base serial: {}", e))?;
    let connected_serials = parts
        .collect::<Result<Vec<usize>, _>>()
        .map_err(|e| format!("Invalid connected serial: {}", e))?;
    Ok((base_serial, connected_serials))
}

fn format_atom_line(
    record_type: &str,
    atom: &Atom,
    residue: &crate::core::models::residue::Residue,
    chain: &crate::core::models::chain::Chain,
    atoms_connected: u8,
    lone_pair: u8,
) -> String {
    format!(
        "{:<6} {:>5} {:<5} {:>3} {:1} {:>5}{:>10.5}{:>10.5}{:>10.5} {:<5}{:>3}{:>2} {:>8.5}",
        record_type,
        atom.serial,
        atom.name,
        residue.name,
        chain.id,
        residue.id,
        atom.position.x,
        atom.position.y,
        atom.position.z,
        atom.force_field_type,
        atoms_connected,
        lone_pair,
        atom.partial_charge
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::models::chain::Chain;
    use crate::core::models::ids::{ChainId, ResidueId};
    use crate::core::models::residue::Residue;
    use std::io::Cursor;

    fn create_test_system() -> MolecularSystem {
        let mut system = MolecularSystem::new();
        let chain_a = system.add_chain('A', ChainType::Protein);
        let res1 = system
            .add_residue(chain_a, 1, "GLY", ResidueType::from_str_optional("GLY"))
            .unwrap();
        let mut atom_n = Atom::new(1, "N", res1, Point3::new(-0.416, -0.535, 0.0));
        atom_n.partial_charge = -0.35;
        atom_n.force_field_type = "N_3".to_string();
        let mut atom_h = Atom::new(2, "H", res1, Point3::new(0.564, -0.535, 0.0));
        atom_h.partial_charge = 0.27;
        atom_h.force_field_type = "H_".to_string();

        let id_n = system.add_atom_to_residue(res1, atom_n).unwrap();
        let id_h = system.add_atom_to_residue(res1, atom_h).unwrap();
        system.add_bond(id_n, id_h, BondOrder::Single);
        system
    }

    const TEST_BGF_DATA: &str = r#"BIOGRF 332
FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)
ATOM       1 N     GLY A     1  -0.41600  -0.53500   0.00000 N_3    1 3 -0.35000
ATOM       2 H     GLY A     1   0.56400  -0.53500   0.00000 H_     1 1  0.27000
CONECT     1     2
END
"#;

    #[test]
    fn read_from_parses_valid_bgf_file() {
        let mut reader = Cursor::new(TEST_BGF_DATA);
        let result = BgfFile::read_from(&mut reader);
        assert!(result.is_ok());
        let (system, metadata) = result.unwrap();

        assert_eq!(system.atoms_iter().count(), 2);
        assert_eq!(system.chains_iter().count(), 1);
        assert_eq!(system.residues_iter().count(), 1);
        assert_eq!(system.bonds().len(), 1);
        assert_eq!(metadata.header_lines.len(), 1);
        assert_eq!(metadata.header_lines[0], "BIOGRF 332");
        let atom_n_id = system.find_atom_by_serial(1).unwrap();
        let atom_h_id = system.find_atom_by_serial(2).unwrap();
        assert_eq!(system.atom(atom_n_id).unwrap().name, "N");
        assert_eq!(system.atom(atom_h_id).unwrap().partial_charge, 0.27);
    }

    #[test]
    fn read_from_handles_hetatm_and_water_chain() {
        let bgf_data = "HETATM     1 O     HOH A     1   0.00000   0.00000   0.00000 O_3    1 2 -0.83400\nEND\n";
        let mut reader = Cursor::new(bgf_data);
        let (system, _) = BgfFile::read_from(&mut reader).unwrap();
        let chain_id = system.find_chain_by_id('A').unwrap();
        assert_eq!(system.chain(chain_id).unwrap().chain_type, ChainType::Water);
    }

    #[test]
    fn read_from_handles_multiple_chains_and_residues() {
        let bgf_data = r#"
ATOM       1 N     GLY A     1  -0.41600  -0.53500   0.00000 N_3    1 3 -0.35000
ATOM       2 CA    ALA A     2   1.00000   1.00000   1.00000 C_3    1 4  0.10000
ATOM       3 N     SER B     1   2.00000   2.00000   2.00000 N_3    1 3 -0.35000
END
"#;
        let mut reader = Cursor::new(bgf_data);
        let (system, _) = BgfFile::read_from(&mut reader).unwrap();
        assert_eq!(system.chains_iter().count(), 2);
        assert_eq!(system.atoms_iter().count(), 3);

        let chain_a_id = system.find_chain_by_id('A').unwrap();
        let chain_b_id = system.find_chain_by_id('B').unwrap();
        assert_eq!(system.chain(chain_a_id).unwrap().residues().len(), 2);
        assert_eq!(system.chain(chain_b_id).unwrap().residues().len(), 1);
    }

    #[test]
    fn read_from_returns_error_for_malformed_atom_line() {
        let bgf_data = "ATOM      X  N    GLY A    1      -0.41600  -0.53500   0.00000 N_3      1 3     -0.35000\n";
        let mut reader = Cursor::new(bgf_data);
        let result = BgfFile::read_from(&mut reader);
        assert!(matches!(result, Err(BgfError::Parse { line_num: 1, .. })));
    }

    #[test]
    fn read_from_returns_error_for_malformed_conect_line() {
        let bgf_data_1 = "ATOM       1 N     GLY A     1       0.0       0.0       0.0 N_3    0 0      0.0\nCONECT     1     X\nEND\n";
        let mut reader_1 = Cursor::new(bgf_data_1);
        let result_1 = BgfFile::read_from(&mut reader_1);
        assert!(matches!(result_1, Err(BgfError::Parse { line_num: 2, .. })));

        let bgf_data_2 = "ATOM       1 N     GLY A     1       0.0       0.0       0.0 N_3    0 0      0.0\nCONECT     1    99\nEND\n";
        let mut reader_2 = Cursor::new(bgf_data_2);
        let result_2 = BgfFile::read_from(&mut reader_2);
        assert!(
            matches!(result_2, Err(BgfError::Parse { line_num: 2, message: m, .. }) if m.contains("Atom serial 99 not found"))
        );
    }

    #[test]
    fn write_to_produces_correct_bgf_format() {
        let system = create_test_system();
        let metadata = BgfMetadata {
            header_lines: vec!["BIOGRF 332".to_string()],
            ..Default::default()
        };
        let mut buffer = Vec::new();
        let result = BgfFile::write_to(&system, &metadata, &mut buffer);
        assert!(result.is_ok());

        let output = String::from_utf8(buffer).unwrap();
        assert!(output.starts_with("BIOGRF 332"));
        assert!(output.contains("ATOM       1 N     GLY A     1"));
        assert!(output.contains("ATOM       2 H     GLY A     1"));
        assert!(output.contains("CONECT     1     2"));
        assert!(output.trim_end().ends_with("END"));
    }

    #[test]
    fn write_system_to_uses_default_metadata() {
        let system = create_test_system();
        let mut buffer = Vec::new();
        let result = BgfFile::write_system_to(&system, &mut buffer);
        assert!(result.is_ok());

        let output = String::from_utf8(buffer).unwrap();
        assert!(output.contains("BIOGRF  332"));
        assert!(output.contains("FORCEFIELD DREIDING"));
    }

    #[test]
    fn read_from_handles_duplicate_conect_records() {
        let bgf_data = r#"ATOM       1 N     GLY A     1       0.0       0.0       0.0 N_3    0 0      0.0
ATOM       2 C     GLY A     1       1.0       0.0       0.0 C_3    0 0      0.0
CONECT     1     2
CONECT     2     1
END
"#;
        let mut reader = Cursor::new(bgf_data);
        let (system, _) = BgfFile::read_from(&mut reader).unwrap();
        assert_eq!(
            system.bonds().len(),
            1,
            "Should only create one bond for duplicate CONECT records"
        );
    }

    #[test]
    fn parse_atom_line_succeeds_on_valid_input() {
        let line =
            "ATOM       1 N     GLY A     1  -0.41600  -0.53500   0.00000 N_3    1 3 -0.35000";
        let result = parse_atom_line(line);
        assert!(result.is_ok());
        let info = result.unwrap();
        assert_eq!(info.serial, 1);
        assert_eq!(info.name, "N");
        assert_eq!(info.charge, -0.35);
        assert_eq!(info.ff_type, "N_3");
    }

    #[test]
    fn parse_atom_line_fails_on_invalid_number() {
        let line = "ATOM      X  N    GLY A    1      -0.41600  -0.53500   0.00000 N_3      1 3     -0.35000";
        let result = parse_atom_line(line);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "Invalid serial: invalid digit found in string".to_string()
        );
    }

    #[test]
    fn parse_conect_line_succeeds_on_valid_input() {
        let line = "CONECT     1     2     3     4";
        let result = parse_conect_line(line);
        assert!(result.is_ok());
        let (base, connected) = result.unwrap();
        assert_eq!(base, 1);
        assert_eq!(connected, vec![2, 3, 4]);
    }

    #[test]
    fn parse_conect_line_fails_on_invalid_format() {
        let line = "CONECT     1     X";
        let result = parse_conect_line(line);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "Invalid connected serial: invalid digit found in string".to_string()
        );
    }

    #[test]
    fn format_atom_line_produces_correct_string() {
        let chain_id = ChainId::default();
        let res_id = ResidueId::default();
        let atom = Atom::new(1, "CA", res_id, Point3::new(1.1, 2.2, 3.3));
        let residue = Residue::new(5, "ALA", Some(ResidueType::Alanine), chain_id);
        let chain = Chain::new('A', ChainType::Protein);

        let mut atom_to_format = atom;
        atom_to_format.force_field_type = "C_3".to_string();
        atom_to_format.partial_charge = 0.123;

        let line = format_atom_line("ATOM  ", &atom_to_format, &residue, &chain, 4, 0);
        let expected =
            "ATOM       1 CA    ALA A     5   1.10000   2.20000   3.30000 C_3    4 0  0.12300";
        assert_eq!(line, expected);
    }
}
