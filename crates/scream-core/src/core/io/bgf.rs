use super::traits::MolecularFile;
use crate::core::models::builder::MolecularSystemBuilder;
use crate::core::models::chain::ChainType;
use crate::core::models::system::MolecularSystem;
use crate::core::models::topology::BondOrder;
use nalgebra::Point3;
use std::collections::HashMap;
use std::io::{self, BufRead, Write};
use thiserror::Error;

#[derive(Debug, Default, Clone)]
pub struct BgfMetadata {
    pub extra_lines: Vec<String>,
}

#[derive(Debug, Default)]
pub struct BgfFile;

#[derive(Debug, Error)]
pub enum BgfError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    #[error("Parse error on line {line_num}: {message}")]
    Parse { line_num: usize, message: String },
}

impl MolecularFile for BgfFile {
    type Metadata = BgfMetadata;
    type Error = BgfError;

    fn read_from(
        reader: &mut impl BufRead,
    ) -> Result<(MolecularSystem, Self::Metadata), Self::Error> {
        let mut builder = MolecularSystemBuilder::new();
        let mut metadata = BgfMetadata::default();
        let mut connectivity_lines: Vec<(String, usize)> = Vec::new();

        let mut current_chain_id = '\0'; // Use a null character as an uninitialized sentinel
        let mut current_res_id = isize::MIN;

        for (line_num, line_result) in reader.lines().enumerate() {
            let line = line_result?;
            let line_num = line_num + 1; // 1-based for error reporting

            if line.starts_with("ATOM") || line.starts_with("HETATM") {
                let (record_type, serial, name, res_name, chain_id, res_id, pos, charge, ff_type) =
                    parse_atom_line(&line).map_err(|msg| BgfError::Parse {
                        line_num,
                        message: msg,
                    })?;

                if chain_id != current_chain_id {
                    let chain_type = match record_type.as_str() {
                        "ATOM" => ChainType::Protein,
                        "HETATM" => match res_name.as_str() {
                            "HOH" | "TIP3" => ChainType::Water,
                            _ => ChainType::Other, // Default for other HETATMs
                        },
                        _ => ChainType::Other,
                    };
                    builder.start_chain(chain_id, chain_type);
                    current_chain_id = chain_id;
                    current_res_id = isize::MIN; // Reset residue on new chain
                }

                if res_id != current_res_id {
                    builder.start_residue(res_id, &res_name);
                    current_res_id = res_id;
                }

                builder.add_atom(serial, &name, &res_name, pos, Some(charge), Some(&ff_type));
            } else if line.starts_with("CONECT") || line.starts_with("ORDER") {
                connectivity_lines.push((line, line_num));
            } else if line.trim().starts_with("FORMAT") {
                // Ignore format lines
            } else if line.trim() == "END" {
                break;
            } else if !line.trim().is_empty() {
                metadata.extra_lines.push(line);
            }
        }

        for (conn_line, line_num) in connectivity_lines {
            if conn_line.starts_with("CONECT") {
                let (base_serial, connected_serials) =
                    parse_conect_line(&conn_line).map_err(|msg| BgfError::Parse {
                        line_num,
                        message: msg,
                    })?;
                for &connected_serial in &connected_serials {
                    builder.add_bond(base_serial, connected_serial, BondOrder::Single);
                }
            } else {
                metadata.extra_lines.push(conn_line);
            }
        }

        Ok((builder.build(), metadata))
    }

    fn write_to(
        system: &MolecularSystem,
        metadata: &Self::Metadata,
        writer: &mut impl Write,
    ) -> Result<(), Self::Error> {
        for line in &metadata.extra_lines {
            writeln!(writer, "{}", line)?;
        }

        for chain in system.chains() {
            let record_type = match chain.chain_type {
                ChainType::Protein | ChainType::DNA | ChainType::RNA => "ATOM  ",
                _ => "HETATM",
            };
            for residue in chain.residues() {
                for &atom_index in &residue.atom_indices {
                    let atom = system.get_atom(atom_index).unwrap();

                    // TODO: Replace with actual logic once implemented in the core.
                    let atoms_connected = 0;
                    let lone_pair = 0;

                    let atom_line = format_atom_line(record_type, atom, atoms_connected, lone_pair);
                    writeln!(writer, "{}", atom_line)?;
                }
            }
        }

        let mut bond_map: HashMap<usize, Vec<usize>> = HashMap::new();
        for bond in system.bonds() {
            bond_map
                .entry(bond.atom1_idx)
                .or_default()
                .push(bond.atom2_idx);
            bond_map
                .entry(bond.atom2_idx)
                .or_default()
                .push(bond.atom1_idx);
        }

        let mut sorted_indices: Vec<_> = bond_map.keys().copied().collect();
        sorted_indices.sort();

        for &atom_index in &sorted_indices {
            let base_serial = system.get_atom(atom_index).unwrap().serial;
            let mut connect_line = format!("CONECT {:>5}", base_serial);

            let mut neighbor_serials: Vec<_> = bond_map[&atom_index]
                .iter()
                .map(|&idx| system.get_atom(idx).unwrap().serial)
                .collect();
            neighbor_serials.sort(); // Sort neighbors for deterministic output.

            for serial in neighbor_serials {
                connect_line.push_str(&format!(" {:>5}", serial));
            }
            writeln!(writer, "{}", connect_line)?;
        }

        writeln!(writer, "END")?;
        Ok(())
    }

    fn write_system_to(
        system: &MolecularSystem,
        writer: &mut impl Write,
    ) -> Result<(), Self::Error> {
        let minimal_metadata = BgfMetadata {
            extra_lines: vec![
                "BIOGRF  332".to_string(),
                "FORCEFIELD DREIDING".to_string(),
                "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)"
                    .to_string(),
                "FORMAT CONECT (a6,14i6)".to_string(),
            ],
        };
        Self::write_to(system, &minimal_metadata, writer)
    }
}

fn parse_atom_line(
    line: &str,
) -> Result<
    (
        String,
        usize,
        String,
        String,
        char,
        isize,
        Point3<f64>,
        f64,
        String,
    ),
    String,
> {
    let record_type = line.get(0..6).unwrap_or("").trim().to_string();
    let serial = line
        .get(7..12)
        .and_then(|s| s.trim().parse().ok())
        .ok_or("Invalid serial number".to_string())?;
    let name = line.get(13..18).unwrap_or("").trim().to_string();
    let res_name = line.get(19..22).unwrap_or("").trim().to_string();
    let chain_id = line
        .get(23..24)
        .and_then(|s| s.chars().next())
        .unwrap_or(' ');
    let res_id = line
        .get(25..30)
        .and_then(|s| s.trim().parse().ok())
        .ok_or("Invalid residue ID".to_string())?;
    let x = line
        .get(30..40)
        .and_then(|s| s.trim().parse().ok())
        .ok_or("Invalid X coordinate".to_string())?;
    let y = line
        .get(40..50)
        .and_then(|s| s.trim().parse().ok())
        .ok_or("Invalid Y coordinate".to_string())?;
    let z = line
        .get(50..60)
        .and_then(|s| s.trim().parse().ok())
        .ok_or("Invalid Z coordinate".to_string())?;
    let ff_type = line.get(61..66).unwrap_or("").trim().to_string();
    let charge = line
        .get(72..80)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0.0);

    Ok((
        record_type,
        serial,
        name,
        res_name,
        chain_id,
        res_id,
        Point3::new(x, y, z),
        charge,
        ff_type,
    ))
}

fn parse_conect_line(line: &str) -> Result<(usize, Vec<usize>), String> {
    let parts: Vec<_> = line.split_whitespace().collect();
    if parts.len() < 2 || parts[0] != "CONECT" {
        return Err("Invalid CONECT line format".to_string());
    }
    let base_serial = parts[1]
        .parse()
        .map_err(|e| format!("Invalid base serial in CONECT: {}", e))?;
    let connected_serials = parts[2..].iter().filter_map(|s| s.parse().ok()).collect();
    Ok((base_serial, connected_serials))
}

fn format_atom_line(
    record_type: &str,
    atom: &crate::core::models::atom::Atom,
    atoms_connected: u8,
    lone_pair: u8,
) -> String {
    format!(
        "{:<6} {:>5} {:<5} {:>3} {:1} {:>5}{:>10.5}{:>10.5}{:>10.5} {:<5}{:>3}{:>2} {:>8.5}",
        record_type,
        atom.serial,
        atom.name,
        atom.res_name,
        atom.chain_id,
        atom.res_id,
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
    use crate::core::models::atom::{Atom, AtomFlags};
    use std::io::Cursor;

    fn create_test_system() -> MolecularSystem {
        let mut builder = MolecularSystemBuilder::new();
        builder
            .start_chain('A', ChainType::Protein)
            .start_residue(1, "GLY")
            .add_atom(
                1,
                "N",
                "GLY",
                Point3::new(-0.416, -0.535, 0.0),
                Some(-0.35),
                Some("N_3"),
            )
            .add_atom(
                2,
                "H",
                "GLY",
                Point3::new(0.564, -0.535, 0.0),
                Some(0.27),
                Some("H_"),
            )
            .add_bond(1, 2, BondOrder::Single);
        builder.build()
    }

    #[test]
    fn read_from_parses_valid_bgf_file() {
        let bgf_data = r#"BIOGRF 332
FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)
ATOM       1 N     GLY A     1  -0.41600  -0.53500   0.00000 N_3    1 3 -0.35000
ATOM       2 H     GLY A     1   0.56400  -0.53500   0.00000 H_     1 1  0.27000
CONECT     1     2
END
"#;
        let mut reader = Cursor::new(bgf_data);
        let result = BgfFile::read_from(&mut reader);
        assert!(result.is_ok());
        let (system, metadata) = result.unwrap();

        assert_eq!(system.atoms().len(), 2);
        assert_eq!(system.chains().len(), 1);
        assert_eq!(system.chains()[0].residues().len(), 1);
        assert_eq!(system.bonds().len(), 1);
        assert_eq!(metadata.extra_lines, vec!["BIOGRF 332"]);
        assert_eq!(system.get_atom_by_serial(1).unwrap().name, "N");
        assert_eq!(system.get_atom_by_serial(2).unwrap().partial_charge, 0.27);
    }

    #[test]
    fn read_from_handles_hetatm_and_water_chain() {
        let bgf_data = "HETATM     1 O     HOH A     1   0.00000   0.00000   0.00000 O_3    1 2 -0.83400\nEND\n";
        let mut reader = Cursor::new(bgf_data);
        let (system, _) = BgfFile::read_from(&mut reader).unwrap();
        assert_eq!(system.chains().len(), 1);
        assert_eq!(system.chains()[0].chain_type, ChainType::Water);
    }

    #[test]
    fn read_from_handles_multiple_chains_and_residues() {
        let bgf_data = r#"
ATOM       1 N     GLY A     1  -0.41600  -0.53500   0.00000 N_3    1 3 -0.35000
ATOM       2 CA    ALA A    2    1.00000   1.00000   1.00000 C_3    1 4  0.10000
ATOM       3 N     SER B    1    2.00000   2.00000   2.00000 N_3    1 3 -0.35000
END
"#;
        let mut reader = Cursor::new(bgf_data);
        let (system, _) = BgfFile::read_from(&mut reader).unwrap();
        assert_eq!(system.chains().len(), 2);
        assert_eq!(system.chains()[0].residues().len(), 2);
        assert_eq!(system.chains()[1].residues().len(), 1);
        assert_eq!(system.atoms().len(), 3);
    }

    #[test]
    fn read_from_returns_error_for_malformed_atom_line() {
        let bgf_data = "ATOM      X  N    GLY A    1      -0.41600  -0.53500   0.00000 N_3      1 3     -0.35000\n";
        let mut reader = Cursor::new(bgf_data);
        let result = BgfFile::read_from(&mut reader);
        assert!(matches!(result, Err(BgfError::Parse { line_num: 1, .. })));
    }

    #[test]
    fn read_from_ignores_malformed_conect_line_parts() {
        let bgf_data = "ATOM       1  N    GLY A    1   -0.41600  -0.53500   0.00000 N_3    1 3 -0.35000\nCONECT     1     X\nEND\n";
        let mut reader = Cursor::new(bgf_data);
        let (system, _) = BgfFile::read_from(&mut reader).unwrap();
        assert!(system.bonds().is_empty());
    }

    #[test]
    fn write_to_produces_correct_bgf_format() {
        let system = create_test_system();
        let metadata = BgfMetadata {
            extra_lines: vec!["BIOGRF 332".to_string()],
        };
        let mut buffer = Vec::new();
        let result = BgfFile::write_to(&system, &metadata, &mut buffer);
        assert!(result.is_ok());

        let output = String::from_utf8(buffer).unwrap();
        assert!(output.starts_with("BIOGRF 332"));
        assert!(output.contains(
            "ATOM       1 N     GLY A     1  -0.41600  -0.53500   0.00000 N_3    0 0 -0.35000"
        ));
        assert!(output.contains("CONECT     1     2"));
        assert!(output.ends_with("END\n"));
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
    fn parse_atom_line_succeeds_on_valid_input() {
        let line =
            "ATOM       1 N     GLY A     1  -0.41600  -0.53500   0.00000 N_3    1 3 -0.35000";
        let result = parse_atom_line(line);
        assert!(result.is_ok());
        let (_, serial, name, _, _, _, _, charge, ff_type) = result.unwrap();
        assert_eq!(serial, 1);
        assert_eq!(name, "N");
        assert_eq!(charge, -0.35);
        assert_eq!(ff_type, "N_3");
    }

    #[test]
    fn parse_atom_line_fails_on_invalid_number() {
        let line = "ATOM      X  N    GLY A    1      -0.41600  -0.53500   0.00000 N_3      1 3     -0.35000";
        let result = parse_atom_line(line);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Invalid serial number".to_string());
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
        let line = "CONNECT     1     2";
        let result = parse_conect_line(line);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "Invalid CONECT line format".to_string()
        );
    }

    #[test]
    fn format_atom_line_produces_correct_string() {
        let atom = Atom {
            index: 1,
            serial: 1,
            name: "CA".to_string(),
            res_name: "ALA".to_string(),
            res_id: 5,
            chain_id: 'A',
            force_field_type: "C_3".to_string(),
            partial_charge: 0.123,
            position: Point3::new(1.1, 2.2, 3.3),
            flags: AtomFlags::empty(),
            delta: 0.0,
            vdw_radius: 1.0,
            vdw_well_depth: 0.1,
            hbond_type_id: 0,
        };
        let line = format_atom_line("ATOM  ", &atom, 4, 0);
        let expected =
            "ATOM       1 CA    ALA A     5   1.10000   2.20000   3.30000 C_3    4 0  0.12300";
        assert_eq!(line, expected);
    }
}
