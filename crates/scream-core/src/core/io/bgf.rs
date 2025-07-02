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
        let mut connectivity_lines = Vec::new();

        let mut current_chain_id = '\0';
        let mut current_res_id = isize::MIN;

        for (line_num, line_result) in reader.lines().enumerate() {
            let line = line_result?;
            let line_num = line_num + 1; // 1-based for error reporting

            if line.starts_with("ATOM") || line.starts_with("HETATM") {
                let (serial, name, res_name, chain_id, res_id, pos, charge, ff_type) =
                    parse_atom_line(&line).map_err(|msg| BgfError::Parse {
                        line_num,
                        message: msg,
                    })?;

                if chain_id != current_chain_id {
                    builder.start_chain(chain_id, ChainType::Protein); // Assuming protein for now
                    current_chain_id = chain_id;
                    current_res_id = isize::MIN; // Reset residue ID when chain changes
                }
                if res_id != current_res_id {
                    builder.start_residue(res_id, &res_name);
                    current_res_id = res_id;
                }

                builder.add_atom(serial, &name, &res_name, pos, Some(charge), Some(&ff_type));
            } else if line.starts_with("CONECT") || line.starts_with("ORDER") {
                connectivity_lines.push(line);
            } else if line.starts_with("END") {
                break;
            } else {
                metadata.extra_lines.push(line);
            }
        }

        for (line_num, conn_line) in connectivity_lines.iter().enumerate() {
            if conn_line.starts_with("CONECT") {
                let (base_serial, connected_serials) =
                    parse_conect_line(conn_line).map_err(|msg| BgfError::Parse {
                        line_num,
                        message: msg,
                    })?;
                for &connected_serial in &connected_serials {
                    builder.add_bond(base_serial, connected_serial, BondOrder::Single);
                }
            }
        }

        let system = builder.build();
        Ok((system, metadata))
    }

    fn write_to(
        system: &MolecularSystem,
        metadata: &Self::Metadata,
        writer: &mut impl Write,
    ) -> Result<(), Self::Error> {
        for line in &metadata.extra_lines {
            writeln!(writer, "{}", line)?;
        }

        for atom in system.atoms() {
            // TODO: Replace these placeholder values with actual calculations later.
            let atoms_connected = 0; // "Dumb" data for now
            let lone_pair = 0; // "Dumb" data for now

            let atom_line = format_atom_line(atom, atoms_connected, lone_pair);
            writeln!(writer, "{}", atom_line)?;
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

        let mut sorted_keys: Vec<_> = bond_map.keys().collect();
        sorted_keys.sort();

        for &key_idx in sorted_keys {
            let base_serial = system.get_atom(key_idx).unwrap().serial;
            let mut connect_line = format!("CONECT {:>5}", base_serial);
            for &neighbor_idx in &bond_map[&key_idx] {
                connect_line.push_str(&format!(
                    " {:>5}",
                    system.get_atom(neighbor_idx).unwrap().serial
                ));
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
) -> Result<(usize, String, String, char, isize, Point3<f64>, f64, String), String> {
    let serial = line
        .get(7..12)
        .and_then(|s| s.trim().parse().ok())
        .ok_or("Invalid serial")?;
    let name = line.get(13..17).unwrap_or("").trim().to_string();
    let res_name = line.get(19..22).unwrap_or("").trim().to_string();
    let chain_id = line
        .get(23..24)
        .and_then(|s| s.chars().next())
        .unwrap_or(' ');
    let res_id = line
        .get(25..29)
        .and_then(|s| s.trim().parse().ok())
        .ok_or("Invalid residue ID")?;
    let x = line
        .get(30..40)
        .and_then(|s| s.trim().parse().ok())
        .ok_or("Invalid X coordinate")?;
    let y = line
        .get(40..50)
        .and_then(|s| s.trim().parse().ok())
        .ok_or("Invalid Y coordinate")?;
    let z = line
        .get(50..60)
        .and_then(|s| s.trim().parse().ok())
        .ok_or("Invalid Z coordinate")?;
    let ff_type = line.get(61..66).unwrap_or("").trim().to_string();
    let charge = line
        .get(71..79)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0.0);

    Ok((
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
    if parts.is_empty() || parts[0] != "CONECT" {
        return Err("Not a CONECT line".to_string());
    }
    let base_serial = parts
        .get(1)
        .and_then(|s| s.parse().ok())
        .ok_or("Invalid base serial in CONECT")?;
    let connected_serials = parts[2..].iter().filter_map(|s| s.parse().ok()).collect();
    Ok((base_serial, connected_serials))
}

fn format_atom_line(
    atom: &crate::core::models::atom::Atom,
    atoms_connected: u8,
    lone_pair: u8,
) -> String {
    format!(
        "ATOM  {:>5} {:<4}  {:<3} {:1}{:>4}    {:>8.5}{:>10.5}{:>10.5} {:<5}   {:>1}  {:>1} {:>8.5}",
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
