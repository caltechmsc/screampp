use super::sorting::sorter::sort_system_atoms;
use super::traits::MolecularFile;
use crate::core::models::atom::Atom;
use crate::core::models::chain::{Chain, ChainType};
use crate::core::models::ids::{AtomId, ChainId, ResidueId};
use crate::core::models::residue::{Residue, ResidueType};
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
    #[error("Logic error during file processing: {0}")]
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
                "ORDER" => metadata.order_lines.push(line.clone()),
                "FORMAT" => {}
                _ => {
                    if !atom_section_started {
                        metadata.header_lines.push(line.clone());
                    } else {
                        metadata.footer_lines.push(line.clone());
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

            let atom = Atom::new(&atom_info.name, active_residue_id, atom_info.pos);

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
            let lone_pairs = 0;

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

        for line in &metadata.footer_lines {
            writeln!(writer, "{}", line)?;
        }

        // --- Phase 4: Write CONECT records using the new serial numbers ---
        let mut adjacency_list = HashMap::<usize, Vec<usize>>::new();
        for bond in system.bonds() {
            if let (Some(&serial1), Some(&serial2)) = (
                atom_id_to_new_serial.get(&bond.atom1_id),
                atom_id_to_new_serial.get(&bond.atom2_id),
            ) {
                adjacency_list.entry(serial1).or_default().push(serial2);
                adjacency_list.entry(serial2).or_default().push(serial1);
            }
        }

        let mut sorted_base_serials: Vec<_> = adjacency_list.keys().copied().collect();
        sorted_base_serials.sort_unstable();

        for base_serial in sorted_base_serials {
            let neighbors = adjacency_list.get_mut(&base_serial).unwrap();
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
