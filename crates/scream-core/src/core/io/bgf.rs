use crate::core::io::traits::MolecularFile;
use crate::core::models::chain::ChainType;
use crate::core::models::system::{MolecularSystem, MolecularSystemBuilder};
use crate::core::models::topology::BondOrder;
use nalgebra::Point3;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{self, BufRead, Write};
use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct RawLine {
    pub content: String,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct BgfAtomIoData {
    pub raw_suffix: String,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct BgfMetadata {
    pub header_lines: BTreeMap<usize, RawLine>,
    pub atom_io_data: HashMap<usize, BgfAtomIoData>,
    pub format_lines: Vec<String>,
}

#[derive(Debug, Error)]
pub enum BgfError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    #[error("Parse error on line {line}: {kind}")]
    Parse {
        line: usize,
        kind: BgfParseErrorKind,
    },
    #[error("Inconsistent data: {0}")]
    Inconsistency(String),
    #[error("Missing required record: {0}")]
    MissingRecord(String),
}

#[derive(Debug, Error)]
pub enum BgfParseErrorKind {
    #[error("Invalid integer format in columns {columns} (value: '{value}')")]
    InvalidInt { columns: String, value: String },
    #[error("Invalid float format in columns {columns} (value: '{value}')")]
    InvalidFloat { columns: String, value: String },
    #[error("Required field in columns {columns} is empty")]
    MissingRequiredField { columns: String },
    #[error("Line is too short for ATOM/HETATM record (must be at least 80 chars)")]
    LineTooShort,
    #[error("CONECT line requires at least two atoms")]
    InvalidConectFormat,
}

fn slice_and_trim(line: &str, start: usize, end: usize) -> &str {
    line.get(start..end).unwrap_or("").trim()
}

pub struct BgfFile;

impl MolecularFile for BgfFile {
    type Metadata = BgfMetadata;
    type Error = BgfError;

    fn read_from(
        reader: &mut impl BufRead,
    ) -> Result<(MolecularSystem, Self::Metadata), Self::Error> {
        let mut builder = MolecularSystemBuilder::new();
        let mut metadata = BgfMetadata::default();
        let mut seen_serials = HashSet::new();

        let mut temp_conect: Vec<(usize, usize)> = Vec::new();
        let mut temp_orders: HashMap<(usize, usize), BondOrder> = HashMap::new();

        let mut current_chain_id = '\0';
        let mut current_residue_id = isize::MIN;
        let mut chain_is_hetero: HashMap<char, bool> = HashMap::new();

        for (line_num, line_res) in reader.lines().enumerate() {
            let line = line_res?;
            let line_num = line_num + 1;

            let record_type = slice_and_trim(&line, 0, 6);
            if record_type.is_empty() {
                if !line.trim().is_empty() {
                    metadata
                        .header_lines
                        .insert(line_num, RawLine { content: line });
                }
                continue;
            }

            match record_type {
                "ATOM" | "HETATM" => {
                    if line.len() < 80 {
                        return Err(BgfError::Parse {
                            line: line_num,
                            kind: BgfParseErrorKind::LineTooShort,
                        });
                    }

                    let serial_str = slice_and_trim(&line, 7, 12);
                    let name_str = slice_and_trim(&line, 13, 18);
                    let res_name_str = slice_and_trim(&line, 19, 22);
                    let chain_id_str = slice_and_trim(&line, 23, 24);
                    let res_id_str = slice_and_trim(&line, 25, 30);
                    let x_str = slice_and_trim(&line, 30, 40);
                    let y_str = slice_and_trim(&line, 40, 50);
                    let z_str = slice_and_trim(&line, 50, 60);
                    let ff_type_str = slice_and_trim(&line, 61, 66);
                    let charge_str = slice_and_trim(&line, 72, 80);

                    if name_str.is_empty() {
                        return Err(BgfError::Parse {
                            line: line_num,
                            kind: BgfParseErrorKind::MissingRequiredField {
                                columns: "14-18".into(),
                            },
                        });
                    }
                    let serial: usize = serial_str.parse().map_err(|_| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidInt {
                            columns: "8-12".into(),
                            value: serial_str.into(),
                        },
                    })?;
                    if !seen_serials.insert(serial) {
                        return Err(BgfError::Inconsistency(format!(
                            "Duplicate atom serial: {}",
                            serial
                        )));
                    }

                    let chain_id: char = chain_id_str.chars().next().unwrap_or('A');
                    let res_id: isize = res_id_str.parse().map_err(|_| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidInt {
                            columns: "26-30".into(),
                            value: res_id_str.into(),
                        },
                    })?;
                    let x: f64 = x_str.parse().map_err(|_| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidFloat {
                            columns: "31-40".into(),
                            value: x_str.into(),
                        },
                    })?;
                    let y: f64 = y_str.parse().map_err(|_| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidFloat {
                            columns: "41-50".into(),
                            value: y_str.into(),
                        },
                    })?;
                    let z: f64 = z_str.parse().map_err(|_| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidFloat {
                            columns: "51-60".into(),
                            value: z_str.into(),
                        },
                    })?;
                    if ff_type_str.is_empty() {
                        return Err(BgfError::Parse {
                            line: line_num,
                            kind: BgfParseErrorKind::MissingRequiredField {
                                columns: "62-66".into(),
                            },
                        });
                    }
                    let charge: f64 = charge_str.parse().map_err(|_| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidFloat {
                            columns: "73-80".into(),
                            value: charge_str.into(),
                        },
                    })?;

                    if chain_id != current_chain_id {
                        let is_hetero = chain_is_hetero.entry(chain_id).or_insert(false);
                        if record_type == "HETATM" {
                            *is_hetero = true;
                        }
                        let chain_type = if *is_hetero {
                            ChainType::Other
                        } else {
                            ChainType::Protein
                        };
                        builder.start_chain(chain_id, chain_type);
                        current_chain_id = chain_id;
                        current_residue_id = isize::MIN;
                    }
                    if res_id != current_residue_id {
                        builder.start_residue(res_id, res_name_str);
                        current_residue_id = res_id;
                    }
                    builder.add_atom(serial, name_str, Point3::new(x, y, z), charge, ff_type_str);

                    metadata.atom_io_data.insert(
                        serial,
                        BgfAtomIoData {
                            raw_suffix: line[80..].to_string(),
                        },
                    );
                }
                "CONECT" | "ORDER" => {
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() < 3 {
                        continue;
                    }
                    let a1_res = parts[1].parse::<usize>();
                    let a2_res = parts[2].parse::<usize>();
                    if let (Ok(a1), Ok(a2)) = (a1_res, a2_res) {
                        let key = (a1.min(a2), a1.max(a2));
                        if record_type == "CONECT" {
                            temp_conect.push(key);
                        } else {
                            let order = if parts.len() > 3 {
                                parts[3].parse().unwrap_or_default()
                            } else {
                                BondOrder::default()
                            };
                            temp_orders.insert(key, order);
                        }
                    }
                }
                "FORMAT" => metadata.format_lines.push(line.clone()),
                "END" => break,
                _ => {
                    metadata.header_lines.insert(
                        line_num,
                        RawLine {
                            content: line.clone(),
                        },
                    );
                }
            }
        }

        temp_conect.sort_unstable();
        temp_conect.dedup();

        for (a1_serial, a2_serial) in temp_conect {
            let order = temp_orders
                .get(&(a1_serial, a2_serial))
                .copied()
                .unwrap_or_default();
            builder.add_bond(a1_serial, a2_serial, order);
        }
        if seen_serials.is_empty() {
            return Err(BgfError::MissingRecord("ATOM/HETATM records".into()));
        }
        Ok((builder.build(), metadata))
    }

    fn write_to(
        system: &MolecularSystem,
        metadata: &Self::Metadata,
        writer: &mut impl Write,
    ) -> Result<(), Self::Error> {
        for (_, line) in &metadata.header_lines {
            writeln!(writer, "{}", line.content)?;
        }
        for line in &metadata.format_lines {
            writeln!(writer, "{}", line)?;
        }

        let mut atom_meta_map = BTreeMap::new();
        for chain in system.chains() {
            for residue in chain.residues() {
                for &atom_idx in &residue.atom_indices {
                    atom_meta_map.insert(
                        atom_idx,
                        (chain.id, residue.id, &residue.name, chain.chain_type),
                    );
                }
            }
        }

        for atom in system.atoms() {
            let (chain_id, res_id, res_name, chain_type) = atom_meta_map.get(&atom.index).unwrap();
            let record_type = if *chain_type == ChainType::Protein {
                "ATOM"
            } else {
                "HETATM"
            };

            let line = format!(
                "{:<6} {:>5} {:<5} {:>3}  {:1}{:>4}    {:>8.3}{:>8.3}{:>8.3}{}{}{:>6}{:>12.5}",
                record_type,
                atom.serial,
                atom.name,
                res_name,
                chain_id,
                res_id,
                atom.position.x,
                atom.position.y,
                atom.position.z,
                "  1.00",
                "  0.00",
                atom.force_field_type,
                atom.partial_charge
            );
            let suffix = metadata
                .atom_io_data
                .get(&atom.serial)
                .map_or("", |d| &d.raw_suffix);
            writeln!(writer, "{}{}", line, suffix)?;
        }

        if !system.bonds().is_empty() {
            let mut bond_map: BTreeMap<usize, Vec<(usize, BondOrder)>> = BTreeMap::new();
            for bond in system.bonds() {
                let atom1 = system.get_atom(bond.atom1_idx).ok_or_else(|| {
                    BgfError::Inconsistency(format!("Bond atom index {} not found", bond.atom1_idx))
                })?;
                let atom2 = system.get_atom(bond.atom2_idx).ok_or_else(|| {
                    BgfError::Inconsistency(format!("Bond atom index {} not found", bond.atom2_idx))
                })?;
                bond_map
                    .entry(atom1.serial)
                    .or_default()
                    .push((atom2.serial, bond.order));
                bond_map
                    .entry(atom2.serial)
                    .or_default()
                    .push((atom1.serial, bond.order));
            }
            for (a1_serial, conns) in &bond_map {
                write!(writer, "CONECT {:>5}", a1_serial)?;
                for (a2_serial, _) in conns {
                    write!(writer, " {:>5}", a2_serial)?;
                }
                writeln!(writer)?;
            }
            for (a1_serial, conns) in &bond_map {
                let non_single: Vec<_> = conns
                    .iter()
                    .filter(|(_, o)| *o != BondOrder::Single)
                    .collect();
                if !non_single.is_empty() {
                    write!(writer, "ORDER  {:>5}", a1_serial)?;
                    for (_, order) in non_single {
                        write!(writer, " {:>5}", format!("{}", *order as u8))?;
                    }
                    writeln!(writer)?;
                }
            }
        }

        writeln!(writer, "END")?;
        Ok(())
    }

    fn write_system_to(
        system: &MolecularSystem,
        writer: &mut impl Write,
    ) -> Result<(), Self::Error> {
        let mut header_lines = BTreeMap::new();
        header_lines.insert(
            1,
            RawLine {
                content: "REMARK Generated by scream-core".to_string(),
            },
        );
        header_lines.insert(
            2,
            RawLine {
                content: "FORCEFIELD DREIDING".to_string(),
            },
        );
        let default_metadata = BgfMetadata {
            header_lines,
            ..Default::default()
        };
        Self::write_to(system, &default_metadata, writer)
    }
}
