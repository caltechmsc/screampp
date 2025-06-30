use super::traits::MolecularFile;
use crate::core::models::atom::Element;
use crate::core::models::chain::ChainType;
use crate::core::models::system::{MolecularSystem, MolecularSystemBuilder};
use crate::core::models::topology::{Bond, BondOrder};
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
    pub extra_columns: BTreeMap<usize, String>,
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
    #[error("Invalid integer format in column {column}: {source}")]
    InvalidInt {
        column: usize,
        #[source]
        source: std::num::ParseIntError,
    },
    #[error("Invalid float format in column {column}: {source}")]
    InvalidFloat {
        column: usize,
        #[source]
        source: std::num::ParseFloatError,
    },
    #[error("Invalid element symbol '{symbol}'")]
    InvalidElement { symbol: String },
    #[error("ATOM/HETATM line has an insufficient number of columns (at least 12 required)")]
    InvalidAtomColumnCount,
    #[error("CONECT line requires at least two atoms")]
    InvalidConectFormat,
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
            let parts: Vec<&str> = line.split_whitespace().collect();

            if parts.is_empty() {
                metadata
                    .header_lines
                    .insert(line_num, RawLine { content: line });
                continue;
            }

            let record_type = parts[0];
            match record_type {
                "ATOM" | "HETATM" => {
                    if parts.len() < 12 {
                        return Err(BgfError::Parse {
                            line: line_num,
                            kind: BgfParseErrorKind::InvalidAtomColumnCount,
                        });
                    }

                    let serial: usize = parts[1].parse().map_err(|e| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidInt {
                            column: 2,
                            source: e,
                        },
                    })?;
                    if !seen_serials.insert(serial) {
                        return Err(BgfError::Inconsistency(format!(
                            "Duplicate atom serial: {}",
                            serial
                        )));
                    }

                    let name = parts[2];
                    let res_name = parts[3];
                    let chain_id: char = parts[4].chars().next().ok_or_else(|| {
                        BgfError::Inconsistency(format!("Missing chain ID on line {}", line_num))
                    })?;
                    let res_id: isize = parts[5].parse().map_err(|e| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidInt {
                            column: 6,
                            source: e,
                        },
                    })?;
                    let x: f64 = parts[6].parse().map_err(|e| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidFloat {
                            column: 7,
                            source: e,
                        },
                    })?;
                    let y: f64 = parts[7].parse().map_err(|e| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidFloat {
                            column: 8,
                            source: e,
                        },
                    })?;
                    let element_str = parts[8];
                    let element: Element = element_str.parse().map_err(|_| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidElement {
                            symbol: element_str.to_string(),
                        },
                    })?;
                    let ff_type = parts[9];
                    let charge: f64 = parts[10].parse().map_err(|e| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidFloat {
                            column: 11,
                            source: e,
                        },
                    })?;
                    let z: f64 = parts[11].parse().map_err(|e| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidFloat {
                            column: 12,
                            source: e,
                        },
                    })?;
                    let position = Point3::new(x, y, z);

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
                    }
                    if res_id != current_residue_id {
                        builder.start_residue(res_id, res_name);
                        current_residue_id = res_id;
                    }
                    builder.add_atom(serial, name, element, position, charge, ff_type);

                    if parts.len() > 12 {
                        let extra_cols = parts[12..]
                            .iter()
                            .enumerate()
                            .map(|(i, &s)| (i + 13, s.to_string()))
                            .collect();
                        metadata.atom_io_data.insert(
                            serial,
                            BgfAtomIoData {
                                extra_columns: extra_cols,
                            },
                        );
                    }
                }
                "CONECT" => {
                    if parts.len() < 3 {
                        return Err(BgfError::Parse {
                            line: line_num,
                            kind: BgfParseErrorKind::InvalidConectFormat,
                        });
                    }
                    let atom1_serial: usize = parts[1].parse().map_err(|e| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidInt {
                            column: 2,
                            source: e,
                        },
                    })?;
                    for part in &parts[2..] {
                        let atom2_serial: usize = part.parse().map_err(|e| BgfError::Parse {
                            line: line_num,
                            kind: BgfParseErrorKind::InvalidInt {
                                column: 0,
                                source: e,
                            },
                        })?;
                        temp_conect.push((
                            atom1_serial.min(atom2_serial),
                            atom1_serial.max(atom2_serial),
                        ));
                    }
                }
                "ORDER" => {
                    if parts.len() < 3 {
                        return Err(BgfError::Inconsistency(format!(
                            "ORDER record on line {} is incomplete",
                            line_num
                        )));
                    }
                    let a1: usize = parts[1].parse().map_err(|e| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidInt {
                            column: 2,
                            source: e,
                        },
                    })?;
                    let a2: usize = parts[2].parse().map_err(|e| BgfError::Parse {
                        line: line_num,
                        kind: BgfParseErrorKind::InvalidInt {
                            column: 3,
                            source: e,
                        },
                    })?;
                    let order = parts
                        .get(3)
                        .unwrap_or(&"1")
                        .parse()
                        .unwrap_or(BondOrder::Single);
                    temp_orders.insert((a1.min(a2), a1.max(a2)), order);
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

        let mut atom_meta_map: HashMap<usize, (char, isize, &str, ChainType)> = HashMap::new();
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
            let (chain_id, res_id, res_name, chain_type) =
                atom_meta_map.get(&atom.index).ok_or_else(|| {
                    BgfError::Inconsistency(format!(
                        "Atom index {} not found in any residue",
                        atom.index
                    ))
                })?;
            let record_type = if *chain_type == ChainType::Protein {
                "ATOM  "
            } else {
                "HETATM"
            };

            let mut line = format!(
                "{}{:>5} {:<5} {:<3} {} {:>5}    {:8.3}{:8.3}{:8.3} {:<2} {:<5} {:8.4}",
                record_type,
                atom.serial,
                atom.name,
                res_name,
                chain_id,
                res_id,
                atom.position.x,
                atom.position.y,
                atom.position.z,
                atom.element,
                atom.force_field_type,
                atom.partial_charge
            );

            if let Some(io_data) = metadata.atom_io_data.get(&atom.serial) {
                for (_, val) in &io_data.extra_columns {
                    line.push_str(&format!(" {}", val));
                }
            }
            writeln!(writer, "{}", line)?;
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
