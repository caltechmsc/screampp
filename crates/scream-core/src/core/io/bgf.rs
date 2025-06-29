use super::traits::MolecularFile;
use crate::core::models::atom::{Atom, Element};
use crate::core::models::chain::ChainType;
use crate::core::models::system::{MolecularSystem, MolecularSystemBuilder};
use crate::core::models::topology::BondOrder;
use nalgebra::Point3;
use std::collections::{BTreeMap, HashMap};
use std::io::{self, BufRead, ErrorKind, Write};
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
    pub header_lines: BTreeMap<usize, RawLine>, // File-level data
    pub atom_io_data: HashMap<usize, BgfAtomIoData>, // Per-atom extra column data
    pub format_lines: Vec<String>,              // Explicit FORMAT lines
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
    #[error("Invalid number format in column {column}: {source}")]
    InvalidNumber {
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

    #[error("ATOM/HETATM line has an invalid number of columns (expected at least 12)")]
    InvalidAtomColumnCount,

    #[error("CONECT line requires at least two atoms")]
    InvalidConectFormat,
}
