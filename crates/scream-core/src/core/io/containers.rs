use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq)]
pub struct FileData {
    pub extra_lines: Vec<String>,
    pub extra_atom_columns: FormatSpecificData,
}

#[derive(Debug, Clone, PartialEq)]
pub enum FormatSpecificData {
    BGF(BgfExtraData), // BGF (Biograph) format specific data
}

#[derive(Debug, Clone, PartialEq)]
pub struct BgfExtraData {
    pub atoms_connected: HashMap<usize, usize>, // Map of atom serials to number of connections
    pub lone_pair: HashMap<usize, usize>,       // Map of atom serials to number of lone pairs
}
