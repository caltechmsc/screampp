use crate::core::models::system::MolecularSystem;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub trait MolecularFile {
    type Metadata;
    type Error: Error + From<io::Error>;

    fn read_from(
        reader: &mut impl BufRead,
    ) -> Result<(MolecularSystem, Self::Metadata), Self::Error>;

    fn write_to(
        system: &MolecularSystem,
        metadata: &Self::Metadata,
        writer: &mut impl Write,
    ) -> Result<(), Self::Error>;

    fn write_system_to(
        system: &MolecularSystem,
        writer: &mut impl Write,
    ) -> Result<(), Self::Error>;

    fn read_from_path<P: AsRef<Path>>(
        path: P,
    ) -> Result<(MolecularSystem, Self::Metadata), Self::Error> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        Self::read_from(&mut reader)
    }

    fn write_to_path<P: AsRef<Path>>(
        system: &MolecularSystem,
        metadata: &Self::Metadata,
        path: P,
    ) -> Result<(), Self::Error> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        Self::write_to(system, metadata, &mut writer)
    }

    fn write_system_to_path<P: AsRef<Path>>(
        system: &MolecularSystem,
        path: P,
    ) -> Result<(), Self::Error> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        Self::write_system_to(system, &mut writer)
    }
}
