use crate::core::models::system::MolecularSystem;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Defines the interface for reading and writing molecular file formats.
///
/// This trait provides a common API for molecular file I/O operations,
/// supporting both reading from and writing to various file formats.
/// Implementors handle format-specific parsing and serialization.
pub trait MolecularFile {
    /// The type of metadata associated with the file format.
    type Metadata;

    /// The error type for I/O operations.
    type Error: Error + From<io::Error>;

    /// Reads a molecular system from a buffered reader.
    ///
    /// # Arguments
    ///
    /// * `reader` - The buffered reader to read from.
    ///
    /// # Return
    ///
    /// Returns the parsed molecular system and associated metadata.
    ///
    /// # Errors
    ///
    /// Returns an error if parsing fails or I/O operations encounter issues.
    fn read_from(
        reader: &mut impl BufRead,
    ) -> Result<(MolecularSystem, Self::Metadata), Self::Error>;

    /// Writes a molecular system and metadata to a writer.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system to write.
    /// * `metadata` - The metadata to include in the output.
    /// * `writer` - The writer to output to.
    ///
    /// # Return
    ///
    /// Returns `Ok(())` on success.
    ///
    /// # Errors
    ///
    /// Returns an error if writing fails or I/O operations encounter issues.
    fn write_to(
        system: &MolecularSystem,
        metadata: &Self::Metadata,
        writer: &mut impl Write,
    ) -> Result<(), Self::Error>;

    /// Writes a molecular system to a writer without metadata.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system to write.
    /// * `writer` - The writer to output to.
    ///
    /// # Return
    ///
    /// Returns `Ok(())` on success.
    ///
    /// # Errors
    ///
    /// Returns an error if writing fails or I/O operations encounter issues.
    fn write_system_to(
        system: &MolecularSystem,
        writer: &mut impl Write,
    ) -> Result<(), Self::Error>;

    /// Reads a molecular system from a file path.
    ///
    /// # Arguments
    ///
    /// * `path` - The path to the file to read.
    ///
    /// # Return
    ///
    /// Returns the parsed molecular system and associated metadata.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened or parsing fails.
    fn read_from_path<P: AsRef<Path>>(
        path: P,
    ) -> Result<(MolecularSystem, Self::Metadata), Self::Error> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        Self::read_from(&mut reader)
    }

    /// Writes a molecular system and metadata to a file path.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system to write.
    /// * `metadata` - The metadata to include in the output.
    /// * `path` - The path to the file to write.
    ///
    /// # Return
    ///
    /// Returns `Ok(())` on success.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created or writing fails.
    fn write_to_path<P: AsRef<Path>>(
        system: &MolecularSystem,
        metadata: &Self::Metadata,
        path: P,
    ) -> Result<(), Self::Error> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        Self::write_to(system, metadata, &mut writer)
    }

    /// Writes a molecular system to a file path without metadata.
    ///
    /// # Arguments
    ///
    /// * `system` - The molecular system to write.
    /// * `path` - The path to the file to write.
    ///
    /// # Return
    ///
    /// Returns `Ok(())` on success.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created or writing fails.
    fn write_system_to_path<P: AsRef<Path>>(
        system: &MolecularSystem,
        path: P,
    ) -> Result<(), Self::Error> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        Self::write_system_to(system, &mut writer)
    }
}
