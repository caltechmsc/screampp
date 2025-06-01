# `scream_rtf.cpp` File Analysis

## File Purpose and Primary Role

This file implements the RTF (Residue Topology File) parser for the SCREAM molecular modeling software. It reads and parses NAMD/CHARMM-compatible RTF files that define the topology and force field parameters for amino acid residues. The primary role is to provide a data structure that maps residue names to their atomic compositions, force field types, and bonding information, which is essential for molecular modeling and energy calculations.

## Key Classes, Structs, and Functions (if any)

### SCREAM_RTF Class

- **Purpose**: Main class that manages the entire RTF file parsing and provides access to residue information
- **Key Methods**:
  - `SCREAM_RTF(string rtf_file)`: Constructor that initializes by parsing an RTF file
  - `get_AminoAcid_RTF(string resName)`: Returns the RTF data for a specific residue
  - `get_ff_type(string resName, string atomLabel)`: Gets the force field type for a specific atom in a residue
  - `_init(string rtf_file)`: Private method that performs the actual file parsing

### AminoAcid_RTF Class

- **Purpose**: Represents the topology information for a single amino acid residue
- **Key Members**:
  - `ff_type`: Map storing atom labels to their force field types
  - `bonds`: Container storing bonding information between atoms
  - `resName`: Name of the residue
- **Key Methods**:
  - `AminoAcid_RTF(stringV& lines)`: Constructor that parses residue data from text lines
  - `get_ff_type(string atom_label)`: Returns the force field type for a given atom
  - `_init(stringV& lines)`: Private method that parses residue-specific data

## Inputs

### Data Structures/Objects

- `stringV` (vector of strings): Used to pass lines of text data for parsing individual residues
- `string`: Residue names and atom labels for lookup operations

### File-Based Inputs

- **RTF files**: NAMD/CHARMM-compatible residue topology files containing:
  - Residue definitions (RESI entries)
  - Atom definitions with force field types and charges
  - Bond connectivity information
  - Other topology parameters (DONOR, ACCEPTOR, IC - currently ignored)

### Environment Variables

- No direct environment variable dependencies identified in this file

### Parameters/Configuration

- RTF file path provided during construction
- Residue names and atom labels used for data retrieval

## Outputs

### Data Structures/Objects

- `SCREAM_RTF` object: Complete RTF database accessible via residue names
- `AminoAcid_RTF*` pointers: Individual residue topology objects
- `string`: Force field type strings for specific atoms

### File-Based Outputs

- No direct file output from this module

### Console Output (stdout/stderr)

- Debug messages showing parsing progress and residue creation
- Error messages for file access failures (sent to stderr)
- Debug output includes parsed lines and newly created residue objects

### Side Effects

- Memory allocation for `AminoAcid_RTF` objects (managed via map in `_rtfTable`)
- Program termination (`exit(2)`) on file access failure

## External Code Dependencies (Libraries/Headers)

### Standard C++ Library

- `<cstdlib>`: For `exit()` function
- `<iostream>`: For console I/O operations
- `<fstream>`: For file reading operations
- `<map>`: For storing residue lookup table (implied from usage)
- `<string>`: For string operations (implied from usage)

### Internal SCREAM Project Headers

- `"scream_rtf.hpp"`: Header file containing class declarations
- Implied dependencies on:
  - Debug class (for debug output)
  - `stringV` type definition
  - `split()` function for string parsing

### External Compiled Libraries

- None identified

## Core Logic/Algorithm Flowchart (Mermaid JS Format)

```mermaid
graph TD
    A[Start: SCREAM_RTF Constructor] --> B[Open RTF file]
    B --> C{File accessible?}
    C -- No --> D[Print error & exit(2)]
    C -- Yes --> E[Initialize parsing variables]
    E --> F[Read line from file]
    F --> G{End of file?}
    G -- Yes --> H[Process final residue]
    G -- No --> I{Comment line?}
    I -- Yes --> F
    I -- No --> J{RESI line?}
    J -- Yes --> K{First residue?}
    K -- Yes --> L[Set current residue name]
    K -- No --> M[Create AminoAcid_RTF for previous residue]
    M --> N[Store in _rtfTable]
    N --> L
    L --> O[Clear residue lines buffer]
    O --> P[Add line to buffer]
    P --> F
    J -- No --> P
    H --> Q[Create final AminoAcid_RTF]
    Q --> R[Store in _rtfTable]
    R --> S[End]

    T[AminoAcid_RTF::_init] --> U[Parse each line]
    U --> V{RESI line?}
    V -- Yes --> W[Extract residue name]
    V -- No --> X{ATOM line?}
    X -- Yes --> Y[Store atom->ff_type mapping]
    X -- No --> Z{BOND/DOUBLE line?}
    Z -- Yes --> AA[Parse bond pairs]
    Z -- No --> BB{PRES line?}
    BB -- Yes --> CC[Break parsing]
    BB -- No --> U
    Y --> U
    AA --> U
    W --> U
    CC --> DD[End]
```

## Potential Areas for Modernization/Refactoring in SCREAM++

### 1. Memory Management and RAII

- **Current Issue**: Manual `new`/`delete` operations with raw pointers in destructor
- **Modernization**: Use `std::unique_ptr` or `std::shared_ptr` for automatic memory management, eliminating the need for manual cleanup in the destructor

### 2. Error Handling and Resource Management

- **Current Issue**: File access failure results in immediate program termination with `exit(2)`
- **Modernization**: Use exception-based error handling with RAII file management (`std::ifstream` with proper exception handling) to allow graceful error recovery and better integration with larger applications

### 3. String Processing and Type Safety

- **Current Issue**: Heavy reliance on C-style string manipulation and custom `stringV` type
- **Modernization**: Utilize modern C++ string processing with `std::string_view`, `std::regex` for parsing, and standard containers (`std::vector<std::string>`) to improve performance and type safety while reducing custom type dependencies
