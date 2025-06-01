# `sc_AASideChain.cpp` File Analysis

## File Purpose and Primary Role

This file implements the `AASideChain` class, which represents amino acid side chains in the SCREAM molecular modeling software. It inherits from the `SideChain` base class and provides specific functionality for handling amino acid side chain atoms, including automatic amino acid type detection, atom retrieval by Greek letter nomenclature, and basic side chain management operations. The class serves as a specialized container for amino acid side chain atoms with knowledge of protein-specific naming conventions.

## Key Classes, Structs, and Functions (if any)

### Classes

- **`AASideChain`**: The primary class that extends `SideChain` to handle amino acid-specific side chain operations
  - Provides constructors for empty initialization and atom list-based initialization
  - Implements copy constructor and assignment operator
  - Contains methods for atom retrieval using Greek letter naming conventions

### Key Methods

- **`AASideChain()`**: Default constructor
- **`AASideChain(const ScreamAtomV& atom_list)`**: Constructor that initializes from atom list and auto-detects amino acid type
- **`get_atom_by_greek_name(string alphabet)`**: Retrieves atoms using Greek letter nomenclature (β, γ, δ, ε, ζ, η)
- **`dummy_assignment_operator(const AASideChain& sc)`**: Private helper for assignment operations

## Inputs

### Data Structures/Objects

- **`ScreamAtomV`**: Vector of `SCREAM_ATOM` pointers used in constructor to initialize the side chain
- **`AASideChain`**: Other instances for copy constructor and assignment operator
- **`string`**: Greek letter identifiers ("B", "G", "D", "E", "Z", "H") for atom retrieval

### File-Based Inputs

- No direct file I/O operations are performed in this file

### Environment Variables

- No direct environment variable usage detected

### Parameters/Configuration

- The class relies on atom naming conventions consistent with PDB format
- Uses residue names stored in `SCREAM_ATOM` objects for amino acid type detection

## Outputs

### Data Structures/Objects

- **`SCREAM_ATOM*`**: Pointers to specific atoms retrieved by Greek letter names
- **`AASideChain`**: New instances created through constructors and assignment
- **`string`**: Amino acid type (SC_name) automatically determined from input atoms

### File-Based Outputs

- No direct file output operations

### Console Output (stdout/stderr)

- No console output in this file

### Side Effects

- Modifies internal `sc_atom_mm` multimap and `SC_name` during construction and assignment
- Memory management through inherited `SideChain` functionality

## External Code Dependencies

### Standard C++ Library

- `<vector>`: For atom list containers
- `<map>`: For atom storage multimap
- `<string>`: For atom names and amino acid type identification
- `<sstream>`: Included but not directly used in this file
- `<iostream>`: Included but not directly used in this file
- `<typeinfo>`: For type information (included but not used)

### Internal SCREAM Project Headers

- `"scream_atom.hpp"`: Core atom class definitions
- `"Rotlib.hpp"`: Rotamer library functionality
- `"defs.hpp"`: Project-wide definitions
- `"sc_SideChain.hpp"`: Base class for side chain functionality
- `"sc_AASideChain.hpp"`: Header for this class
- `"scream_tools.hpp"`: Utility functions for string processing and atom checking

### External Compiled Libraries

- None detected in this file

## Core Logic/Algorithm Flowchart (Mermaid JS Format)

```mermaid
graph TD
    A[AASideChain Constructor Called] --> B{Empty Constructor?}
    B -- Yes --> C[Initialize with SideChain default]
    B -- No --> D[Receive ScreamAtomV atom_list]

    D --> E{atom_list.size() == 1?}
    E -- Yes --> F{Is HCA atom?}
    F -- Yes --> G[Set SC_name = "GLY"]
    F -- No --> H[Get resName from first atom]
    E -- No --> H

    G --> I[Initialize sc_atom_mm]
    H --> I
    I --> J[Construction Complete]

    K[get_atom_by_greek_name Called] --> L{Which Greek Letter?}
    L -- B --> M[Return CB atom]
    L -- G --> N[Search for OG/OG1/SG/CG/CG1]
    L -- D --> O[Search for OD/OD1/SD/ND1/CD/CD1]
    L -- E --> P[Search for OE/OE1/NE/CE/CE1]
    L -- Z --> Q[Search for NZ/CZ]
    L -- H --> R[Search for NH1]
    L -- Other --> S[Return NULL]

    N --> T{Found?}
    O --> T
    P --> T
    Q --> T
    R --> T
    T -- Yes --> U[Return SCREAM_ATOM*]
    T -- No --> S

    C --> J
    M --> U
    S --> V[Return NULL]
```

## Potential Areas for Modernization/Refactoring in SCREAM++

### 1. **Smart Pointer Usage and Memory Safety**

The current code uses raw pointers (`SCREAM_ATOM*`) extensively, which can lead to memory management issues. Modern C++ should utilize `std::shared_ptr<ScreamAtom>` or `std::unique_ptr<ScreamAtom>` to ensure automatic memory management and prevent memory leaks. The multimap storage should be updated to `std::multimap<std::string, std::shared_ptr<ScreamAtom>>`.

### 2. **Exception Handling and Error Management**

The current implementation returns `NULL` pointers when atoms are not found and uses comments like "change this to exception handling later." This should be replaced with a proper exception hierarchy (e.g., `AtomNotFoundException`) and/or `std::optional<ScreamAtom*>` return types to make error conditions explicit and force calling code to handle them appropriately.

### 3. **Modern C++ String and Container Usage**

The code should leverage modern C++ features like `std::string_view` for string parameters to avoid unnecessary copies, `constexpr` for compile-time constants, and range-based for loops. The Greek letter matching logic could be refactored using `std::unordered_set` for O(1) lookups or even better, an `enum class` for type-safe Greek letter identification instead of string comparisons.
