# `GenericRotamer_defs.hpp` File Analysis

## File Purpose and Primary Role

This header file serves as a global definition and initialization module for generic amino acid rotamers in the SCREAM project. Its primary responsibility is to define a namespace containing pre-instantiated `AARotamer` objects for specific amino acids (currently Alanine 'A' and Cysteine 'C'). The file establishes a centralized location for rotamer library access, providing globally available rotamer objects that can be used throughout the SCREAM application for protein side-chain placement operations.

## Key Classes, Structs, and Functions (if any)

**Namespace:**

- `GenericRotamer`: Contains global rotamer definitions and path configuration

**Global Variables:**

- `GenericRotamerPath`: A string defining the base path to rotamer library files ("../lib/GenericRotamers/")
- `A`: A const pointer to an `AARotamer` object for Alanine, initialized with corresponding .lib and .cnn files
- `C`: A const pointer to an `AARotamer` object for Cysteine, initialized with corresponding .lib and .cnn files

**Note:** No classes, structs, or functions are defined in this file - it only contains variable declarations and object instantiations.

## Inputs

**Data Structures/Objects:**

- Depends on the `AARotamer` class constructor which likely takes file paths as input parameters

**File-Based Inputs:**

- `A.lib`: Rotamer library file for Alanine (likely contains rotamer conformations, chi angles, probabilities)
- `A.cnn`: Connection/connectivity file for Alanine (likely contains atomic connectivity information)
- `C.lib`: Rotamer library file for Cysteine
- `C.cnn`: Connection/connectivity file for Cysteine

**Environment Variables:**

- None directly used in this file

**Parameters/Configuration:**

- `GenericRotamerPath`: Hardcoded relative path ("../lib/GenericRotamers/") that determines where rotamer files are located
- File naming convention assumes single-letter amino acid codes with .lib and .cnn extensions

## Outputs

**Data Structures/Objects:**

- `const AARotamer*` objects: Global pointers to initialized rotamer objects for amino acids A and C
- These objects are available for use by other parts of the SCREAM system

**File-Based Outputs:**

- None - this file only reads data

**Console Output (stdout/stderr):**

- No direct console output from this file
- Potential error messages may come from `AARotamer` constructor if file loading fails

**Side Effects:**

- Creates global objects in memory at program startup
- Establishes global state that persists throughout program execution

## External Code Dependencies (Libraries/Headers)

**Standard C++ Library:**

- `<string>` (implied by use of `string` type, though not explicitly included - likely included via AARotamer.hpp)

**Internal SCREAM Project Headers:**

- `AARotamer.hpp`: Contains the `AARotamer` class definition and implementation

**External Compiled Libraries:**

- None apparent from this file

## Core Logic/Algorithm Flowchart (Mermaid JS Format)

```mermaid
graph TD
    A[Program Start] --> B[Include AARotamer.hpp]
    B --> C[Define GenericRotamer namespace]
    C --> D[Set GenericRotamerPath = "../lib/GenericRotamers/"]
    D --> E[Create AARotamer object for Alanine]
    E --> F[Load A.lib and A.cnn files]
    F --> G[Create AARotamer object for Cysteine]
    G --> H[Load C.lib and C.cnn files]
    H --> I[Objects available for global use]
    I --> J[Program continues with rotamer objects accessible]
```

## Potential Areas for Modernization/Refactoring in SCREAM++

1. **Replace Raw Pointers with Smart Pointers**: The current use of `new` for creating `AARotamer` objects with raw pointers creates potential memory leaks and unclear ownership semantics. Modern C++ should use `std::unique_ptr` or `std::shared_ptr` for automatic memory management.

2. **Eliminate Global Variables and Use Dependency Injection**: The current design relies on global variables which makes testing difficult and creates tight coupling. A better approach would be to use a rotamer factory or registry pattern with dependency injection, allowing for better testability and flexibility.

3. **Improve Path Configuration Management**: The hardcoded relative path should be replaced with a more flexible configuration system, possibly using environment variables, command-line arguments, or a configuration file. Consider using `std::filesystem::path` for cross-platform path handling and validation of file existence at initialization.
