# Error Handling & Exceptions

genepie maps GENESIS Fortran return codes onto a typed Python exception hierarchy,
and captures any message the Fortran layer writes to stderr. This means failures
are **catchable** — a bad input file or an unsupported option raises a Python
exception instead of aborting the interpreter.

## Exception hierarchy

```{mermaid}
flowchart TD
  GE["GenesisError"] --> GVE["GenesisValidationError<br/>(input checked in Python)"]
  GE --> GME["GenesisMemoryError"]
  GE --> GOE["GenesisOverflowError"]
  GE --> GLE["GenesisLibraryLoadError"]
  GE --> GFE["GenesisFortranError<br/>(has .code and .stderr_output)"]
  GFE --> M["GenesisFortranMemoryError (100-199)"]
  GFE --> F["GenesisFortranFileError (200-299)"]
  GFE --> V["GenesisFortranValidationError (300-399)"]
  GFE --> D["GenesisFortranDataError (400-499)"]
  GFE --> N["GenesisFortranNotSupportedError (500-599)"]
  GFE --> I["GenesisFortranInternalError (600-699)"]
```

- `GenesisValidationError` is raised **before** calling Fortran (e.g. a missing
  `cvfile`, an invalid selection, a lazy-mode restriction).
- `GenesisFortranError` and its subclasses come **from** Fortran; they carry a
  numeric `.code` and the captured `.stderr_output`.

## Error codes

| Code range | Category | Exception | Example codes |
|-----------:|----------|-----------|---------------|
| 100–199 | Memory | `GenesisFortranMemoryError` | `ERROR_ALLOC` (101), `ERROR_DEALLOC` (102) |
| 200–299 | File I/O | `GenesisFortranFileError` | `ERROR_FILE_NOT_FOUND` (201), `ERROR_FILE_FORMAT` (202) |
| 300–399 | Validation | `GenesisFortranValidationError` | `ERROR_INVALID_PARAM` (301), `ERROR_ATOM_COUNT` (304), `ERROR_SELECTION` (307) |
| 400–499 | Data | `GenesisFortranDataError` | `ERROR_DATA_MISMATCH` (401), `ERROR_NO_DATA` (402) |
| 500–599 | Not supported | `GenesisFortranNotSupportedError` | `ERROR_NOT_SUPPORTED` (501), `ERROR_DIM_NOT_SUPP` (502) |
| 600–699 | Internal | `GenesisFortranInternalError` | `ERROR_INTERNAL` (601), `ERROR_SYNTAX` (602) |

The numeric constants are available as `genepie.ErrorCode` attributes.

## Usage

```python
from genepie import genesis_exe, SMolecule
from genepie.exceptions import (
    GenesisError,
    GenesisFortranError,
    GenesisFortranNotSupportedError,
    GenesisValidationError,
    ErrorCode,
)

try:
    result = genesis_exe.rmsd_analysis(mol, traj, analysis_selection="an:CA")
except GenesisValidationError as e:
    # Bad input caught before Fortran ran.
    print("invalid input:", e)
except GenesisFortranNotSupportedError as e:
    print(f"unsupported feature (code {e.code}): {e}")
except GenesisFortranError as e:
    # Anything raised by the Fortran layer.
    print(f"Fortran error (code {e.code}): {e}")
    print("Fortran stderr:", e.stderr_output)
    if e.code == ErrorCode.ERROR_ATOM_COUNT:
        print("atom-count mismatch between molecule and trajectory")
except GenesisError as e:
    print("other GENESIS error:", e)
```

### A real example: rejecting an unsupported option

`wham_analysis`/`mbar_analysis` only accept precomputed reaction coordinates, so
passing a DCD file raises a *not supported* error rather than failing obscurely:

```python
from genepie.exceptions import GenesisFortranNotSupportedError

try:
    genesis_exe.mbar_analysis(cvfile="run{}.dat", dcdfile="whatever.dcd", ...)
except GenesisFortranNotSupportedError as e:
    print(f"MBAR does not support DCD input (code {e.code})")
```

### A real example: missing input files are caught, not fatal

A nonexistent `cvfile` used to reach Fortran's `open_file()`, which called
`exit(1)` and killed the process. genepie now validates the path up front
(`GenesisValidationError`), and even if that guard is bypassed, the library-mode
error guard converts the would-be `exit(1)` into a catchable
`GenesisFortranFileError` with `code == 201` (`ERROR_FILE_NOT_FOUND`).
