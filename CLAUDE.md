# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GENESIS (GENeralized-Ensemble SImulation System) is a molecular dynamics simulation software for biomolecular systems. This repository is focused on developing the **Python interface** for GENESIS analysis tools, located in `src/genepie/`.

## Build Commands

### Quick Start (Python Interface Development)

```bash
# Set up Python environment
uv venv --python=python3.11
source .venv/bin/activate
uv pip install torch torchvision torchaudio nglview numpy mdtraj MDAnalysis plotly jupyterlab py3Dmol scikit-learn gdown

# Build GENESIS (first time or after configure.ac changes)
autoreconf -fi
./configure --disable-mpi CC=gcc-14 FC=gfortran LAPACK_LIBS="-L/usr/local/lib -llapack -lblas"
make -j$(sysctl -n hw.ncpu)

# Install Python package (editable mode)
uv pip install -e .
```

### Rebuilding After Changes

| Change Type | Required Action |
|-------------|-----------------|
| Python files (.py) | No rebuild needed (instant reflection) |
| Fortran files (.fpp) | `make` |
| New Fortran files added | `make clean && make` |
| configure.ac changes | `autoreconf -fi && ./configure ... && make` |

### Demo Notebook

```bash
uv run jupyter lab demo/demo.ipynb
```

## PyPI Package (genepie)

The package is distributed via PyPI as `genepie`. Users can install with:

```bash
pip install genepie
```

### Package Contents

| Component | Description |
|-----------|-------------|
| `genepie` Python package | Python interface with `libpython_interface.so` |
| `atdyn` | MD engine (CLI command) |
| 43 analysis tools | CLI commands (rmsd_analysis, trj_analysis, etc.) |

**Note**: `spdyn` requires MPI and is NOT included in PyPI. Users needing spdyn should build from source.

### Key Package Files

| File | Purpose |
|------|---------|
| `pyproject.toml` | Package metadata, dependencies, CLI entry points, package-data (test data) |
| `setup.py` | Forces platform-specific wheel (not pure Python) |
| `src/genepie/cli.py` | CLI entry point functions for all binaries |
| `MANIFEST.in` | Specifies files to include in source distribution (including test data) |

### CLI Commands Available After Installation

```bash
# MD Engine
atdyn INP

# Analysis Tools (examples)
rmsd_analysis INP
trj_analysis INP
mbar_analysis INP
wham_analysis INP
crd_convert INP
kmeans_clustering INP
# ... and 37 more
```

## GitHub Actions (CI/CD)

### Workflow: `.github/workflows/build-wheels.yml`

Multi-stage build process for PyPI wheel distribution:

| Stage | Purpose | Platforms |
|-------|---------|-----------|
| `build-genesis` | Build Fortran binaries (atdyn, analysis tools, libpython_interface) | Linux x86_64, macOS arm64, macOS x86_64 |
| `build-wheels` | Create Python wheels with bundled binaries | Same as above |
| `build-sdist` | Create source distribution | Linux |
| `publish` | Upload to PyPI (on tag push `v*`) | Linux |
| `publish-testpypi` | Upload to TestPyPI (on workflow_dispatch) | Linux |

### Triggers

- **Push tags `v*`**: Build and publish to PyPI
- **Pull requests to main**: Build and test (no publish)
- **workflow_dispatch**: Build and publish to TestPyPI

### Publishing to PyPI

```bash
# Create and push a version tag
git tag v0.1.0
git push origin v0.1.0
# GitHub Actions will automatically build and publish
```

### Testing on TestPyPI

Use "Run workflow" button on GitHub Actions page to trigger `publish-testpypi`.

### Workflow: `.github/workflows/tests.yml`

Automated test suite that runs on every push and pull request to main:

| Job | Platform | Tests |
|-----|----------|-------|
| `build-and-test-linux` | Linux x86_64 (manylinux_2_28 container) | All tests |
| `build-and-test-macos` | macOS arm64 (Apple Silicon) | All tests |

**Tests executed:**
- Basic tests: test_rmsd, test_rg, test_drms, test_crd_convert, test_avecrd
- Regression tests: test_trj, test_wham, test_mbar_1d, test_mbar_block, test_hb_atom, test_hb_snap, test_kmeans
- Error handling tests: test_error_handling
- ATDYN tests: test_atdyn

**Note**: test_msd and test_diffusion are skipped in CI due to high memory usage. They can be run locally.

## Python Interface Architecture

### Directory: `src/genepie/`

The Python interface uses ctypes to call Fortran functions compiled into `libpython_interface.so`.

### Core Components

| File | Purpose |
|------|---------|
| `genesis_exe.py` | Main API - Python functions calling GENESIS analysis tools |
| `libgenesis.py` | C type definitions for Fortran subroutine bindings |
| `libloader.py` | Loads `libpython_interface.so` shared library |
| `s_molecule.py` | `SMolecule` class - Python representation of molecular structure |
| `s_molecule_c.py` | `SMoleculeC` class - ctypes Structure for C interface |
| `s_trajectories.py` | `STrajectories` class - trajectory data handling |
| `s_trajectories_c.py` | `STrajectoriesC` class - ctypes Structure for C interface |
| `c2py_util.py` | C data → numpy array conversion utilities (with validation) |
| `py2c_util.py` | numpy array → C data conversion utilities |
| `ctrl_files.py` | Generates temporary control files for GENESIS functions |
| `exceptions.py` | Custom exception classes (GenesisError hierarchy) |
| `file_validators.py` | File path validation (existence, patterns, topology) |
| `param_validators.py` | Parameter validation (enums, ranges, MD parameters) |
| `validation.py` | Input validation utilities for sizes/pointers |
| `output_capture.py` | Context managers for stdout/stderr capture |

### Fortran Interface Files (.fpp)

| Pattern | Purpose |
|---------|---------|
| `*_c_mod.fpp` | C-callable wrappers for GENESIS Fortran routines |
| `*_analysis.fpp` | Analysis algorithm implementations |
| `conv_f_c_util.fpp` | Fortran ↔ C data conversion utilities |
| `s_molecule_c_mod.fpp` | s_molecule C structure allocation/deallocation |
| `s_trajectories_c_mod.fpp` | s_trajectories C structure handling |
| `atdyn_c_mod.fpp` | ATDYN MD/minimization C wrappers with state reset |

#### Timer Reset for Library Mode

The `src/lib/timers.fpp` module includes a `reset_timers()` subroutine for library mode usage. This resets accumulated timer state between sequential atdyn runs:

```fortran
! In atdyn_c_mod.fpp - called at start of each MD/minimization run
call reset_timers()
```

The `reset_atdyn_state_c()` function is also exposed to Python for explicit state cleanup.

### Available Analysis Functions (in `genesis_exe.py`)

- `crd_convert()` - Coordinate/trajectory conversion
- `trj_analysis()` - Distance, angle, dihedral analysis
- `rmsd_analysis()` - RMSD calculation
- `drms_analysis()` - Distance RMSD calculation
- `rg_analysis()` - Radius of gyration
- `msd_analysis()` - Mean squared displacement
- `diffusion_analysis()` - Diffusion coefficient calculation
- `hb_analysis()` - Hydrogen bond analysis
- `avecrd_analysis()` - Average coordinate calculation
- `wham_analysis()` - WHAM free energy analysis
- `mbar_analysis()` - MBAR free energy analysis
- `kmeans_clustering()` - K-means trajectory clustering

### ATDYN MD Engine Functions (in `genesis_exe.py`)

The Python interface also provides functions to run molecular dynamics and energy minimization:

- `run_atdyn_md()` - Run MD simulation, returns energies and final coordinates
- `run_atdyn_min()` - Run energy minimization, returns energies and minimized coordinates

#### Supported File Formats

| Format | Topology | Coordinates | Parameters |
|--------|----------|-------------|------------|
| AMBER | `prmtopfile` | `ambcrdfile` | (in prmtop) |
| GROMACS | `grotopfile` | `grocrdfile` | (in grotop) |
| CHARMM | `psffile` | `pdbfile`/`crdfile` | `parfile`, `strfile` |

#### Example Usage

```python
from genepie import genesis_exe

# AMBER format
energies, coords = genesis_exe.run_atdyn_md(
    prmtopfile="protein.prmtop",
    ambcrdfile="protein.inpcrd",
    nsteps=1000,
    timestep=0.002,
    eneout_period=100,
    ensemble="NVE",
)

# GROMACS format
energies, coords = genesis_exe.run_atdyn_md(
    grotopfile="protein.top",
    grocrdfile="protein.gro",
    nsteps=1000,
    timestep=0.002,
)

# CHARMM format with parameter files
energies, coords = genesis_exe.run_atdyn_md(
    psffile="protein.psf",
    pdbfile="protein.pdb",
    parfile=["par_all36_prot.prm", "par_all36_lipid.prm"],
    strfile=["toppar_water.str"],
    nsteps=1000,
)

# Energy minimization
energies, coords, converged, final_grad = genesis_exe.run_atdyn_min(
    prmtopfile="protein.prmtop",
    ambcrdfile="protein.inpcrd",
    method="SD",  # Steepest Descent
    nsteps=500,
)
```

#### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `nsteps` | Number of MD steps or minimization steps | Required |
| `timestep` | Time step in ps | 0.001 |
| `ensemble` | NVE, NVT, or NPT | "NVE" |
| `temperature` | Target temperature in K | 300.0 |
| `eneout_period` | Energy output frequency | 100 |
| `electrostatic` | PME or CUTOFF | "PME" |
| `pme_alpha` | PME Ewald coefficient | Auto |
| `pme_ngrid_x/y/z` | PME grid dimensions | Auto |
| `switchdist` | Switch distance in Å | 10.0 |
| `cutoffdist` | Cutoff distance in Å | 12.0 |
| `pairlistdist` | Pair list distance in Å | 13.5 |
| `rigid_bond` | Enable SHAKE | False |
| `shake_iteration` | Max SHAKE iterations | 500 |
| `shake_tolerance` | SHAKE convergence | 1.0e-10 |

#### Fortran State Management

When running multiple atdyn calls in the same Python process, the Fortran library maintains global state (timers, FFT plans, etc.). For sequential runs:

- **2-3 runs**: Generally work in the same process
- **6+ runs**: May cause segfaults due to accumulated Fortran state

For reliable testing with many sequential runs, use the subprocess-isolated variants:

```python
from genepie import genesis_exe

# Isolated functions run in separate subprocess (crash-safe)
result = genesis_exe.run_atdyn_md_isolated(
    prmtopfile="protein.prmtop",
    ambcrdfile="protein.inpcrd",
    nsteps=1000,
    timeout=300.0,  # Optional timeout in seconds
)
print(result.energies, result.final_coords)

result = genesis_exe.run_atdyn_min_isolated(
    prmtopfile="protein.prmtop",
    ambcrdfile="protein.inpcrd",
    nsteps=500,
)
print(result.energies, result.final_coords, result.converged)
```

### Integration with MDTraj and MDAnalysis

`SMolecule` has conversion methods:
- `to_mdtraj_topology()` / `from_mdtraj_topology()`
- `to_mdanalysis_universe()` / `from_mdanalysis_universe()`

### Error Handling

Custom exception hierarchy for better error diagnosis:

| Exception | Purpose |
|-----------|---------|
| `GenesisError` | Base exception for all GENESIS errors |
| `GenesisFortranError` | Errors from Fortran code (includes `code` and `stderr_output` attributes) |
| `GenesisFortranMemoryError` | Memory allocation/deallocation errors (code 100-199) |
| `GenesisFortranFileError` | File I/O errors (code 200-299) |
| `GenesisFortranValidationError` | Parameter validation errors from Fortran (code 300-399) |
| `GenesisFortranDataError` | Data inconsistency errors (code 400-499) |
| `GenesisFortranNotSupportedError` | Unsupported feature errors (code 500-599) |
| `GenesisValidationError` | Input validation errors (before Fortran call) |
| `GenesisMemoryError` | Memory/pointer operation errors |
| `GenesisOverflowError` | Integer overflow in size calculations |

#### Error Code Categories

| Category | Code Range | Examples |
|----------|------------|----------|
| MEMORY_ERROR | 100-199 | ERROR_ALLOC (101), ERROR_DEALLOC (102) |
| FILE_ERROR | 200-299 | ERROR_FILE_NOT_FOUND (201), ERROR_FILE_FORMAT (202) |
| VALIDATION_ERROR | 300-399 | ERROR_INVALID_PARAM (301), ERROR_MISSING_PARAM (302) |
| DATA_ERROR | 400-499 | ERROR_DATA_MISMATCH (401), ERROR_NO_DATA (402) |
| NOT_SUPPORTED_ERROR | 500-599 | ERROR_NOT_SUPPORTED (501), ERROR_DIMENSION (502) |
| INTERNAL_ERROR | 600-699 | ERROR_INTERNAL (601), ERROR_SYNTAX (602) |

Usage:
```python
from genepie import (
    GenesisError,
    GenesisFortranError,
    GenesisFortranNotSupportedError,
    GenesisValidationError,
    ErrorCode,
)

try:
    result = genesis_exe.rmsd_analysis(...)
except GenesisFortranNotSupportedError as e:
    print(f"Unsupported feature (code {e.code}): {e}")
except GenesisFortranError as e:
    print(f"Fortran error (code {e.code}): {e}")
    print(f"Fortran stderr: {e.stderr_output}")
    # Check specific error codes
    if e.code == ErrorCode.ERROR_ATOM_COUNT:
        print("Atom count mismatch")
except GenesisValidationError as e:
    print(f"Invalid input: {e}")
except GenesisError as e:
    print(f"GENESIS error: {e}")
```

## Testing

### Test Directory Structure

```
src/genepie/tests/
├── data/                          # Test data (included in PyPI package)
│   ├── bpti/                      # BPTI test system
│   │   ├── BPTI_ionize.pdb
│   │   ├── BPTI_ionize.psf
│   │   └── BPTI_run.dcd
│   ├── ralp_dppc/                 # RALP-DPPC test system
│   │   ├── RALP_DPPC_run.pdb
│   │   ├── RALP_DPPC.psf
│   │   └── RALP_DPPC_run.dcd
│   ├── chignolin/                 # Downloaded from Google Drive (.gitignore)
│   │   ├── chignolin.pdb
│   │   ├── chignolin.psf
│   │   └── chignolin.dcd
│   ├── molecule.pdb
│   └── msd.data
├── conftest.py                    # Path constants for test data
├── download_test_data.py          # Download chignolin data from Google Drive
├── test_integration.py            # Comprehensive integration tests (42 tests)
├── test_*.py                      # Individual analysis tests
└── all_run.sh                     # Test runner script

tests/regression_test/             # Reference data for regression tests (source repo only)
└── test_analysis/                 # Expected output values
```

**Test Categories:**
- **Basic tests** (test_rmsd, test_rg, etc.): Validate output ranges, work with PyPI package
- **Regression tests** (test_trj, test_wham, test_mbar_*, test_atdyn): Compare against reference values, require source repo
- **Integration tests** (test_integration): Comprehensive tests, require chignolin download

### Testing Strategy

| Change Type | Tests to Run |
|-------------|--------------|
| Python API changes | Basic tests + `test_integration` |
| Fortran interface | All tests including regression tests |
| Bug fixes | Relevant test + add new regression test |
| New analysis function | Create `test_<name>.py` + add to `all_run.sh` |

### Running Tests

```bash
# First, build the Fortran shared library
make

# Install genepie in editable mode (if not already done)
uv pip install -e .

# Run individual tests (basic validation)
python -m genepie.tests.test_rmsd
python -m genepie.tests.test_crd_convert
python -m genepie.tests.test_rg

# Run all basic tests (18 tests)
cd src/genepie/tests
./all_run.sh

# Regression tests (compare with reference values in tests/regression_test/)
python -m genepie.tests.test_trj
python -m genepie.tests.test_wham
python -m genepie.tests.test_mbar_1d
python -m genepie.tests.test_mbar_block
python -m genepie.tests.test_atdyn
```

**Note**: Basic tests (test_rmsd, test_rg, etc.) validate output ranges only. Regression tests compare against reference values and require the full source repository.

### Integration Tests

Comprehensive test script (42 tests) covering all major functionality. Requires chignolin test data:

```bash
# Download chignolin data from Google Drive (first time only)
python -m genepie.tests.download_test_data

# Run integration tests
python -m genepie.tests.test_integration
```

Tests include:
- SMolecule loading and manipulation
- crd_convert trajectory loading with selections
- Analysis functions (trj_analysis, rg_analysis, rmsd_analysis, etc.)
- Free energy (WHAM) analysis
- MDTraj/MDAnalysis integration
- scikit-learn/PyTorch integration
- Error handling and exception classes

### ATDYN MD Engine Tests

Tests for the atdyn Python interface covering AMBER, GROMACS, and CHARMM file formats:

```bash
python -m genepie.tests.test_atdyn
```

Currently includes 6 tests:
- `test_atdyn_min_glycam` - AMBER format minimization
- `test_atdyn_md_glycam` - AMBER format MD with PME
- `test_atdyn_md_bpti` - GROMACS format MD with PME
- `test_atdyn_md_jac_param27` - CHARMM format MD with PME
- `test_atdyn_md_dppc_nvt` - CHARMM format NVT with Langevin
- `test_atdyn_min_dppc` - CHARMM format minimization

**Note**: These tests use subprocess isolation to avoid Fortran global state issues. Each test runs in a separate Python subprocess to ensure clean library state.

### Error Handling Tests

```bash
python -m genepie.tests.test_error_handling
```

## Adding a New Analysis Function

1. Create Fortran wrapper in `src/analysis/interface/python_interface/*_c_mod.fpp` with `bind(C)` interface
2. Create analysis implementation in `*_analysis.fpp` (or reuse existing GENESIS code)
3. Add function signature to `src/genepie/libgenesis.py`
4. Create Python wrapper function in `src/genepie/genesis_exe.py`
5. Write regression test as `src/genepie/tests/test_<analysis_name>.py`
6. Add test to `src/genepie/tests/all_run.sh`
7. Update `src/analysis/interface/python_interface/Makefile.am` to include new `.fpp` files

## Key Patterns

### Calling Fortran from Python
```python
# In genesis_exe.py
mol_c = molecule.to_SMoleculeC()  # Convert Python → C structure
try:
    with suppress_stdout_capture_stderr() as captured:  # Suppress stdout, capture stderr
        LibGenesis().lib.some_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(result_c),
            # ... other args
        )
    result = c2py_util.conv_double_ndarray(result_c, size)  # C → numpy
except Exception as e:
    # captured.stderr contains Fortran error messages
    raise GenesisFortranError(str(e), stderr_output=captured.stderr)
finally:
    LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))
```

### Control File Generation
Analysis functions use temporary control files (GENESIS INP format):
```python
with tempfile.NamedTemporaryFile(dir=os.getcwd(), delete=True) as ctrl:
    ctrl_files.write_ctrl_output(ctrl, rmsfile="dummy.rms")
    ctrl_files.write_ctrl_selection(ctrl, selection_group, selection_mole_name)
    # ... call Fortran function with ctrl.name
```

### Zerocopy Memory Management

**Zerocopy は STrajectories (座標データ) と解析結果のみに適用。SMolecule は従来のコピー方式。**

#### Zerocopy 対象

| コンポーネント | 方式 | メモリ所有者 |
|---------------|------|-------------|
| STrajectories (座標) | ✅ Zerocopy | Python GC |
| 解析結果 (RMSD等) | ✅ Zerocopy | Python GC |
| SMolecule (分子構造) | ❌ コピー | Python/Fortran両方 |

#### STrajectories の Zerocopy パターン

```
Python: np.zeros() → .ctypes.data_as(c_void_p) → Fortran: C_F_POINTER()
        ↑ GC owns memory                         ↑ Creates alias (no copy)
```

```python
coords = np.zeros((n_frame, n_atom, 3), dtype=np.float64)
lib.analysis_c(coords.ctypes.data_as(ctypes.c_void_p), ...)
# coords now contains Fortran's output
```

#### SMolecule のメモリ管理 (コピー方式)

```
Python SMolecule (GC管理)
    ↓ to_SMoleculeC() でコピー
C SMoleculeC (Fortran allocate)
    ↓ 解析関数で使用
    ↓
deallocate_s_molecule_c() で明示解放 ← 必須！
```

```python
mol_c = molecule.to_SMoleculeC()  # Python → C コピー
try:
    lib.some_analysis(ctypes.byref(mol_c), ...)
finally:
    lib.deallocate_s_molecule_c(ctypes.byref(mol_c))  # 明示解放
```

**SMolecule が zerocopy できない理由:**
- 文字列フィールド (6件): Unicode → ASCII エンコード必須
- 整数型不一致 (15+件): numpy int64 → Fortran int32
- メモリレイアウト (16+件): row-major ↔ column-major

#### Array Layout

| Language | Layout | Example Access |
|----------|--------|----------------|
| Python | `(nframe, natom, 3)` C-order | `coords[frame, atom, 0]` |
| Fortran | `(3, natom, nframe)` F-order | `coords_f(1, atom+1, frame+1)` |

#### Zerocopy 対応状況

| Zerocopy ✅ | Legacy (copy-based) |
|-------------|---------------------|
| `crd_convert`, `rmsd_analysis`, `rg_analysis` | `msd_analysis`, `hb_analysis` |
| `trj_analysis`, `drms_analysis`, `diffusion_analysis` | `wham_analysis`, `mbar_analysis`, `avecrd_analysis`, `kmeans_clustering` |

#### Safety Rules

1. numpy array must outlive Fortran call
2. Never resize array while Fortran holds pointer
3. Use `np.ascontiguousarray()` for non-contiguous arrays

## Code Structure (Full Repository)

- `src/spdyn/` - Domain decomposition MD engine (MPI parallel)
- `src/atdyn/` - Atom decomposition MD engine
- `src/lib/` - Shared Fortran library
- `src/analysis/` - Analysis tools
  - `libana/` - Analysis library
  - `interface/python_interface/` - Fortran interface files (`.fpp`) for Python bindings
  - `trj_analysis/`, `free_energy/`, `mode_analysis/`, etc.
- `src/genepie/` - **Main Python package** (PyPI distribution)
  - `tests/` - Test files for genepie

## Known Issues

### Ubuntu Version Compatibility

**Supported Ubuntu versions:**

| Ubuntu Version | glibc | Status |
|---------------|-------|--------|
| 24.04 LTS | 2.39 | Supported |
| 22.04 LTS | 2.35 | Supported |
| 20.04 LTS | 2.31 | Supported |
| 18.04 LTS | 2.27 | Not supported (glibc too old) |

The PyPI wheels require glibc 2.28+ (manylinux_2_28). Ubuntu 18.04 users must build from source.

### Linux: `undefined symbol: _ZGVdN4v_cos` (Fixed)

**Status**: This issue has been fixed in the default build configuration.

**Root cause**: The `-ffast-math` compiler flag was triggering GCC's auto-vectorization, which generated calls to glibc's libmvec library (`_ZGVdN4v_cos` and similar symbols).

**Fix**: The `-ffast-math` flag has been removed from the default compiler flags in `configure.ac`. PyPI wheels built after this change will not have this issue.

**If you still encounter this issue** (e.g., with older wheel versions or custom builds with `-ffast-math`):

```bash
# Workaround for affected systems
ARCH_DIR=$(gcc -print-multiarch 2>/dev/null || echo "x86_64-linux-gnu")
export LD_PRELOAD=/lib/${ARCH_DIR}/libmvec.so.1:/lib/${ARCH_DIR}/libm.so.6${LD_PRELOAD:+:$LD_PRELOAD}
python your_script.py
```

### Sequential atdyn runs cause segfaults

**Symptom**: Running 6+ atdyn MD simulations in the same Python process causes segfaults.

**Cause**: Fortran global state (PME FFT plans, memory allocations) accumulates across runs.

**Workaround**: Use subprocess isolation for tests (see test_atdyn.py pattern).
