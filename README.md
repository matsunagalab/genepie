# GENESIS

**GENeralized-Ensemble SImulation System** - A molecular dynamics simulation software for biomolecular systems.

## Python Interface (genepie)

The `genepie` package provides a Python interface to GENESIS analysis tools and the ATDYN MD engine.

---

## For Users

### Installation

```bash
# Create virtual environment with uv
uv venv --python=python3.11
source .venv/bin/activate

# From PyPI (coming soon)
uv pip install genepie

# Currently available from TestPyPI:
uv pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ genepie
```

**Requirements:**
- Python 3.9+
- Linux (x86_64) or macOS (arm64, x86_64)
- glibc 2.28+ for Linux (see table below)

| Ubuntu Version | glibc | Status |
|---------------|-------|--------|
| 24.04 LTS | 2.39 | Supported |
| 22.04 LTS | 2.35 | Supported |
| 20.04 LTS | 2.31 | Supported |
| 18.04 LTS | 2.27 | Not supported (build from source) |

### Testing Your Installation

```bash
# Run individual tests (no additional data required)
python -m genepie.tests.test_rmsd
python -m genepie.tests.test_crd_convert
python -m genepie.tests.test_rg
python -m genepie.tests.test_drms
python -m genepie.tests.test_avecrd

# Integration tests (requires ~500 MB download)
uv pip install gdown mdtraj MDAnalysis
python -m genepie.tests.download_test_data
python -m genepie.tests.test_integration
```

Note: Some tests (test_trj, test_wham, test_mbar_*, test_atdyn) require the full source repository and are intended for developers only.

### Quick Start

```python
from genepie import genesis_exe, SMolecule

# Load molecular structure
mol = SMolecule.from_file(pdbfile="protein.pdb", psffile="protein.psf")
print(f"Loaded {mol.num_atoms} atoms")

# Load trajectory and calculate RMSD
traj = genesis_exe.crd_convert(
    psffile="protein.psf",
    pdbfile="protein.pdb",
    dcdfile="trajectory.dcd",
    selection_group="an:CA",
)
rmsd = genesis_exe.rmsd_analysis(molecule=mol, trajectories=traj)
print(f"RMSD: {rmsd.mean():.2f} Å")

# Run MD simulation
energies, coords = genesis_exe.run_atdyn_md(
    prmtopfile="protein.prmtop",
    ambcrdfile="protein.inpcrd",
    nsteps=1000,
    ensemble="NVT",
    temperature=300.0,
)
```

### Available Analysis Functions

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

### MD Engine Functions

- `run_atdyn_md()` - Run MD simulation
- `run_atdyn_min()` - Run energy minimization
- `run_atdyn_md_isolated()` - Run MD in subprocess (crash-safe)
- `run_atdyn_min_isolated()` - Run minimization in subprocess

### Supported File Formats

| Format | Topology | Coordinates | Parameters |
|--------|----------|-------------|------------|
| AMBER | `prmtopfile` | `ambcrdfile` | (in prmtop) |
| GROMACS | `grotopfile` | `grocrdfile` | (in grotop) |
| CHARMM | `psffile` | `pdbfile`/`crdfile` | `parfile`, `strfile` |

---

## For Developers

### Installation from Source

```bash
# Clone repository
git clone https://github.com/matsunagalab/genepie.git
cd genepie

# Set up Python environment with uv
uv venv --python=python3.11
source .venv/bin/activate
uv pip install numpy

# Install build dependencies
# Linux (Ubuntu/Debian):
#   sudo apt install gfortran liblapack-dev libblas-dev autoconf automake libtool
# macOS:
#   brew install gcc lapack autoconf automake libtool

# Build GENESIS
autoreconf -fi

# Linux:
./configure --disable-mpi CC=gcc FC=gfortran LAPACK_LIBS="-llapack -lblas"

# macOS:
./configure --disable-mpi CC=gcc-14 FC=gfortran \
    LAPACK_LIBS="-L$(brew --prefix lapack)/lib -llapack -lblas"

make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)

# Install in editable mode
uv pip install -e .
```

### Running Tests

```bash
# Run individual tests
python -m genepie.tests.test_rmsd
python -m genepie.tests.test_crd_convert

# Run all basic tests (18 tests)
cd src/genepie/tests
./all_run.sh

# Integration tests (requires ~500 MB download)
uv pip install gdown mdtraj MDAnalysis
python -m genepie.tests.download_test_data
python -m genepie.tests.test_integration    # 42 tests

# Error handling tests
python -m genepie.tests.test_error_handling # 64 tests

# Regression tests (compare with reference values)
# These tests use data in tests/regression_test/
python -m genepie.tests.test_trj
python -m genepie.tests.test_wham
python -m genepie.tests.test_mbar_1d
python -m genepie.tests.test_mbar_block
python -m genepie.tests.test_atdyn
```

### Developer Workflow

#### 1. Local Development Cycle

| Step | Action | Notes |
|------|--------|-------|
| 1 | Edit code | Python: instant; Fortran: requires `make` |
| 2 | Run tests | See testing strategy below |
| 3 | Commit changes | Use descriptive commit messages |
| 4 | Push to branch | Create feature branch from main |
| 5 | Create PR | Target: main branch |

#### 2. Rebuilding After Changes

| Change Type | Required Action |
|-------------|-----------------|
| Python files (.py) | No rebuild needed |
| Fortran files (.fpp) | `make` |
| New Fortran files | `make clean && make` |
| configure.ac changes | `autoreconf -fi && ./configure ... && make` |

#### 3. Testing Strategy

| Change Type | Tests to Run |
|-------------|--------------|
| Python API changes | Basic tests + `test_integration` |
| Fortran interface | All tests including regression tests |
| Bug fixes | Relevant test + add new regression test |
| New analysis function | Create `test_<name>.py` + add to `all_run.sh` |

#### 4. Pull Request Process

1. Push changes to a feature branch
2. Create PR targeting `main` branch
3. GitHub Actions automatically builds wheels for all platforms
4. Review CI build results (Linux x86_64, macOS arm64, macOS x86_64)
5. Merge after approval and CI passes

#### 5. Release Workflow

```
Feature Branch → PR → main → TestPyPI (optional) → Tag → PyPI
```

**Step 1: Test on TestPyPI (optional but recommended)**

```bash
# Trigger manually via GitHub Actions UI:
# Actions → "Build and publish wheels" → Run workflow (workflow_dispatch)

# Then test the package:
pip install --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple/ genepie
python -c "import genepie; print(genepie.__version__)"
```

**Step 2: Release to PyPI**

```bash
# 1. Update version in pyproject.toml
# 2. Commit the version bump
git add pyproject.toml
git commit -m "Bump version to 0.1.x"
git push origin main

# 3. Create and push version tag
git tag v0.1.x
git push origin v0.1.x
# GitHub Actions automatically builds and publishes to PyPI
```

#### 6. CI/CD Pipeline Overview

| Workflow | Trigger | Action |
|----------|---------|--------|
| `tests.yml` | Push/PR to main | Run test suite on Linux & macOS |
| `build-wheels.yml` | Push/PR to main | Build wheels for all platforms |
| `build-wheels.yml` | `workflow_dispatch` | Build + publish to TestPyPI |
| `build-wheels.yml` | Push tag `v*` | Build + publish to PyPI |

**Test Platforms (`tests.yml`):**
- Linux x86_64 (manylinux_2_28 container)
- macOS arm64 (Apple Silicon)

**Build Platforms (`build-wheels.yml`):**
- Linux x86_64 (manylinux_2_28)
- macOS arm64 (Apple Silicon)
- macOS x86_64 (Intel)

### Project Structure

```
genepie/
├── src/
│   ├── genepie/           # Python interface (main package)
│   │   ├── genesis_exe.py # Analysis function wrappers
│   │   ├── libloader.py   # Shared library loader
│   │   └── tests/         # Test files and data
│   ├── atdyn/             # MD engine
│   └── analysis/          # Analysis tools
├── CLAUDE.md              # Developer guide for Claude Code
└── pyproject.toml         # Package configuration
```

### Adding a New Analysis Tool

**Recommended: Unified Architecture** (for RMSD, RG, DRMS pattern)

When CLI already exists with `trj_source_mod` support:
1. Export `analyze_*_unified()` from CLI's `*_analyze.fpp` with primitive arguments
2. Create `*_c_mod.fpp` that calls unified function via `init_source_memory()` + `init_sink_array()`
3. Add function signature to `libgenesis.py` and wrapper to `genesis_exe.py`
4. Add test as `tests/test_<name>.py`

**Alternative: Separate Implementation** (for HB, WHAM, MBAR pattern)

When unified pattern doesn't fit:
1. Create `*_c_mod.fpp` with `bind(C)` interface
2. Create `*_impl.fpp` for analysis implementation
3. Update `Makefile.am` to include new `.fpp` files
4. Add function signature to `libgenesis.py` and wrapper to `genesis_exe.py`
5. Add test as `tests/test_<name>.py`

See [CLAUDE.md](CLAUDE.md) for detailed instructions.

---

## Documentation

- [GENESIS Website](https://www.r-ccs.riken.jp/labs/cbrt/)
- [CLAUDE.md](CLAUDE.md) - Developer guide

## License

LGPL-3.0-or-later. See LICENSE file for details.
