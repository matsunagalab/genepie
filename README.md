# GENESIS

**GENeralized-Ensemble SImulation System** - A molecular dynamics simulation software for biomolecular systems.

## Python Interface (genepie)

The `genepie` package provides a Python interface to GENESIS analysis tools and the ATDYN MD engine.

---

## For Users

### Installation

```bash
# From PyPI (coming soon)
pip install genepie

# Currently available from TestPyPI:
pip install -i https://test.pypi.org/simple/ genepie
```

**Requirements:**
- Python 3.9+
- Linux (x86_64) or macOS (arm64, x86_64)
- glibc 2.28+ for Linux (Ubuntu 20.04+)

### Testing Your Installation

```bash
# Run individual tests
python -m genepie.tests.test_rmsd
python -m genepie.tests.test_crd_convert
python -m genepie.tests.test_trj

# Run all basic tests (18 tests)
cd $(python -c "import genepie; from pathlib import Path; print(Path(genepie.__file__).parent / 'tests')")
./all_run.sh

# Integration tests (requires ~500 MB download)
pip install gdown mdtraj MDAnalysis
python -m genepie.tests.download_test_data
python -m genepie.tests.test_integration
```

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
git clone https://github.com/matsunagalab/genesis.git
cd genesis

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
python -m genepie.tests.test_wham

# Run all basic tests (18 tests)
cd src/genepie/tests
./all_run.sh

# Integration tests (requires ~500 MB download)
uv pip install gdown mdtraj MDAnalysis
python -m genepie.tests.download_test_data
python -m genepie.tests.test_integration    # 42 tests

# Error handling tests
python -m genepie.tests.test_error_handling # 64 tests
```

### Project Structure

```
genesis/
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

**1. Fortran side** (`src/analysis/interface/python_interface/`):
- Create `*_c_mod.fpp` with `bind(C)` interface wrapping GENESIS routines
- Create `*_impl.fpp` for analysis implementation (or reuse existing code)
- Update `Makefile.am` to include new `.fpp` files

**2. Python side** (`src/genepie/`):
- Add function signature to `libgenesis.py`
- Create wrapper function in `genesis_exe.py`
- Add test as `tests/test_<name>.py`

See [CLAUDE.md](CLAUDE.md) for detailed instructions.

---

## Documentation

- [GENESIS Website](https://www.r-ccs.riken.jp/labs/cbrt/)
- [CLAUDE.md](CLAUDE.md) - Developer guide

## License

LGPL-3.0-or-later. See LICENSE file for details.
