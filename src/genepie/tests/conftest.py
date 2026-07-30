"""Test configuration and shared fixtures for genepie tests."""
import pathlib

# Test data directory
TEST_DIR = pathlib.Path(__file__).parent
DATA_DIR = TEST_DIR / "data"

# BPTI system paths
BPTI_PDB = DATA_DIR / "bpti" / "BPTI_ionize.pdb"
BPTI_PSF = DATA_DIR / "bpti" / "BPTI_ionize.psf"
BPTI_DCD = DATA_DIR / "bpti" / "BPTI_run.dcd"

# RALP-DPPC system paths
RALP_PDB = DATA_DIR / "ralp_dppc" / "RALP_DPPC_run.pdb"
RALP_PSF = DATA_DIR / "ralp_dppc" / "RALP_DPPC.psf"
RALP_DCD = DATA_DIR / "ralp_dppc" / "RALP_DPPC_run.dcd"

# Chignolin system paths (for integration tests - downloaded from Google Drive)
CHIGNOLIN_PDB = DATA_DIR / "chignolin" / "chignolin.pdb"
CHIGNOLIN_PSF = DATA_DIR / "chignolin" / "chignolin.psf"
CHIGNOLIN_DCD = DATA_DIR / "chignolin" / "chignolin.dcd"

# T-REMD trialanine (Ala3) tutorial data (downloaded from Google Drive).
# Populated by ``python -m genepie.tests.download_tremd_data``; see that module.
TREMD_DIR = DATA_DIR / "remd_alat_tutorial"
TREMD_NREPLICA = 20
TREMD_PDB = TREMD_DIR / "trialanine.pdb"
TREMD_PSF = TREMD_DIR / "ala3.psf"
# Per-state (parameter-sorted) inputs; {} expands to the 1-based state index.
TREMD_POT_PATTERN = str(TREMD_DIR / "remd_paramID{}.pot")
TREMD_DCD_PATTERN = str(TREMD_DIR / "remd_paramID{}_trialanine.dcd")
# Reference MBAR results produced in the original paper.
TREMD_REF_DIR = TREMD_DIR / "reference"
TREMD_FENE_REF = TREMD_REF_DIR / "fene.dat"
TREMD_WEIGHT_PATTERN = str(TREMD_REF_DIR / "weight{}.dat")
TREMD_TOR_PATTERN = str(TREMD_REF_DIR / "remd_paramID{}.tor")
# Temperature ladder (K) of the 20 parameter-sorted states.
TREMD_TEMPERATURES = [
    300.00, 302.53, 305.09, 307.65, 310.24,
    312.85, 315.47, 318.12, 320.78, 323.46,
    326.16, 328.87, 331.61, 334.37, 337.14,
    339.94, 342.75, 345.59, 348.44, 351.26,
]
TREMD_TARGET_TEMPERATURE = 300.0

# Other test data
MOLECULE_PDB = DATA_DIR / "molecule.pdb"
MSD_DATA = DATA_DIR / "msd.data"
