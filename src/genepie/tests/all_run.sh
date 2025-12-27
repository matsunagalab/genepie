#!/bin/bash

# Change to the script's directory (genepie/tests/)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Add src/ directory to PYTHONPATH so python can find the genepie module
export PYTHONPATH="${SCRIPT_DIR}/../..:${PYTHONPATH}"

python -m genepie.tests.test_crd_convert "$@"
python -m genepie.tests.test_trj "$@"
python -m genepie.tests.test_wham "$@"
python -m genepie.tests.test_mbar_1d "$@"
python -m genepie.tests.test_mbar_block "$@"
python -m genepie.tests.test_avecrd "$@"
python -m genepie.tests.test_kmeans "$@"
python -m genepie.tests.test_hb_atom "$@"
python -m genepie.tests.test_hb_snap "$@"
python -m genepie.tests.test_rmsd "$@"
python -m genepie.tests.test_drms "$@"
python -m genepie.tests.test_rg "$@"
python -m genepie.tests.test_msd "$@"
python -m genepie.tests.test_diffusion "$@"
python -m genepie.tests.test_mdanalysis "$@"
python -m genepie.tests.test_mdtraj "$@"
python -m genepie.tests.test_error_handling "$@"
# Note: test_atdyn uses subprocess isolation to avoid Fortran global state issues
python -m genepie.tests.test_atdyn "$@"
