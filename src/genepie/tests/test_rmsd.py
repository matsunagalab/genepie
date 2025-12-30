# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import os
import subprocess
import sys
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..ctrl_files import TrajectoryParameters
from ..s_molecule import SMolecule
from .. import genesis_exe


def test_rmsd_no_fitting():
    """Test RMSD analysis without fitting (direct comparison to reference)."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
    trajs, subset_mol = genesis_exe.crd_convert(
            mol,
            traj_params=[
                TrajectoryParameters(
                    trjfile=str(BPTI_DCD),
                    md_step=10,
                    mdout_period=1,
                    ana_period=1,
                    repeat=1,
                ),
            ],
            trj_format="DCD",
            trj_type="COOR+BOX",
            trj_natom=0,
            selection_group=["all"],
            fitting_method="NO",
            fitting_atom=1,
            check_only=False,
            pbc_correct="NO",
    )
    _ = subset_mol

    try:
        for t in trajs:
            # Test RMSD without fitting (fitting_selection=None)
            result = genesis_exe.rmsd_analysis(
                mol, t,
                analysis_selection="an:CA",
                fitting_selection=None,  # No fitting
                ana_period=1,
                mass_weighted=False,
            )

            # Validate results
            assert result.rmsd is not None, "RMSD result should not be None"
            assert len(result.rmsd) > 0, "RMSD result should have at least one frame"
            assert all(r >= 0 for r in result.rmsd), "RMSD values should be non-negative"
            assert all(r < 50.0 for r in result.rmsd), "RMSD values should be reasonable (< 50 A)"

            print(f"RMSD no fitting (n={len(result.rmsd)}): "
                  f"min={min(result.rmsd):.5f}, max={max(result.rmsd):.5f}")

    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def test_rmsd_with_fitting():
    """Test RMSD analysis with TR+ROT fitting."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
    trajs, subset_mol = genesis_exe.crd_convert(
            mol,
            traj_params=[
                TrajectoryParameters(
                    trjfile=str(BPTI_DCD),
                    md_step=10,
                    mdout_period=1,
                    ana_period=1,
                    repeat=1,
                ),
            ],
            trj_format="DCD",
            trj_type="COOR+BOX",
            trj_natom=0,
            selection_group=["all"],
            fitting_method="NO",
            fitting_atom=1,
            check_only=False,
            pbc_correct="NO",
    )
    _ = subset_mol

    try:
        for t in trajs:
            # Test RMSD with TR+ROT fitting
            result = genesis_exe.rmsd_analysis(
                mol, t,
                analysis_selection="sid:BPTI and an:CA",
                fitting_selection="sid:BPTI and an:CA",  # Enable fitting
                fitting_method="TR+ROT",
                ana_period=1,
                mass_weighted=False,
            )

            # Validate results
            assert result.rmsd is not None, "RMSD with fitting should not be None"
            assert len(result.rmsd) > 0, "RMSD should have at least one frame"
            assert all(r >= 0 for r in result.rmsd), "RMSD values should be non-negative"
            assert all(r < 50.0 for r in result.rmsd), "RMSD values should be reasonable (< 50 A)"

            print(f"RMSD with fitting (n={len(result.rmsd)}): "
                  f"min={min(result.rmsd):.5f}, max={max(result.rmsd):.5f}")

    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def _run_test_in_subprocess(test_name: str) -> bool:
    """Run a single test function in isolated subprocess to avoid Fortran state issues."""
    code = f'''
import sys
# Handle package imports when run as subprocess
if __name__ == "__main__":
    import pathlib
    pkg_dir = pathlib.Path("{__file__}").resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))

from genepie.tests.test_rmsd import {test_name}
{test_name}()
'''
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        timeout=120
    )

    if result.returncode == 0:
        # Print stdout (test output)
        if result.stdout:
            print(result.stdout, end='')
        return True
    else:
        print(f"stderr: {result.stderr}" if result.stderr else "")
        return False


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")

    # Run each test in separate subprocess to avoid Fortran global state accumulation
    # This prevents crashes on macOS ARM64 when multiple tests run in same process
    tests = [
        "test_rmsd_no_fitting",
        "test_rmsd_with_fitting",
    ]

    failed = []
    for test_name in tests:
        if _run_test_in_subprocess(test_name):
            print(f"\n{test_name}: PASSED")
        else:
            print(f"\n{test_name}: FAILED")
            failed.append(test_name)

    if failed:
        raise RuntimeError(f"Tests failed: {', '.join(failed)}")

    print("\nAll RMSD tests passed!")


if __name__ == "__main__":
    main()
