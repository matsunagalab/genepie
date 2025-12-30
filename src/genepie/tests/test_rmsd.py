# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import os
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..ctrl_files import TrajectoryParameters
from ..s_molecule import SMolecule
from .. import genesis_exe


def test_rmsd_analysis():
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)
    trajs, subset_mol =  genesis_exe.crd_convert(
            mol,
            traj_params = [
                TrajectoryParameters(
                    trjfile = str(BPTI_DCD),
                    md_step = 10,
                    mdout_period = 1,
                    ana_period = 1,
                    repeat = 1,
                    ),
                ],
            trj_format = "DCD",
            trj_type = "COOR+BOX",
            trj_natom = 0,
            selection_group = ["all", ],
            fitting_method = "NO",
            fitting_atom = 1,
            check_only = False,
            pbc_correct = "NO",
    )

    _ = subset_mol

    try:
        for t in trajs:
            d = genesis_exe.rmsd_analysis(
                    mol, t,
                    selection_group = ["sid:BPTI and an:CA", ],
                    fitting_method = "TR+ROT",
                    fitting_atom = 1,
                    check_only = False,
                    analysis_atom  = 1,
                    )
            # Validate RMSD results
            assert d.rmsd is not None, "RMSD result should not be None"
            assert len(d.rmsd) > 0, "RMSD result should have at least one frame"
            assert all(r >= 0 for r in d.rmsd), "RMSD values should be non-negative"
            # RMSD should be in reasonable range (0-50 Angstroms for proteins)
            assert all(r < 50.0 for r in d.rmsd), "RMSD values should be reasonable (< 50 Å)"
            print(f"RMSD values (n={len(d.rmsd)}): min={min(d.rmsd):.3f}, max={max(d.rmsd):.3f}", flush=True)

    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def test_rmsd_zerocopy():
    """Test RMSD zerocopy analysis (no fitting version)."""
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
            # Test zerocopy RMSD (no fitting)
            result = genesis_exe.rmsd_analysis_zerocopy(
                mol, t,
                analysis_selection="an:CA",
                ana_period=1,
                mass_weighted=False,
            )

            # Validate results
            assert result.rmsd is not None, "RMSD result should not be None"
            assert len(result.rmsd) > 0, "RMSD result should have at least one frame"
            assert all(r >= 0 for r in result.rmsd), "RMSD values should be non-negative"
            assert all(r < 50.0 for r in result.rmsd), "RMSD values should be reasonable (< 50 Å)"

            print(f"Zerocopy RMSD (n={len(result.rmsd)}): "
                  f"min={min(result.rmsd):.5f}, max={max(result.rmsd):.5f}")

    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def test_rmsd_zerocopy_vs_legacy_no_fitting():
    """Compare zerocopy and legacy RMSD when fitting is disabled."""
    import numpy as np

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
            # Run legacy version with NO fitting
            result_legacy = genesis_exe.rmsd_analysis(
                mol, t,
                selection_group=["sid:BPTI and an:CA"],
                fitting_method="NO",
                fitting_atom=1,
                check_only=False,
                analysis_atom=1,
            )

            # Run zerocopy version
            result_zerocopy = genesis_exe.rmsd_analysis_zerocopy(
                mol, t,
                analysis_selection="sid:BPTI and an:CA",
                ana_period=1,
                mass_weighted=False,
            )

            print(f"Legacy RMSD (no fitting): min={min(result_legacy.rmsd):.5f}, "
                  f"max={max(result_legacy.rmsd):.5f}")
            print(f"Zerocopy RMSD (no fitting): min={min(result_zerocopy.rmsd):.5f}, "
                  f"max={max(result_zerocopy.rmsd):.5f}")

            # Check lengths match
            assert len(result_legacy.rmsd) == len(result_zerocopy.rmsd), \
                "Result lengths should match"

            # Check values are close (within 0.001 Angstrom)
            diff = np.abs(np.array(result_legacy.rmsd) - np.array(result_zerocopy.rmsd))
            max_diff = np.max(diff)
            print(f"Max difference: {max_diff:.6f} Å")
            assert max_diff < 0.001, f"Results should match within 0.001 Å, got {max_diff}"

    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def test_rmsd_zerocopy_with_fitting():
    """Test RMSD zerocopy with TR+ROT fitting vs legacy implementation."""
    import numpy as np

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
            # Test zerocopy with fitting (TR+ROT)
            result_zc = genesis_exe.rmsd_analysis_zerocopy_with_fitting(
                mol, t,
                fitting_selection="sid:BPTI and an:CA",
                analysis_selection="sid:BPTI and an:CA",
                fitting_method="TR+ROT",
                ana_period=1,
                mass_weighted=False,
            )

            # Validate results
            assert result_zc.rmsd is not None, "Zerocopy fitting RMSD should not be None"
            assert len(result_zc.rmsd) > 0, "Zerocopy fitting RMSD should have values"
            assert all(r >= 0 for r in result_zc.rmsd), "RMSD values should be non-negative"
            assert all(r < 50.0 for r in result_zc.rmsd), "RMSD values should be reasonable (< 50 Å)"

            print(f"Zerocopy with fitting RMSD (n={len(result_zc.rmsd)}): "
                  f"min={min(result_zc.rmsd):.5f}, max={max(result_zc.rmsd):.5f}")

            # Compare with legacy implementation
            result_legacy = genesis_exe.rmsd_analysis(
                mol, t,
                selection_group=["sid:BPTI and an:CA"],
                fitting_method="TR+ROT",
                fitting_atom=1,
                check_only=False,
                analysis_atom=1,
            )

            print(f"Legacy RMSD (n={len(result_legacy.rmsd)}): "
                  f"min={min(result_legacy.rmsd):.5f}, max={max(result_legacy.rmsd):.5f}")

            # Check that results match
            assert len(result_zc.rmsd) == len(result_legacy.rmsd), \
                f"Result lengths should match: {len(result_zc.rmsd)} vs {len(result_legacy.rmsd)}"

            # Compare values (should be very close)
            diff = np.abs(np.array(result_zc.rmsd) - np.array(result_legacy.rmsd))
            max_diff = np.max(diff)
            print(f"Max difference between zerocopy and legacy: {max_diff:.6f} Å")

            # Allow small numerical tolerance
            assert max_diff < 0.001, f"Results should match within 0.001 Å, got {max_diff}"

    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def _run_test_in_subprocess(test_name: str) -> bool:
    """Run a single test function in isolated subprocess to avoid Fortran state issues."""
    import subprocess
    import sys

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
        "test_rmsd_zerocopy",
        "test_rmsd_zerocopy_with_fitting",
        "test_rmsd_analysis",
    ]

    failed = []
    for test_name in tests:
        if _run_test_in_subprocess(test_name):
            print(f"\n✓ {test_name}: PASSED")
        else:
            print(f"\n✗ {test_name}: FAILED")
            failed.append(test_name)

    if failed:
        raise RuntimeError(f"Tests failed: {', '.join(failed)}")


if __name__ == "__main__":
    main()
