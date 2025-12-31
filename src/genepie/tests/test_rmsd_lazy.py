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
import numpy as np
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..s_molecule import SMolecule
from .. import genesis_exe


def test_rmsd_lazy_no_fitting():
    """Test lazy RMSD analysis without fitting."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)

    # Run lazy RMSD analysis
    result = genesis_exe.rmsd_analysis_lazy(
        mol,
        str(BPTI_DCD),
        analysis_selection="an:CA",
        fitting_selection=None,  # No fitting
        ana_period=1,
        mass_weighted=False,
        has_box=True,  # BPTI_DCD has box info
    )

    # Validate DCD info
    assert result.dcd_nframe > 0, "DCD should have at least one frame"
    assert result.dcd_natom > 0, "DCD should have atoms"
    assert result.dcd_natom == mol.num_atoms, "DCD atom count should match molecule"

    # Validate RMSD results
    assert result.rmsd is not None, "RMSD result should not be None"
    assert len(result.rmsd) > 0, "RMSD result should have at least one frame"
    assert all(r >= 0 for r in result.rmsd), "RMSD values should be non-negative"
    assert all(r < 50.0 for r in result.rmsd), "RMSD values should be reasonable (< 50 A)"

    print(f"Lazy RMSD no fitting (n={len(result.rmsd)}): "
          f"min={min(result.rmsd):.5f}, max={max(result.rmsd):.5f}")
    print(f"DCD info: nframe={result.dcd_nframe}, natom={result.dcd_natom}")


def test_rmsd_lazy_with_fitting():
    """Test lazy RMSD analysis with TR+ROT fitting."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)

    # Run lazy RMSD analysis with fitting
    result = genesis_exe.rmsd_analysis_lazy(
        mol,
        str(BPTI_DCD),
        analysis_selection="sid:BPTI and an:CA",
        fitting_selection="sid:BPTI and an:CA",
        fitting_method="TR+ROT",
        ana_period=1,
        mass_weighted=False,
        has_box=True,
    )

    # Validate results
    assert result.rmsd is not None, "RMSD with fitting should not be None"
    assert len(result.rmsd) > 0, "RMSD should have at least one frame"
    assert all(r >= 0 for r in result.rmsd), "RMSD values should be non-negative"
    assert all(r < 50.0 for r in result.rmsd), "RMSD values should be reasonable (< 50 A)"

    print(f"Lazy RMSD with fitting (n={len(result.rmsd)}): "
          f"min={min(result.rmsd):.5f}, max={max(result.rmsd):.5f}")


def _run_memory_rmsd():
    """Helper function to run memory-based RMSD analysis (for subprocess)."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
    trajs, _ = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
    )
    for t in trajs:
        result = genesis_exe.rmsd_analysis(
            mol, t,
            analysis_selection="sid:BPTI and an:CA",
            fitting_selection="sid:BPTI and an:CA",
            fitting_method="TR+ROT",
            ana_period=1,
            mass_weighted=False,
        )
    # Print as JSON for parsing
    import json
    print("RMSD_RESULT:" + json.dumps(list(result.rmsd)))


def _run_lazy_rmsd():
    """Helper function to run lazy-based RMSD analysis (for subprocess)."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
    result = genesis_exe.rmsd_analysis_lazy(
        mol,
        str(BPTI_DCD),
        analysis_selection="sid:BPTI and an:CA",
        fitting_selection="sid:BPTI and an:CA",
        fitting_method="TR+ROT",
        ana_period=1,
        mass_weighted=False,
        has_box=True,
    )
    # Print as JSON for parsing
    import json
    print("RMSD_RESULT:" + json.dumps(list(result.rmsd)))


def test_rmsd_lazy_vs_memory():
    """Compare lazy RMSD analysis with memory-based RMSD analysis.

    Both methods should produce identical results.
    Runs each method in separate subprocess to avoid Fortran state issues.
    """
    import json

    # Run memory-based in subprocess
    memory_code = f'''
import sys, pathlib
pkg_dir = pathlib.Path("{__file__}").resolve().parent
sys.path.insert(0, str(pkg_dir.parent.parent))
from genepie.tests.test_rmsd_lazy import _run_memory_rmsd
_run_memory_rmsd()
'''
    memory_proc = subprocess.run(
        [sys.executable, "-c", memory_code],
        capture_output=True, text=True, timeout=120
    )
    if memory_proc.returncode != 0:
        print(f"Memory stderr: {memory_proc.stderr}")
        raise RuntimeError(f"Memory-based RMSD failed: {memory_proc.stderr}")

    memory_result = None
    for line in memory_proc.stdout.split('\n'):
        if line.startswith("RMSD_RESULT:"):
            memory_result = json.loads(line[len("RMSD_RESULT:"):])
            break
    if memory_result is None:
        raise RuntimeError(f"No RMSD_RESULT in memory output: {memory_proc.stdout}")

    # Run lazy-based in subprocess
    lazy_code = f'''
import sys, pathlib
pkg_dir = pathlib.Path("{__file__}").resolve().parent
sys.path.insert(0, str(pkg_dir.parent.parent))
from genepie.tests.test_rmsd_lazy import _run_lazy_rmsd
_run_lazy_rmsd()
'''
    lazy_proc = subprocess.run(
        [sys.executable, "-c", lazy_code],
        capture_output=True, text=True, timeout=120
    )
    if lazy_proc.returncode != 0:
        print(f"Lazy stderr: {lazy_proc.stderr}")
        raise RuntimeError(f"Lazy RMSD failed: {lazy_proc.stderr}")

    lazy_result = None
    for line in lazy_proc.stdout.split('\n'):
        if line.startswith("RMSD_RESULT:"):
            lazy_result = json.loads(line[len("RMSD_RESULT:"):])
            break
    if lazy_result is None:
        raise RuntimeError(f"No RMSD_RESULT in lazy output: {lazy_proc.stdout}")

    # Compare results
    assert len(memory_result) == len(lazy_result), \
        f"Frame count mismatch: memory={len(memory_result)}, lazy={len(lazy_result)}"

    # Check values are close (allowing small floating point differences)
    for i, (mem_val, lazy_val) in enumerate(zip(memory_result, lazy_result)):
        assert np.isclose(mem_val, lazy_val, rtol=1e-4, atol=1e-6), \
            f"Frame {i}: memory={mem_val}, lazy={lazy_val}"

    print(f"Memory vs Lazy comparison passed: {len(memory_result)} frames")
    print(f"  Memory: min={min(memory_result):.5f}, max={max(memory_result):.5f}")
    print(f"  Lazy:   min={min(lazy_result):.5f}, max={max(lazy_result):.5f}")


def test_rmsd_lazy_ana_period():
    """Test lazy RMSD analysis with ana_period > 1."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)

    # Run with ana_period=2
    result_p2 = genesis_exe.rmsd_analysis_lazy(
        mol,
        str(BPTI_DCD),
        analysis_selection="an:CA",
        fitting_selection=None,
        ana_period=2,
        mass_weighted=False,
        has_box=True,
    )

    # Run with ana_period=1
    result_p1 = genesis_exe.rmsd_analysis_lazy(
        mol,
        str(BPTI_DCD),
        analysis_selection="an:CA",
        fitting_selection=None,
        ana_period=1,
        mass_weighted=False,
        has_box=True,
    )

    # With ana_period=2, we should have about half the frames
    expected_frames = result_p1.dcd_nframe // 2
    assert len(result_p2.rmsd) == expected_frames, \
        f"Expected {expected_frames} frames with ana_period=2, got {len(result_p2.rmsd)}"

    print(f"ana_period test passed:")
    print(f"  ana_period=1: {len(result_p1.rmsd)} frames")
    print(f"  ana_period=2: {len(result_p2.rmsd)} frames")


def _run_test_in_subprocess(test_name: str) -> bool:
    """Run a single test function in isolated subprocess to avoid Fortran state issues."""
    code = f'''
import sys
# Handle package imports when run as subprocess
if __name__ == "__main__":
    import pathlib
    pkg_dir = pathlib.Path("{__file__}").resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))

from genepie.tests.test_rmsd_lazy import {test_name}
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
        print(f"stdout: {result.stdout}" if result.stdout else "")
        print(f"stderr: {result.stderr}" if result.stderr else "")
        return False


def main():
    # Run each test in separate subprocess to avoid Fortran global state accumulation
    tests = [
        "test_rmsd_lazy_no_fitting",
        "test_rmsd_lazy_with_fitting",
        "test_rmsd_lazy_vs_memory",
        "test_rmsd_lazy_ana_period",
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

    print("\nAll lazy RMSD tests passed!")


if __name__ == "__main__":
    main()
