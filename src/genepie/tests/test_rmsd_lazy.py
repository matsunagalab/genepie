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
import tempfile
from pathlib import Path
import numpy as np
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..s_molecule import SMolecule
from ..s_trajectories import STrajectories, TRJ_TYPE_COOR_BOX
from .. import genesis_exe
from ._dcd_test_utils import rewrite_dcd


def test_rmsd_lazy_no_fitting():
    """Test lazy RMSD analysis without fitting using unified API."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)

    # Create lazy trajectory via crd_convert(lazy=True)
    lazy_trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=True,
    )

    lazy_traj = lazy_trajs[0]
    assert lazy_traj.is_lazy, "Trajectory should be lazy"

    # Run RMSD analysis with lazy trajectory
    result = genesis_exe.rmsd_analysis(
        subset_mol,
        lazy_traj,
        analysis_selection="an:CA",
        fitting_selection=None,  # No fitting
        ana_period=1,
        mass_weighted=False,
    )

    # Validate RMSD results
    assert result.rmsd is not None, "RMSD result should not be None"
    assert len(result.rmsd) > 0, "RMSD result should have at least one frame"
    assert all(r >= 0 for r in result.rmsd), "RMSD values should be non-negative"
    assert all(r < 50.0 for r in result.rmsd), "RMSD values should be reasonable (< 50 A)"

    print(f"Lazy RMSD no fitting (n={len(result.rmsd)}): "
          f"min={min(result.rmsd):.5f}, max={max(result.rmsd):.5f}")


def test_rmsd_lazy_with_fitting():
    """Test lazy RMSD analysis with TR+ROT fitting using unified API."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)

    # Create lazy trajectory via crd_convert(lazy=True)
    lazy_trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=True,
    )

    lazy_traj = lazy_trajs[0]

    # Run RMSD analysis with fitting
    result = genesis_exe.rmsd_analysis(
        subset_mol,
        lazy_traj,
        analysis_selection="sid:BPTI and an:CA",
        fitting_selection="sid:BPTI and an:CA",
        fitting_method="TR+ROT",
        ana_period=1,
        mass_weighted=False,
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
    trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=False,  # Memory mode
    )
    result = genesis_exe.rmsd_analysis(
        subset_mol, trajs[0],
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
    lazy_trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=True,  # Lazy mode
    )
    result = genesis_exe.rmsd_analysis(
        subset_mol, lazy_trajs[0],
        analysis_selection="sid:BPTI and an:CA",
        fitting_selection="sid:BPTI and an:CA",
        fitting_method="TR+ROT",
        ana_period=1,
        mass_weighted=False,
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

    # Create lazy trajectories
    lazy_trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=True,
    )
    lazy_traj = lazy_trajs[0]

    # Run with ana_period=1
    result_p1 = genesis_exe.rmsd_analysis(
        subset_mol,
        lazy_traj,
        analysis_selection="an:CA",
        fitting_selection=None,
        ana_period=1,
        mass_weighted=False,
    )

    # For ana_period=2, we need a fresh lazy trajectory (DCD file is read sequentially)
    lazy_trajs2, subset_mol2 = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=True,
    )

    result_p2 = genesis_exe.rmsd_analysis(
        subset_mol2,
        lazy_trajs2[0],
        analysis_selection="an:CA",
        fitting_selection=None,
        ana_period=2,
        mass_weighted=False,
    )

    # With ana_period=2, we should have about half the frames
    expected_frames = lazy_traj.nframe // 2
    assert len(result_p2.rmsd) == expected_frames, \
        f"Expected {expected_frames} frames with ana_period=2, got {len(result_p2.rmsd)}"

    print(f"ana_period test passed:")
    print(f"  ana_period=1: {len(result_p1.rmsd)} frames")
    print(f"  ana_period=2: {len(result_p2.rmsd)} frames")


def test_rmsd_lazy_big_endian():
    """Lazy coordinates must honor the DCD endian for every frame."""
    mol = SMolecule.from_file(
        pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        big_endian_dcd = Path(tmpdir) / "bpti_big_endian.dcd"
        rewrite_dcd(BPTI_DCD, big_endian_dcd, output_endian=">")

        little_trajs, little_mol = genesis_exe.crd_convert(
            mol, [str(BPTI_DCD)], trj_type="COOR+BOX", lazy=True
        )
        big_trajs, big_mol = genesis_exe.crd_convert(
            mol, [str(big_endian_dcd)], trj_type="COOR+BOX", lazy=True
        )
        little = genesis_exe.rmsd_analysis(
            little_mol, little_trajs[0], analysis_selection="an:CA"
        ).rmsd
        big = genesis_exe.rmsd_analysis(
            big_mol, big_trajs[0], analysis_selection="an:CA"
        ).rmsd
        np.testing.assert_allclose(big, little, rtol=1e-5, atol=1e-6)


def test_rmsd_lazy_compatibility_wrapper():
    """The legacy one-shot API delegates to the canonical lazy path."""
    mol = SMolecule.from_file(
        pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB
    )
    lazy_trajs, subset_mol = genesis_exe.crd_convert(
        mol, [str(BPTI_DCD)], trj_type="COOR+BOX", lazy=True
    )
    canonical = genesis_exe.rmsd_analysis(
        subset_mol,
        lazy_trajs[0],
        analysis_selection="an:CA",
        fitting_selection="an:CA",
    ).rmsd
    compatibility = genesis_exe.rmsd_analysis_lazy(
        mol,
        str(BPTI_DCD),
        analysis_selection="an:CA",
        fitting_selection="an:CA",
        has_box=True,
    )
    assert compatibility.dcd_nframe == lazy_trajs[0].nframe
    assert compatibility.dcd_natom == mol.num_atoms
    np.testing.assert_allclose(
        compatibility.rmsd, canonical, rtol=1e-5, atol=1e-6
    )


def test_rmsd_lazy_truncated_dcd():
    """A truncated frame must raise a catchable file-format exception."""
    from ..exceptions import GenesisFortranFileError

    mol = SMolecule.from_file(
        pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB
    )
    source = Path(BPTI_DCD).read_bytes()

    def assert_rejected(path):
        lazy_traj = STrajectories.from_lazy(
            str(path),
            TRJ_TYPE_COOR_BOX,
            nframe=10,
            natom=mol.num_atoms,
        )
        try:
            genesis_exe.rmsd_analysis(
                mol, lazy_traj, analysis_selection="an:CA"
            )
        except GenesisFortranFileError as exc:
            assert exc.code == 202
        else:
            raise AssertionError("Truncated DCD should raise a file error")

    with tempfile.TemporaryDirectory() as tmpdir:
        truncated = Path(tmpdir) / "truncated.dcd"
        truncated.write_bytes(source[:-17])
        assert_rejected(truncated)

        whole_frame = 3 * (8 + 4 * mol.num_atoms) + 56
        frame_truncated = Path(tmpdir) / "frame_truncated.dcd"
        frame_truncated.write_bytes(source[:-whole_frame])
        assert_rejected(frame_truncated)


def test_rmsd_lazy_missing_dcd_fortran_catch():
    """A missing DCD must surface as a catchable error, never crash the process.

    The Python wrapper normally guards with os.path.exists, but the real
    protection also lives in Fortran. Here we deliberately bypass the Python
    existence check on the canonical lazy-trajectory path.
    """
    from ..analysis import _common
    from ..exceptions import GenesisError, GenesisFortranFileError

    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
    missing = "/no/such/dir/missing.dcd"
    lazy_traj = STrajectories.from_lazy(
        missing,
        TRJ_TYPE_COOR_BOX,
        nframe=1,
        natom=mol.num_atoms,
    )

    real_exists = _common.os.path.exists
    _common.os.path.exists = \
        lambda p: True if str(p) == missing else real_exists(p)
    try:
        raised = None
        try:
            genesis_exe.rmsd_analysis(
                mol, lazy_traj, analysis_selection="an:CA"
            )
        except GenesisError as e:
            raised = e
    finally:
        _common.os.path.exists = real_exists

    assert raised is not None, "missing DCD must raise (not silently pass)"
    assert isinstance(raised, GenesisFortranFileError), \
        f"expected GenesisFortranFileError, got {type(raised).__name__}: {raised}"
    assert raised.code == 201, f"expected ERROR_FILE_NOT_FOUND (201), got {raised.code}"
    print(f"Fortran-level catch OK: {type(raised).__name__} code={raised.code}: {raised}")


def test_rmsd_lazy_dcd_atom_count_error():
    """A topology/DCD mismatch keeps the dedicated atom-count error code."""
    from ..exceptions import GenesisFortranValidationError

    mol = SMolecule.from_file(
        pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB
    )
    selection_indices = np.arange(
        1, mol.num_atoms + 1, dtype=np.int32
    )
    lazy_traj = STrajectories.from_lazy(
        str(BPTI_DCD),
        TRJ_TYPE_COOR_BOX,
        nframe=10,
        natom=mol.num_atoms + 1,
        selection_indices=selection_indices,
    )
    try:
        genesis_exe.rmsd_analysis(
            mol, lazy_traj, analysis_selection="an:CA"
        )
    except GenesisFortranValidationError as exc:
        assert exc.code == 304
    else:
        raise AssertionError("DCD atom-count mismatch should raise code 304")


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
        "test_rmsd_lazy_big_endian",
        "test_rmsd_lazy_compatibility_wrapper",
        "test_rmsd_lazy_truncated_dcd",
        "test_rmsd_lazy_missing_dcd_fortran_catch",
        "test_rmsd_lazy_dcd_atom_count_error",
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
