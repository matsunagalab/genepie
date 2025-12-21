#!/usr/bin/env python3
"""
Test script for atdyn Python interface.
Covers all regression test categories from test.py "mpirun -np 1 atdyn".

Test categories:
1. AMBER format (glycam) - MIN and MD
2. GROMACS format (bpti) - MD
3. CHARMM format (jac_param27, dppc) - MIN and MD

Note: Each test is run in a separate subprocess to isolate Fortran library
global state (timers, PME FFT plans, memory allocations) between test runs.
This is necessary because the Fortran GENESIS library maintains internal state
that cannot be fully reset between consecutive runs in the same process.
The subprocess isolation ensures reliable, reproducible test results.
"""

import sys
import os
import subprocess
import tempfile

# Add path for development
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def get_base_dir():
    """Get the base directory of the GENESIS repository."""
    return os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))


def get_test_dir(system_name):
    """Get the test directory for a given system."""
    return os.path.join(get_base_dir(), 'tests', 'regression_test', 'build', system_name)


def get_param_dir():
    """Get the parameter file directory."""
    return os.path.join(get_base_dir(), 'tests', 'regression_test', 'param')


def run_test_in_subprocess(test_func_name):
    """Run a test function in a separate subprocess to isolate Fortran state."""
    script = f'''
import sys
import os
sys.path.insert(0, "{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}")
from python_interface import test_atdyn
result = test_atdyn.{test_func_name}()
sys.exit(0 if result else 1)
'''
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '1'

    result = subprocess.run(
        [sys.executable, '-c', script],
        capture_output=True,
        text=True,
        env=env,
        cwd=os.path.dirname(os.path.abspath(__file__))
    )

    # Print output from subprocess
    if result.stdout:
        print(result.stdout, end='')
    if result.stderr:
        print(result.stderr, end='', file=sys.stderr)

    return result.returncode == 0


# =============================================================================
# AMBER format tests (glycam)
# =============================================================================

def test_atdyn_min_glycam():
    """Test energy minimization with AMBER/glycam system."""
    from python_interface import genesis_exe

    print("=" * 60)
    print("Testing run_atdyn_min (AMBER format: glycam/MINIMIZE)")
    print("=" * 60)

    test_dir = get_test_dir('glycam')
    prmtopfile = os.path.join(test_dir, 'glycam.top')
    ambcrdfile = os.path.join(test_dir, 'glycam.rst')
    rstfile = os.path.join(test_dir, 'rst')

    print(f"  prmtopfile: {prmtopfile}")

    if not os.path.exists(prmtopfile):
        print(f"Test files not found: {prmtopfile}")
        print("SKIPPED")
        return True

    try:
        result = genesis_exe.run_atdyn_min(
            prmtopfile=prmtopfile,
            ambcrdfile=ambcrdfile,
            rstfile=rstfile,
            forcefield="AMBER",
            electrostatic="PME",
            switchdist=12.0,
            cutoffdist=12.0,
            pairlistdist=14.0,
            pme_alpha=0.34,
            pme_ngrid_x=64,
            pme_ngrid_y=64,
            pme_ngrid_z=64,
            pme_nspline=4,
            dispersion_corr="epress",
            method="SD",
            nsteps=20,
            eneout_period=2,
            nbupdate_period=4,
            rigid_bond=False,
            boundary_type="PBC",
            box_size_x=69.5294360,
            box_size_y=68.0597930,
            box_size_z=56.2256950,
        )

        print(f"  Number of atoms: {result.final_coords.shape[1]}")
        print(f"  Final energy (total): {result.energies[0, 0]:.4f}")

        if result.energies[0, 0] > 0:
            print("ERROR: Total energy should be negative")
            return False

        print("PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_atdyn_md_glycam():
    """Test MD with AMBER/glycam system (glycam/PME)."""
    from python_interface import genesis_exe

    print("=" * 60)
    print("Testing run_atdyn_md (AMBER format: glycam/PME)")
    print("=" * 60)

    test_dir = get_test_dir('glycam')
    prmtopfile = os.path.join(test_dir, 'glycam.top')
    ambcrdfile = os.path.join(test_dir, 'glycam.rst')
    rstfile = os.path.join(test_dir, 'rst')

    print(f"  prmtopfile: {prmtopfile}")

    if not os.path.exists(prmtopfile):
        print(f"Test files not found: {prmtopfile}")
        print("SKIPPED")
        return True

    try:
        result = genesis_exe.run_atdyn_md(
            prmtopfile=prmtopfile,
            ambcrdfile=ambcrdfile,
            rstfile=rstfile,
            forcefield="AMBER",
            electrostatic="PME",
            switchdist=12.0,
            cutoffdist=12.0,
            pairlistdist=14.0,
            pme_alpha=0.34,
            pme_ngrid_x=64,
            pme_ngrid_y=64,
            pme_ngrid_z=64,
            pme_nspline=4,
            dispersion_corr="epress",
            integrator="VVER",
            nsteps=20,
            timestep=0.001,
            eneout_period=2,
            nbupdate_period=5,
            iseed=314159,
            verbose=True,
            rigid_bond=True,
            shake_iteration=500,
            shake_tolerance=1.0e-10,
            water_model="WAT",
            ensemble="NVE",
            tpcontrol="NO",
            temperature=0,
            boundary_type="PBC",
            box_size_x=69.5294360,
            box_size_y=68.0597930,
            box_size_z=56.2256950,
        )

        print(f"  Number of atoms: {result.final_coords.shape[1]}")
        print(f"  Energy shape: {result.energies.shape}")

        if result.energies[0, 0] > 0:
            print("ERROR: Total energy should be negative")
            return False

        print("PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


# =============================================================================
# GROMACS format tests (bpti)
# =============================================================================

def test_atdyn_md_bpti():
    """Test MD with GROMACS/bpti system (bpti/PME_WATER)."""
    from python_interface import genesis_exe

    print("=" * 60)
    print("Testing run_atdyn_md (GROMACS format: bpti/PME_WATER)")
    print("=" * 60)

    test_dir = get_test_dir('bpti')
    grotopfile = os.path.join(test_dir, 'bpti.top')
    grocrdfile = os.path.join(test_dir, 'bpti.gro')
    rstfile = os.path.join(test_dir, 'rst')

    print(f"  grotopfile: {grotopfile}")

    if not os.path.exists(grotopfile):
        print(f"Test files not found: {grotopfile}")
        print("SKIPPED")
        return True

    try:
        result = genesis_exe.run_atdyn_md(
            grotopfile=grotopfile,
            grocrdfile=grocrdfile,
            rstfile=rstfile,
            forcefield="GROAMBER",
            electrostatic="PME",
            switchdist=12.0,
            cutoffdist=12.0,
            pairlistdist=14.0,
            pme_alpha=0.34,
            pme_ngrid_x=64,
            pme_ngrid_y=64,
            pme_ngrid_z=64,
            pme_nspline=4,
            output_style="GENESIS",
            integrator="VVER",
            nsteps=20,
            timestep=0.001,
            eneout_period=2,
            nbupdate_period=5,
            iseed=314159,
            verbose=True,
            rigid_bond=True,
            shake_iteration=500,
            shake_tolerance=1.0e-10,
            water_model="SOL",
            ensemble="NVE",
            tpcontrol="NO",
            temperature=0,
            boundary_type="PBC",
            box_size_x=65.3318,
            box_size_y=65.3318,
            box_size_z=65.3318,
        )

        print(f"  Number of atoms: {result.final_coords.shape[1]}")
        print(f"  Energy shape: {result.energies.shape}")

        if result.energies[0, 0] > 0:
            print("ERROR: Total energy should be negative")
            return False

        print("PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


# =============================================================================
# CHARMM format tests (jac_param27, dppc)
# =============================================================================

def test_atdyn_md_jac_param27():
    """Test MD with CHARMM/jac_param27 system (jac_param27/PME_WATER)."""
    from python_interface import genesis_exe

    print("=" * 60)
    print("Testing run_atdyn_md (CHARMM format: jac_param27/PME_WATER)")
    print("=" * 60)

    test_dir = get_test_dir('jac_param27')
    param_dir = get_param_dir()
    topfile = os.path.join(param_dir, 'top_all27_prot_lipid.rtf')
    parfile = os.path.join(param_dir, 'par_all27_prot_lipid.prm')
    psffile = os.path.join(test_dir, 'jac_param27.psf')
    pdbfile = os.path.join(test_dir, 'jac_param27.pdb')
    rstfile = os.path.join(test_dir, 'rst')

    print(f"  psffile: {psffile}")

    if not os.path.exists(psffile):
        print(f"Test files not found: {psffile}")
        print("SKIPPED")
        return True

    try:
        result = genesis_exe.run_atdyn_md(
            topfile=topfile,
            parfile=parfile,
            psffile=psffile,
            pdbfile=pdbfile,
            rstfile=rstfile,
            forcefield="CHARMM",
            electrostatic="PME",
            switchdist=8.0,
            cutoffdist=10.0,
            pairlistdist=12.0,
            pme_alpha=0.34,
            pme_ngrid_x=64,
            pme_ngrid_y=64,
            pme_ngrid_z=64,
            pme_nspline=4,
            vdw_force_switch=False,
            output_style="GENESIS",
            integrator="LEAP",
            nsteps=20,
            timestep=0.001,
            eneout_period=2,
            nbupdate_period=5,
            iseed=314159,
            verbose=True,
            rigid_bond=True,
            shake_iteration=500,
            shake_tolerance=1.0e-10,
            water_model="TIP3",
            ensemble="NVE",
            tpcontrol="NO",
            temperature=0,
            boundary_type="PBC",
            box_size_x=65.5,
            box_size_y=65.5,
            box_size_z=65.5,
        )

        print(f"  Number of atoms: {result.final_coords.shape[1]}")
        print(f"  Energy shape: {result.energies.shape}")

        if result.energies[0, 0] > 0:
            print("ERROR: Total energy should be negative")
            return False

        print("PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_atdyn_md_dppc_nvt():
    """Test MD with CHARMM/dppc system with NVT ensemble (dppc/VVER_NVT_LANGEVIN)."""
    from python_interface import genesis_exe

    print("=" * 60)
    print("Testing run_atdyn_md (CHARMM format: dppc/VVER_NVT_LANGEVIN)")
    print("=" * 60)

    test_dir = get_test_dir('dppc')
    param_dir = get_param_dir()
    topfile = os.path.join(param_dir, 'top_all36_lipid.rtf')
    parfile = os.path.join(param_dir, 'par_all36_lipid.prm')
    strfile = os.path.join(param_dir, 'toppar_water_ions.str')
    psffile = os.path.join(test_dir, 'dppc.psf')
    pdbfile = os.path.join(test_dir, 'dppc.pdb')
    rstfile = os.path.join(test_dir, 'rst')

    print(f"  psffile: {psffile}")

    if not os.path.exists(psffile):
        print(f"Test files not found: {psffile}")
        print("SKIPPED")
        return True

    try:
        result = genesis_exe.run_atdyn_md(
            topfile=topfile,
            parfile=parfile,
            strfile=strfile,
            psffile=psffile,
            pdbfile=pdbfile,
            rstfile=rstfile,
            forcefield="CHARMM",
            electrostatic="PME",
            switchdist=10.0,
            cutoffdist=12.0,
            pairlistdist=13.5,
            pme_alpha=0.34,
            pme_ngrid_x=72,
            pme_ngrid_y=72,
            pme_ngrid_z=72,
            pme_nspline=4,
            output_style="GENESIS",
            integrator="VVER",
            nsteps=20,
            timestep=0.001,
            eneout_period=2,
            nbupdate_period=5,
            iseed=314159,
            verbose=True,
            rigid_bond=True,
            shake_iteration=500,
            shake_tolerance=1.0e-10,
            ensemble="NVT",
            tpcontrol="LANGEVIN",
            temperature=300.0,
            pressure=1.0,
            boundary_type="PBC",
            box_size_x=69.4792,
            box_size_y=69.4792,
            box_size_z=71.6508,
        )

        print(f"  Number of atoms: {result.final_coords.shape[1]}")
        print(f"  Energy shape: {result.energies.shape}")

        if result.energies[0, 0] > 0:
            print("ERROR: Total energy should be negative")
            return False

        print("PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_atdyn_min_dppc():
    """Test minimization with CHARMM/dppc system (dppc/MINIMIZE)."""
    from python_interface import genesis_exe

    print("=" * 60)
    print("Testing run_atdyn_min (CHARMM format: dppc/MINIMIZE)")
    print("=" * 60)

    test_dir = get_test_dir('dppc')
    param_dir = get_param_dir()
    topfile = os.path.join(param_dir, 'top_all36_lipid.rtf')
    parfile = os.path.join(param_dir, 'par_all36_lipid.prm')
    strfile = os.path.join(param_dir, 'toppar_water_ions.str')
    psffile = os.path.join(test_dir, 'dppc.psf')
    pdbfile = os.path.join(test_dir, 'dppc.pdb')
    rstfile = os.path.join(test_dir, 'rst')

    print(f"  psffile: {psffile}")

    if not os.path.exists(psffile):
        print(f"Test files not found: {psffile}")
        print("SKIPPED")
        return True

    try:
        result = genesis_exe.run_atdyn_min(
            topfile=topfile,
            parfile=parfile,
            strfile=strfile,
            psffile=psffile,
            pdbfile=pdbfile,
            rstfile=rstfile,
            forcefield="CHARMM",
            electrostatic="PME",
            switchdist=10.0,
            cutoffdist=12.0,
            pairlistdist=13.5,
            pme_alpha=0.34,
            pme_ngrid_x=72,
            pme_ngrid_y=72,
            pme_ngrid_z=72,
            pme_nspline=4,
            output_style="GENESIS",
            method="SD",
            nsteps=20,
            eneout_period=2,
            nbupdate_period=4,
            rigid_bond=False,
            boundary_type="PBC",
            box_size_x=69.4792,
            box_size_y=69.4792,
            box_size_z=71.6508,
        )

        print(f"  Number of atoms: {result.final_coords.shape[1]}")
        print(f"  Final energy (total): {result.energies[0, 0]:.4f}")

        if result.energies[0, 0] > 0:
            print("ERROR: Total energy should be negative")
            return False

        print("PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


# =============================================================================
# Main
# =============================================================================

def run_all_tests_isolated():
    """Run all tests in separate subprocesses to isolate Fortran state."""
    tests = [
        ("AMBER MIN (glycam/MINIMIZE)", "test_atdyn_min_glycam"),
        ("AMBER MD (glycam/PME)", "test_atdyn_md_glycam"),
        ("GROMACS MD (bpti/PME_WATER)", "test_atdyn_md_bpti"),
        ("CHARMM MD (jac_param27/PME_WATER)", "test_atdyn_md_jac_param27"),
        ("CHARMM MD NVT (dppc/VVER_NVT_LANGEVIN)", "test_atdyn_md_dppc_nvt"),
        ("CHARMM MIN (dppc/MINIMIZE)", "test_atdyn_min_dppc"),
    ]

    results = []
    for name, func_name in tests:
        passed = run_test_in_subprocess(func_name)
        results.append((name, passed))

    print()
    print("=" * 60)
    print("Summary")
    print("=" * 60)

    all_passed = True
    passed_count = 0
    total_count = len(results)
    for name, passed in results:
        status = "PASSED" if passed else "FAILED"
        print(f"  {name}: {status}")
        if passed:
            passed_count += 1
        else:
            all_passed = False

    print()
    print(f"Total: {passed_count}/{total_count} passed")

    return all_passed


def run_single_test(test_name):
    """Run a single test directly (for subprocess invocation)."""
    test_funcs = {
        "test_atdyn_min_glycam": test_atdyn_min_glycam,
        "test_atdyn_md_glycam": test_atdyn_md_glycam,
        "test_atdyn_md_bpti": test_atdyn_md_bpti,
        "test_atdyn_md_jac_param27": test_atdyn_md_jac_param27,
        "test_atdyn_md_dppc_nvt": test_atdyn_md_dppc_nvt,
        "test_atdyn_min_dppc": test_atdyn_min_dppc,
    }

    if test_name in test_funcs:
        return test_funcs[test_name]()
    else:
        print(f"Unknown test: {test_name}")
        return False


if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Run single test (for subprocess or direct invocation)
        test_name = sys.argv[1]
        result = run_single_test(test_name)
        sys.exit(0 if result else 1)
    else:
        # Run all tests in isolated subprocesses
        all_passed = run_all_tests_isolated()
        sys.exit(0 if all_passed else 1)
