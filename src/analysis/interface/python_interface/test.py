import ctypes
import os
import pathlib
import numpy as np
import numpy.typing as npt
from libgenesis import LibGenesis
from s_molecule import SMolecule, py2c_s_molecule
from s_trajectories import STrajectories
import genesis_exe
import py2c_util
import c2py_util


def test():
    # 関数を呼び出す
    pdb_filename = pathlib.Path("molecule.pdb")
    with SMolecule.from_pdb_file(pdb_filename) as mol:
        # 結果を処理する
        print("num_atoms = ", mol.num_atoms)
        for i in range(max(0, mol.num_atoms - 5), mol.num_atoms):
            print(mol.atom_coord[i])
            print(mol.atom_no[i], mol.segment_name[i], mol.atom_name[i])

        print("num_atoms = ", mol.num_atoms)
        mol_c = py2c_s_molecule(mol)
        LibGenesis().lib.test_conv_c2f(ctypes.byref(mol_c))


def test_crd():
    # 関数を呼び出す
    pdb_path = pathlib.Path(
            "../../../../tests/regression_test/test_analysis/trajectories/BPTI_charmm/BPTI_ionize.pdb")
    psf_path = pathlib.Path(
            "../../../../tests/regression_test/test_analysis/trajectories/BPTI_charmm/BPTI_ionize.psf")
    ctrl_path = pathlib.Path("./test_crd_inp")
            # "./test_crd_inp../../../../tests/regression_test/test_analysis/test_crd_convert/BPTI/inp")
    # traj_path = pathlib.Path(
    #         "../../../../tests/regression_test/test_analysis/trajectories/BPTI_charmm/BPTI_run.dcd")
    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, ctrl_path) as trajs:
            pass


def test_trj_analysis():
    # 関数を呼び出す
    pdb_path = pathlib.Path(
            "../../../../tests/regression_test/test_analysis/trajectories/BPTI_charmm/BPTI_ionize.pdb")
    psf_path = pathlib.Path(
            "../../../../tests/regression_test/test_analysis/trajectories/BPTI_charmm/BPTI_ionize.psf")
    crd_ctrl_path = pathlib.Path("./test_crd_inp")
    trj_analysis_ctrl_path = pathlib.Path("./test_trj_analysis_inp")

    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                d = genesis_exe.trj_analysis(
                        mol, t, 1, trj_analysis_ctrl_path)
                print(d)


def test_rg_analysis():
    # 関数を呼び出す
    pdb_path = pathlib.Path(
            "../../../../tests/regression_test/test_analysis/trajectories/BPTI_charmm/BPTI_ionize.pdb")
    psf_path = pathlib.Path(
            "../../../../tests/regression_test/test_analysis/trajectories/BPTI_charmm/BPTI_ionize.psf")
    crd_ctrl_path = pathlib.Path("./test_crd_inp")
    rg_analysis_ctrl_path = pathlib.Path("./test_rg_analysis_inp")

    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                rg = genesis_exe.rg_analysis(mol, t, 1, rg_analysis_ctrl_path)
                print(rg)


if __name__ == "__main__":
    # test()
    # test_crd()
    # test_trj_analysis()
    test_rg_analysis()
