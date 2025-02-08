import ctypes
import os
import pathlib
import numpy as np
import numpy.typing as npt
from libgenesis import LibGenesis
from s_molecule import SMolecule, py2c_s_molecule
from s_trajectories import STrajectoriesArray
from s_trajectories_c import STrajectoriesC
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


def crd_convert(molecule: SMolecule,
                ctrl_filename: str | bytes | os.PathLike) -> STrajectoriesArray:
    buf = ctypes.c_void_p(None)
    num_trajs_c = ctypes.c_int(0)
    mol_c = py2c_s_molecule(molecule)
    LibGenesis().lib.crd_convert_c(
            ctypes.byref(mol_c),
            py2c_util.pathlike_to_byte(ctrl_filename),
            ctypes.byref(buf),
            ctypes.byref(num_trajs_c))
    LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))
    return STrajectoriesArray(buf, num_trajs_c.value)


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
        with crd_convert(mol, ctrl_path) as trajs:
            pass


def trj_analysis(molecule: SMolecule, trajs :STrajectoriesC,
                 ana_period: int,
                 ctrl_path: str | bytes | os.PathLike
                 ) -> tuple(npt.NDArray[np.float64]):
    mol_c = py2c_s_molecule(molecule)

    num_distance = ctypes.c_int(0)
    ana_period_c = ctypes.c_int(ana_period)
    result_distance_c = ctypes.c_void_p(None)

    LibGenesis().lib.trj_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(trajs),
            ctypes.byref(ana_period_c),
            py2c_util.pathlike_to_byte(ctrl_path),
            ctypes.byref(result_distance_c),
            ctypes.byref(num_distance),
            )

    n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
    result_distance = c2py_util.conv_double_ndarray(
            result_distance_c, [n_frame_c.value, num_distance.value])
    LibGenesis().lib.deallocate_double2(
            ctypes.byref(result_distance_c),
            ctypes.byref(n_frame_c), ctypes.byref(num_distance))
    return result_distance


def test_trj_analysis():
    # 関数を呼び出す
    pdb_path = pathlib.Path(
            "../../../../tests/regression_test/test_analysis/trajectories/BPTI_charmm/BPTI_ionize.pdb")
    psf_path = pathlib.Path(
            "../../../../tests/regression_test/test_analysis/trajectories/BPTI_charmm/BPTI_ionize.psf")
    crd_ctrl_path = pathlib.Path("./test_crd_inp")
    trj_analysis_ctrl_path = pathlib.Path("./test_trj_analysis_inp")

    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with crd_convert(mol, crd_ctrl_path) as trajs:
            d = trj_analysis(mol, trajs.traj_p[0], 1, trj_analysis_ctrl_path)
            print(d)


if __name__ == "__main__":
    # test()
    # test_crd()
    test_trj_analysis()
