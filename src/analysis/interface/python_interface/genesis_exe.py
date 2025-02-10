import ctypes
import os
import numpy as np
import numpy.typing as npt
from libgenesis import LibGenesis
from s_molecule import SMolecule, py2c_s_molecule
from s_trajectories import STrajectories, STrajectoriesArray
import c2py_util
import py2c_util


def crd_convert(molecule: SMolecule,
                ctrl_filename: str | bytes | os.PathLike) -> STrajectoriesArray:
    """
    Executes crd_convert.

    Args:
        molecule:
        ctrl_filename:

    Returns:
        Array of STrajectories
    """
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


def trj_analysis(molecule: SMolecule, trajs :STrajectories,
                 ana_period: int,
                 ctrl_path: str | bytes | os.PathLike
                 ) -> tuple(npt.NDArray[np.float64]):
    """
    Executes trj_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
        ctrl_path:

    Returns:
        distance
    """
    mol_c = py2c_s_molecule(molecule)
    num_distance = ctypes.c_int(0)
    ana_period_c = ctypes.c_int(ana_period)
    result_distance_c = ctypes.c_void_p(None)
    LibGenesis().lib.trj_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(trajs.get_c_obj()),
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


def rg_analysis(molecule: SMolecule, trajs :STrajectories,
                ana_period: int,
                ctrl_path: str | bytes | os.PathLike
                ) -> npt.NDArray[np.float64]:
    mol_c = py2c_s_molecule(molecule)

    ana_period_c = ctypes.c_int(ana_period)
    result_rg_c = ctypes.c_void_p(None)

    LibGenesis().lib.rg_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(trajs.get_c_obj()),
            ctypes.byref(ana_period_c),
            py2c_util.pathlike_to_byte(ctrl_path),
            ctypes.byref(result_rg_c),
            )

    n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
    result_rg = c2py_util.conv_double_ndarray(
            result_rg_c, n_frame_c.value)
    LibGenesis().lib.deallocate_double(
            ctypes.byref(result_rg_c),
            ctypes.byref(n_frame_c))
    return result_rg


def rmsd_analysis(molecule: SMolecule, trajs :STrajectories,
                ana_period: int,
                ctrl_path: str | bytes | os.PathLike
                ) -> npt.NDArray[np.float64]:
    mol_c = py2c_s_molecule(molecule)

    ana_period_c = ctypes.c_int(ana_period)
    result_ra_c = ctypes.c_void_p(None)

    LibGenesis().lib.ra_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(trajs.get_c_obj()),
            ctypes.byref(ana_period_c),
            py2c_util.pathlike_to_byte(ctrl_path),
            ctypes.byref(result_ra_c),
            )

    n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
    result_ra = c2py_util.conv_double_ndarray(
            result_ra_c, n_frame_c.value)
    LibGenesis().lib.deallocate_double(
            ctypes.byref(result_ra_c),
            ctypes.byref(n_frame_c))
    return result_ra
