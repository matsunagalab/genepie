import ctypes
import os
import numpy as np
import numpy.typing as npt
from libgenesis import LibGenesis
from s_molecule import SMolecule, py2c_s_molecule
from s_trajectories_c import STrajectoriesC
import c2py_util
import py2c_util


def trj_analysis(molecule: SMolecule, trajs :STrajectoriesC,
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
