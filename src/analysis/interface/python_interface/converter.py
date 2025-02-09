import ctypes
import os
from libgenesis import LibGenesis
from s_molecule import SMolecule, py2c_s_molecule
from s_trajectories import STrajectoriesArray
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
