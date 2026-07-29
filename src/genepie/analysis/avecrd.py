"""Average coordinates over a trajectory."""

import ctypes
import io
from collections import namedtuple
from typing import (
    Iterable,
    Optional,
)

from ..libgenesis import LibGenesis
from ..s_molecule import SMolecule
from ..s_trajectories import STrajectories
from .. import ctrl_files
from .. import c2py_util
from .._fortran import (
    ctrl_to_bytes,
    fortran_status,
    molecule_c,
)


AvecrdAnalysisResult = namedtuple(
        'AvecrdAnalysisResult',
        ['pdb'])


def avecrd_analysis(
        molecule: SMolecule, trajs: STrajectories,
        ana_period: Optional[int] = 1,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        fitting_method: Optional[str] = None,
        fitting_atom: Optional[int] = None,
        check_only: Optional[bool] = None,
        num_iterations: Optional[int] = None,
        analysis_atom: Optional[int] = None,
        ):
    """
    Executes aa_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
    Returns:
        (pdb,)
    """
    if ana_period is None:
        ana_period = 1
    pdb_ave_c = ctypes.c_void_p(None)
    ana_period_c = ctypes.c_int(ana_period)

    try:
        with molecule_c(molecule) as mol_c:
            ctrl = io.BytesIO()
            ctrl_files.write_ctrl_output(
                    ctrl,
                    pdb_avefile="dummy.pdb")
            ctrl_files.write_ctrl_selection(
                    ctrl, selection_group, selection_mole_name)
            ctrl_files.write_ctrl_fitting(
                    ctrl, fitting_method, fitting_atom)
            ctrl.write(b"[OPTION]\n")
            ctrl_files.write_kwargs(
                    ctrl,
                    check_only=check_only,
                    num_iterations=num_iterations,
                    analysis_atom=analysis_atom,
                    )

            ctrl_bytes, ctrl_len = ctrl_to_bytes(ctrl)
            with fortran_status() as (status, msgbuf, msglen):
                LibGenesis().lib.aa_analysis_c(
                        ctypes.byref(mol_c),
                        ctypes.byref(trajs.get_c_obj()),
                        ctypes.byref(ana_period_c),
                        ctrl_bytes,
                        ctrl_len,
                        ctypes.byref(pdb_ave_c),
                        ctypes.byref(status),
                        msgbuf,
                        ctypes.c_int(msglen),
                        )

        pdb_ave = c2py_util.conv_string(pdb_ave_c) if pdb_ave_c else None
        return AvecrdAnalysisResult(pdb_ave)
    finally:
        if pdb_ave_c:
            LibGenesis().lib.deallocate_c_string(ctypes.byref(pdb_ave_c))
