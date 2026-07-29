"""Hydrogen bond analysis."""

import ctypes
import io
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


def hb_analysis(molecule: SMolecule, trajs: STrajectories,
                ana_period: Optional[int] = 1,
                selection_group: Optional[Iterable[str]] = None,
                selection_mole_name: Optional[Iterable[str]] = None,
                check_only: Optional[bool] = None,
                output_type: Optional[str] = None,
                solvent_list: Optional[str] = None,
                analysis_atom: Optional[int] = None,
                target_atom: Optional[int] = None,
                boundary_type: Optional[str] = None,
                hb_distance: Optional[float] = None,
                dha_angle: Optional[float] = None,
                hda_angle: Optional[float] = None,
                ) -> str:
    """
    Executes hb_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
    Returns:
        result
    """
    if ana_period is None:
        ana_period = 1
    result = ctypes.c_void_p(None)
    ana_period_c = ctypes.c_int(ana_period)

    try:
        with molecule_c(molecule) as mol_c:
            ctrl = io.BytesIO()
            ctrl_files.write_ctrl_output(
                    ctrl,
                    outfile="dummy.out")
            ctrl_files.write_ctrl_selection(
                    ctrl, selection_group, selection_mole_name)
            ctrl.write(b"[OPTION]\n")
            ctrl_files.write_kwargs(
                    ctrl,
                    check_only=check_only,
                    output_type=output_type,
                    solvent_list=solvent_list,
                    analysis_atom=analysis_atom,
                    target_atom=target_atom,
                    boundary_type=boundary_type,
                    hb_distance=hb_distance,
                    dha_angle=dha_angle,
                    hda_angle=hda_angle,
                    )

            ctrl_bytes, ctrl_len = ctrl_to_bytes(ctrl)
            with fortran_status() as (status, msgbuf, msglen):
                LibGenesis().lib.hb_analysis_c(
                        ctypes.byref(mol_c),
                        ctypes.byref(trajs.get_c_obj()),
                        ctypes.byref(ana_period_c),
                        ctrl_bytes,
                        ctrl_len,
                        ctypes.byref(result),
                        ctypes.byref(status),
                        msgbuf,
                        ctypes.c_int(msglen),
                        )

        return c2py_util.conv_string(result)

    finally:
        if result:
            try:
                LibGenesis().lib.deallocate_c_string(ctypes.byref(result))
            except Exception:
                pass
