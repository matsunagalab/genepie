"""Average coordinates over a trajectory."""

import ctypes
import io
from collections import namedtuple
from typing import (
    Iterable,
    Optional,
)

import numpy as np

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
from ._common import _SourceArgs


AvecrdAnalysisResult = namedtuple(
        'AvecrdAnalysisResult',
        ['pdb', 'coord'])
AvecrdAnalysisResult.__doc__ = """Result of :func:`avecrd_analysis`.

Attributes:
    pdb: Averaged structure of the ``analysis_atom`` group as a PDB string.
    coord: Averaged coordinates of every atom of ``molecule``, shape
        ``(natom, 3)`` in angstrom (float64). Unlike ``pdb`` it is not rounded
        to three decimals, and it covers all atoms, so it can be fed directly
        to :func:`~genepie.analysis.pca.flccrd_analysis` and
        :func:`~genepie.analysis.pca.prjcrd_analysis`.
"""


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
        mass_weight: Optional[bool] = None,
        ):
    """
    Executes aa_analysis (iterative average structure).

    Every iteration fits each frame onto the current reference (initially the
    coordinates of ``molecule``), averages the fitted frames and makes the
    average the next reference, as ``avecrd_analysis`` does on the command
    line. Lazy trajectories (``crd_convert(..., lazy=True)``) are read from
    the DCD file once per iteration.

    Args:
        molecule: Molecular structure (initial reference coordinates)
        trajs: Trajectory to average (in memory or lazy)
        ana_period: Analysis period (default: 1)
        selection_group: ``[SELECTION]`` groups, e.g. ``["an:CA"]``
        selection_mole_name: ``[SELECTION]`` mole_name entries
        fitting_method: ``NO``, ``TR+ROT`` (typical), ``TR``, ``TR+ZROT``,
            ``XYTR`` or ``XYTR+ZROT``
        fitting_atom: 1-based index of the selection group used for fitting
        mass_weight: Mass-weighted fitting (default: ``[FITTING]`` default)
        check_only: ``[OPTION] check_only``
        num_iterations: Number of fit-and-average iterations
        analysis_atom: 1-based index of the selection group exported to
            ``pdb``
    Returns:
        AvecrdAnalysisResult(pdb, coord)
    """
    if ana_period is None:
        ana_period = 1
    pdb_ave_c = ctypes.c_void_p(None)

    src = _SourceArgs(trajs, ana_period)
    # Averaged coordinates of every atom: Python (natom, 3) C-order aliases
    # the Fortran (3, natom) buffer the wrapper fills.
    coord = np.zeros((molecule.num_atoms, 3), dtype=np.float64)

    try:
        with molecule_c(molecule) as mol_c:
            ctrl = io.BytesIO()
            ctrl_files.write_ctrl_output(
                    ctrl,
                    pdb_avefile="dummy.pdb")
            ctrl_files.write_ctrl_selection(
                    ctrl, selection_group, selection_mole_name)
            ctrl_files.write_ctrl_fitting(
                    ctrl, fitting_method, fitting_atom,
                    mass_weight=mass_weight)
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
                        *src.args,
                        ctrl_bytes,
                        ctrl_len,
                        coord.ctypes.data_as(ctypes.c_void_p),
                        ctypes.byref(pdb_ave_c),
                        ctypes.byref(status),
                        msgbuf,
                        ctypes.c_int(msglen),
                        )

        pdb_ave = c2py_util.conv_string(pdb_ave_c) if pdb_ave_c else None
        return AvecrdAnalysisResult(pdb_ave, coord)
    finally:
        if pdb_ave_c:
            LibGenesis().lib.deallocate_c_string(ctypes.byref(pdb_ave_c))
