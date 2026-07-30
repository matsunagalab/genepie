"""Radius of gyration."""

import ctypes
from collections import namedtuple
import numpy as np

from ..libgenesis import LibGenesis
from ..s_molecule import SMolecule
from ..s_trajectories import STrajectories
from ..exceptions import GenesisValidationError
from .._fortran import fortran_status
from .converter import selection
from ._common import _prepare_lazy_trajectory


RgAnalysisResult = namedtuple(
        'RgAnalysisResult',
        ['rg'])


def _rg_analysis_lazy(
        molecule: SMolecule,
        trajs: STrajectories,
        analysis_selection: str,
        ana_period: int = 1,
        mass_weighted: bool = True,
        ) -> "RgAnalysisResult":
    """
    Private implementation: RG analysis with lazy DCD loading.

    Called by rg_analysis() when trajs.is_lazy is True.
    """
    lib = LibGenesis().lib

    dcd_filename_bytes, source_selection, effective_period, n_result = \
        _prepare_lazy_trajectory(trajs, ana_period)
    trj_type = trajs.lazy_trj_type
    if molecule.num_atoms != trajs.natom:
        raise GenesisValidationError(
            "Lazy trajectory and molecule must contain the same selected atoms"
        )

    # Get atom indices using GENESIS selection
    analysis_indices = selection(molecule, analysis_selection)
    n_analysis = len(analysis_indices)

    # Ensure mass array is contiguous and correct dtype
    mass = np.ascontiguousarray(molecule.mass, dtype=np.float64)

    # Get pointer to mass array (zero-copy)
    mass_ptr = mass.ctypes.data_as(ctypes.c_void_p)

    # Pre-allocate exactly the number of analyzed frames.
    result_rg = np.zeros(n_result, dtype=np.float64)
    result_ptr = result_rg.ctypes.data_as(ctypes.c_void_p)

    filename_len = len(dcd_filename_bytes)

    # Output variables
    nstru_out = ctypes.c_int()
    dcd_nframe_out = ctypes.c_int()
    dcd_natom_out = ctypes.c_int()

    with fortran_status() as (status, msg, msglen):
        lib.rg_analysis_lazy_c(
            dcd_filename_bytes,
            ctypes.c_int(filename_len),
            ctypes.c_int(trj_type),
            ctypes.c_int(trajs.lazy_dcd_natom),
            source_selection.ctypes.data_as(ctypes.c_void_p),
            ctypes.c_int(trajs.natom),
            mass_ptr,
            ctypes.c_int(molecule.num_atoms),
            ctypes.c_int(effective_period),
            analysis_indices.ctypes.data_as(ctypes.c_void_p),
            ctypes.c_int(n_analysis),
            ctypes.c_int(1 if mass_weighted else 0),
            result_ptr,
            ctypes.c_int(n_result),
            ctypes.byref(nstru_out),
            ctypes.byref(dcd_nframe_out),
            ctypes.byref(dcd_natom_out),
            ctypes.byref(status),
            msg,
            ctypes.c_int(msglen),
        )


    # Return RgAnalysisResult (same as non-lazy version)
    return RgAnalysisResult(result_rg[:nstru_out.value])


def rg_analysis(
        molecule: SMolecule,
        trajs: STrajectories,
        analysis_selection: str,
        ana_period: int = 1,
        mass_weighted: bool = True,
        ) -> RgAnalysisResult:
    """
    Executes radius of gyration (RG) analysis.

    This function calculates the radius of gyration for selected atoms
    across all trajectory frames.

    Args:
        molecule: Molecular structure
        trajs: Trajectories to analyze (memory or lazy)
        analysis_selection: GENESIS selection string (e.g., "an:CA", "heavy")
        ana_period: Analysis period (default: 1)
        mass_weighted: Use mass weighting for RG calculation (default: True)

    Returns:
        RgAnalysisResult containing the radius of gyration array

    Examples:
        >>> # Memory-based (standard)
        >>> result = rg_analysis(mol, trajs, "an:CA")
        >>> print(result.rg)

        >>> # Lazy mode: works with lazy trajectories from crd_convert(lazy=True)
        >>> lazy_trajs, mol = crd_convert(mol, ["traj.dcd"], lazy=True)
        >>> result = rg_analysis(mol, lazy_trajs[0], "an:CA")
    """
    # Handle lazy trajectories
    if trajs.is_lazy:
        return _rg_analysis_lazy(
            molecule=molecule,
            trajs=trajs,
            analysis_selection=analysis_selection,
            ana_period=ana_period,
            mass_weighted=mass_weighted,
        )

    lib = LibGenesis().lib

    # Get atom indices using GENESIS selection
    analysis_indices = selection(molecule, analysis_selection)
    n_analysis = len(analysis_indices)

    # Ensure mass array is contiguous and correct dtype
    mass = np.ascontiguousarray(molecule.mass, dtype=np.float64)

    # Get pointer to mass array (zero-copy)
    mass_ptr = mass.ctypes.data_as(ctypes.c_void_p)

    # Pre-allocate result array (zero-copy)
    n_frame = int(trajs.nframe / ana_period)
    result_rg = np.zeros(n_frame, dtype=np.float64)
    result_ptr = result_rg.ctypes.data_as(ctypes.c_void_p)

    # Output variables
    nstru_out = ctypes.c_int()

    with fortran_status() as (status, msg, msglen):
        lib.rg_analysis_c(
            mass_ptr,
            ctypes.c_int(molecule.num_atoms),
            ctypes.byref(trajs.get_c_obj()),
            ctypes.c_int(ana_period),
            analysis_indices.ctypes.data_as(ctypes.c_void_p),
            ctypes.c_int(n_analysis),
            ctypes.c_int(1 if mass_weighted else 0),
            result_ptr,
            ctypes.c_int(n_frame),
            ctypes.byref(nstru_out),
            ctypes.byref(status),
            msg,
            ctypes.c_int(msglen),
        )


    # Return the pre-allocated result array (already filled by Fortran)
    return RgAnalysisResult(result_rg[:nstru_out.value])
