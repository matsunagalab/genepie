"""Distance root mean square deviation."""

import ctypes
import os
from collections import namedtuple
import numpy as np

from ..libgenesis import LibGenesis
from ..s_trajectories import STrajectories
from ..exceptions import GenesisValidationError
from .._fortran import fortran_status


DrmsAnalysisResult = namedtuple(
        'DrmsAnalysisResult',
        ['drms'])


def _drms_analysis_lazy(
        trajs: STrajectories,
        contact_list: np.ndarray,
        contact_dist: np.ndarray,
        ana_period: int = 1,
        pbc_correct: bool = False,
        ) -> "DrmsAnalysisResult":
    """
    Private implementation: DRMS analysis with lazy DCD loading.

    Called by drms_analysis() when trajs.is_lazy is True.
    """
    lib = LibGenesis().lib

    # Extract lazy DCD info from STrajectories
    dcd_file = trajs.lazy_dcd_file
    trj_type = trajs.lazy_trj_type  # 1=COOR, 2=COOR+BOX
    n_atoms = trajs.natom

    # Validate DCD file exists
    if not os.path.exists(dcd_file):
        raise GenesisValidationError(f"DCD file not found: {dcd_file}")

    # Ensure arrays are contiguous and correct dtype
    contact_list_f = np.asfortranarray(contact_list, dtype=np.int32)
    contact_dist_f = np.ascontiguousarray(contact_dist, dtype=np.float64)

    # Validate shapes
    if contact_list_f.ndim != 2 or contact_list_f.shape[0] != 2:
        raise GenesisValidationError(
            f"contact_list must be shape (2, n_contact), got {contact_list.shape}"
        )
    n_contact = contact_list_f.shape[1]
    if contact_dist_f.shape[0] != n_contact:
        raise GenesisValidationError(
            f"contact_dist must have {n_contact} elements, got {contact_dist.shape[0]}"
        )

    # Get pointers (zero-copy input)
    contact_list_ptr = contact_list_f.ctypes.data_as(ctypes.c_void_p)
    contact_dist_ptr = contact_dist_f.ctypes.data_as(ctypes.c_void_p)

    # Pre-allocate result array using nframe from lazy STrajectories
    max_frames = trajs.nframe
    result_drms = np.zeros(max_frames, dtype=np.float64)
    result_ptr = result_drms.ctypes.data_as(ctypes.c_void_p)

    # Convert filename to C string
    dcd_filename_bytes = dcd_file.encode('utf-8')
    filename_len = len(dcd_filename_bytes)

    # Output variables
    nstru_out = ctypes.c_int()
    dcd_nframe_out = ctypes.c_int()
    dcd_natom_out = ctypes.c_int()

    with fortran_status() as (status, msg, msglen):
        lib.drms_analysis_lazy_c(
            dcd_filename_bytes,
            ctypes.c_int(filename_len),
            ctypes.c_int(trj_type),
            contact_list_ptr,
            contact_dist_ptr,
            ctypes.c_int(n_contact),
            ctypes.c_int(n_atoms),
            ctypes.c_int(ana_period),
            ctypes.c_int(1 if pbc_correct else 0),
            result_ptr,
            ctypes.c_int(max_frames),
            ctypes.byref(nstru_out),
            ctypes.byref(dcd_nframe_out),
            ctypes.byref(dcd_natom_out),
            ctypes.byref(status),
            msg,
            ctypes.c_int(msglen),
        )


    # Return DrmsAnalysisResult (same as non-lazy version)
    return DrmsAnalysisResult(result_drms[:nstru_out.value])


def drms_analysis(
        trajs: STrajectories,
        contact_list: np.ndarray,
        contact_dist: np.ndarray,
        ana_period: int = 1,
        pbc_correct: bool = False,
        ) -> DrmsAnalysisResult:
    """
    Executes distance RMSD (DRMS) analysis.

    This function calculates the root-mean-square deviation of distances
    between predefined atom contact pairs compared to reference distances.

    Args:
        trajs: Trajectories to analyze (memory or lazy)
        contact_list: Contact atom pairs as (2, n_contact) array with 1-indexed
                      atom indices. Each column is [atom1_idx, atom2_idx].
        contact_dist: Reference distances for each contact pair (n_contact,)
        ana_period: Analysis period (default: 1)
        pbc_correct: Apply PBC correction for distances (default: False)

    Returns:
        DrmsAnalysisResult containing the DRMS array

    Examples:
        >>> # Memory-based (standard)
        >>> contact_list = np.array([[1, 2, 3], [10, 11, 12]], dtype=np.int32)
        >>> contact_dist = np.array([5.0, 6.0, 7.0], dtype=np.float64)
        >>> result = drms_analysis(trajs, contact_list, contact_dist)
        >>> print(result.drms)

        >>> # Lazy mode: works with lazy trajectories from crd_convert(lazy=True)
        >>> lazy_trajs, mol = crd_convert(mol, ["traj.dcd"], lazy=True)
        >>> result = drms_analysis(lazy_trajs[0], contact_list, contact_dist)
    """
    # Handle lazy trajectories
    if trajs.is_lazy:
        return _drms_analysis_lazy(
            trajs=trajs,
            contact_list=contact_list,
            contact_dist=contact_dist,
            ana_period=ana_period,
            pbc_correct=pbc_correct,
        )

    lib = LibGenesis().lib

    # Ensure arrays are contiguous and correct dtype
    contact_list_f = np.asfortranarray(contact_list, dtype=np.int32)
    contact_dist_f = np.ascontiguousarray(contact_dist, dtype=np.float64)

    # Validate shapes
    if contact_list_f.ndim != 2 or contact_list_f.shape[0] != 2:
        raise GenesisValidationError(
            f"contact_list must be shape (2, n_contact), got {contact_list.shape}"
        )
    n_contact = contact_list_f.shape[1]
    if contact_dist_f.shape[0] != n_contact:
        raise GenesisValidationError(
            f"contact_dist must have {n_contact} elements, got {contact_dist.shape[0]}"
        )

    # Get pointers (zero-copy input)
    contact_list_ptr = contact_list_f.ctypes.data_as(ctypes.c_void_p)
    contact_dist_ptr = contact_dist_f.ctypes.data_as(ctypes.c_void_p)

    # Pre-allocate result array (full zero-copy)
    n_frame = int(trajs.nframe / ana_period)
    result_drms = np.zeros(n_frame, dtype=np.float64)
    result_ptr = result_drms.ctypes.data_as(ctypes.c_void_p)

    # Output variables
    nstru_out = ctypes.c_int()

    with fortran_status() as (status, msg, msglen):
        lib.drms_analysis_c(
            contact_list_ptr,
            contact_dist_ptr,
            ctypes.c_int(n_contact),
            ctypes.byref(trajs.get_c_obj()),
            ctypes.c_int(ana_period),
            ctypes.c_int(1 if pbc_correct else 0),
            result_ptr,
            ctypes.c_int(n_frame),
            ctypes.byref(nstru_out),
            ctypes.byref(status),
            msg,
            ctypes.c_int(msglen),
        )


    # Return result sliced to actual size
    return DrmsAnalysisResult(result_drms[:nstru_out.value])
