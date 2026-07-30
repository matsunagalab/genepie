"""Root mean square deviation, in memory and via lazy DCD reads."""

import ctypes
import os
from collections import namedtuple
import numpy as np
from typing import Optional

from ..libgenesis import LibGenesis
from ..s_molecule import SMolecule
from ..s_trajectories import STrajectories
from ..exceptions import GenesisValidationError
from .._fortran import fortran_status
from ._common import (
    _resolve_enum,
    _FITTING_METHOD_MAP,
    _prepare_lazy_trajectory,
)
from .converter import selection, crd_convert_info


RmsdAnalysisResult = namedtuple(
        'RmsdAnalysisResult',
        ['rmsd'])


def _rmsd_analysis_lazy(
        molecule: SMolecule,
        trajs: STrajectories,
        analysis_selection: str,
        fitting_selection: Optional[str] = None,
        fitting_method: str = "TR+ROT",
        ana_period: int = 1,
        mass_weighted: bool = False,
        ref_coord: Optional[np.ndarray] = None,
        ) -> "RmsdAnalysisResult":
    """
    Private implementation: RMSD analysis with lazy DCD loading.

    Called by rmsd_analysis() when trajs.is_lazy is True.
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

    # Ensure arrays are contiguous and correct dtype
    mass = np.ascontiguousarray(molecule.mass, dtype=np.float64)

    # Reference coordinates: Fortran expects (3, n_atoms)
    if ref_coord is None:
        ref_coord_arr = molecule.atom_coord
    else:
        ref_coord_arr = ref_coord
    ref_coord_f = np.asfortranarray(ref_coord_arr.T, dtype=np.float64)

    # Get pointers (zero-copy)
    mass_ptr = mass.ctypes.data_as(ctypes.c_void_p)
    ref_coord_ptr = ref_coord_f.ctypes.data_as(ctypes.c_void_p)

    # Pre-allocate exactly the number of analyzed frames.
    result_rmsd = np.zeros(n_result, dtype=np.float64)
    result_ptr = result_rmsd.ctypes.data_as(ctypes.c_void_p)
    filename_len = len(dcd_filename_bytes)

    # Output variables
    nstru_out = ctypes.c_int()
    dcd_nframe_out = ctypes.c_int()
    dcd_natom_out = ctypes.c_int()

    # Fitting parameters
    fitting_idx_ptr = ctypes.c_void_p(0)
    n_fitting = 0
    method_int = 0

    if fitting_selection is not None:
        method_int = _resolve_enum(
            fitting_method, _FITTING_METHOD_MAP, "fitting_method")

        # Get fitting indices
        fitting_indices = selection(molecule, fitting_selection)
        n_fitting = len(fitting_indices)
        fitting_idx_ptr = fitting_indices.ctypes.data_as(ctypes.c_void_p)

    with fortran_status() as (status, msg, msglen):
        lib.rmsd_analysis_lazy_c(
            dcd_filename_bytes,
            ctypes.c_int(filename_len),
            ctypes.c_int(trj_type),
            ctypes.c_int(trajs.lazy_dcd_natom),
            source_selection.ctypes.data_as(ctypes.c_void_p),
            ctypes.c_int(trajs.natom),
            mass_ptr,
            ref_coord_ptr,
            ctypes.c_int(molecule.num_atoms),
            ctypes.c_int(effective_period),
            fitting_idx_ptr,
            ctypes.c_int(n_fitting),
            analysis_indices.ctypes.data_as(ctypes.c_void_p),
            ctypes.c_int(n_analysis),
            ctypes.c_int(method_int),
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


    # Return RmsdAnalysisResult (same as non-lazy version)
    return RmsdAnalysisResult(result_rmsd[:nstru_out.value])


def rmsd_analysis(
        molecule: SMolecule,
        trajs: STrajectories,
        analysis_selection: str,
        fitting_selection: Optional[str] = None,
        fitting_method: str = "TR+ROT",
        ana_period: int = 1,
        mass_weighted: bool = False,
        ref_coord: Optional[np.ndarray] = None,
        ) -> RmsdAnalysisResult:
    """
    Executes RMSD analysis with optional structural fitting.

    This function calculates the root-mean-square deviation (RMSD) between
    trajectory frames and a reference structure. When fitting_selection is
    provided, structural alignment is performed before RMSD calculation.

    Args:
        molecule: Molecular structure (provides mass and reference coordinates)
        trajs: Trajectories to analyze
        analysis_selection: GENESIS selection for RMSD calculation atoms (e.g., "an:CA")
        fitting_selection: GENESIS selection for fitting atoms (e.g., "an:CA").
            If None, no fitting is performed (use when trajectory is pre-aligned).
        fitting_method: Fitting method (only used when fitting_selection is provided):
            - "NO": No fitting
            - "TR+ROT": Translation + rotation (default, most common)
            - "TR": Translation only
            - "TR+ZROT": Translation + Z-axis rotation
            - "XYTR": XY-plane translation only
            - "XYTR+ZROT": XY translation + Z-axis rotation
        ana_period: Analysis period (default: 1)
        mass_weighted: Use mass weighting for both fitting and RMSD (default: False)
        ref_coord: Reference coordinates (default: molecule.atom_coord).
                   Shape should be (n_atoms, 3).

    Returns:
        RmsdAnalysisResult containing the RMSD array

    Examples:
        >>> # Calculate RMSD without fitting (pre-aligned trajectory)
        >>> result = rmsd_analysis(mol, trajs, analysis_selection="an:CA")
        >>> print(result.rmsd)

        >>> # Calculate RMSD with TR+ROT fitting using CA atoms
        >>> result = rmsd_analysis(
        ...     mol, trajs,
        ...     analysis_selection="an:CA",
        ...     fitting_selection="an:CA",
        ...     fitting_method="TR+ROT"
        ... )
        >>> print(result.rmsd)

        >>> # Fit on CA atoms, but calculate RMSD for all heavy atoms
        >>> result = rmsd_analysis(
        ...     mol, trajs,
        ...     analysis_selection="heavy",
        ...     fitting_selection="an:CA",
        ...     fitting_method="TR+ROT"
        ... )

        >>> # Lazy mode: works with lazy trajectories from crd_convert(lazy=True)
        >>> lazy_trajs, mol = crd_convert(mol, ["traj.dcd"], lazy=True)
        >>> result = rmsd_analysis(mol, lazy_trajs[0], analysis_selection="an:CA")
    """
    # Handle lazy trajectories
    if trajs.is_lazy:
        return _rmsd_analysis_lazy(
            molecule=molecule,
            trajs=trajs,
            analysis_selection=analysis_selection,
            fitting_selection=fitting_selection,
            fitting_method=fitting_method,
            ana_period=ana_period,
            mass_weighted=mass_weighted,
            ref_coord=ref_coord,
        )

    lib = LibGenesis().lib

    # Get atom indices using GENESIS selection
    analysis_indices = selection(molecule, analysis_selection)
    n_analysis = len(analysis_indices)

    # Ensure arrays are contiguous and correct dtype
    mass = np.ascontiguousarray(molecule.mass, dtype=np.float64)

    # Reference coordinates: Fortran expects (3, n_atoms)
    if ref_coord is None:
        ref_coord_arr = molecule.atom_coord
    else:
        ref_coord_arr = ref_coord
    ref_coord_f = np.asfortranarray(ref_coord_arr.T, dtype=np.float64)

    # Get pointers (zero-copy)
    mass_ptr = mass.ctypes.data_as(ctypes.c_void_p)
    ref_coord_ptr = ref_coord_f.ctypes.data_as(ctypes.c_void_p)

    # Pre-allocate result array (full zero-copy)
    n_frame = int(trajs.nframe / ana_period)
    result_rmsd = np.zeros(n_frame, dtype=np.float64)
    result_ptr = result_rmsd.ctypes.data_as(ctypes.c_void_p)

    # Output variables
    nstru_out = ctypes.c_int()

    if fitting_selection is None:
        # No fitting - call simple RMSD analysis
        with fortran_status() as (status, msg, msglen):
            lib.rmsd_analysis_c(
                mass_ptr,
                ref_coord_ptr,
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
    else:
        # With fitting
        method_int = _resolve_enum(
            fitting_method, _FITTING_METHOD_MAP, "fitting_method")

        # Get fitting indices
        fitting_indices = selection(molecule, fitting_selection)
        n_fitting = len(fitting_indices)

        with fortran_status() as (status, msg, msglen):
            lib.rmsd_analysis_fitting_c(
                mass_ptr,
                ref_coord_ptr,
                ctypes.c_int(molecule.num_atoms),
                ctypes.byref(trajs.get_c_obj()),
                ctypes.c_int(ana_period),
                fitting_indices.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(n_fitting),
                analysis_indices.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(n_analysis),
                ctypes.c_int(method_int),
                ctypes.c_int(1 if mass_weighted else 0),
                result_ptr,
                ctypes.c_int(n_frame),
                ctypes.byref(nstru_out),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen),
            )


    # Return the pre-allocated result array (already filled by Fortran)
    return RmsdAnalysisResult(result_rmsd[:nstru_out.value])


RmsdLazyAnalysisResult = namedtuple(
        'RmsdLazyAnalysisResult',
        ['rmsd', 'dcd_nframe', 'dcd_natom'])


def rmsd_analysis_lazy(
        molecule: SMolecule,
        dcd_file: str,
        analysis_selection: str,
        fitting_selection: Optional[str] = None,
        fitting_method: str = "TR+ROT",
        ana_period: int = 1,
        mass_weighted: bool = False,
        ref_coord: Optional[np.ndarray] = None,
        has_box: bool = False,
        max_frames: int = 100000,
        ) -> RmsdLazyAnalysisResult:
    """
    Executes RMSD analysis with lazy DCD loading (memory efficient).

    This function calculates RMSD without loading the entire DCD trajectory
    into memory. Instead, frames are read on-demand directly from the DCD file.
    This is particularly useful for large trajectories that don't fit in memory.

    Args:
        molecule: Molecular structure (provides mass and reference coordinates)
        dcd_file: Path to DCD trajectory file
        analysis_selection: GENESIS selection for RMSD calculation atoms (e.g., "an:CA")
        fitting_selection: GENESIS selection for fitting atoms (e.g., "an:CA").
            If None, no fitting is performed (use when trajectory is pre-aligned).
        fitting_method: Fitting method (only used when fitting_selection is provided):
            - "NO": No fitting
            - "TR+ROT": Translation + rotation (default, most common)
            - "TR": Translation only
            - "TR+ZROT": Translation + Z-axis rotation
            - "XYTR": XY-plane translation only
            - "XYTR+ZROT": XY translation + Z-axis rotation
        ana_period: Analysis period (default: 1)
        mass_weighted: Use mass weighting for both fitting and RMSD (default: False)
        ref_coord: Reference coordinates (default: molecule.atom_coord).
                   Shape should be (n_atoms, 3).
        has_box: Whether DCD file contains box information (default: False)
        max_frames: Compatibility safety limit for analyzed frames. The result
            size is determined from the DCD header; an explicit smaller limit
            raises GenesisValidationError (default: 100000).

    Returns:
        RmsdLazyAnalysisResult containing:
            - rmsd: RMSD array
            - dcd_nframe: Total frames in DCD file
            - dcd_natom: Atoms per frame in DCD file

    Examples:
        >>> # Calculate RMSD without fitting (lazy loading)
        >>> result = rmsd_analysis_lazy(mol, "trajectory.dcd",
        ...                             analysis_selection="an:CA")
        >>> print(f"DCD has {result.dcd_nframe} frames, {result.dcd_natom} atoms")
        >>> print(result.rmsd)

        >>> # Calculate RMSD with TR+ROT fitting using CA atoms
        >>> result = rmsd_analysis_lazy(
        ...     mol, "large_trajectory.dcd",
        ...     analysis_selection="an:CA",
        ...     fitting_selection="an:CA",
        ...     fitting_method="TR+ROT"
        ... )
    """
    # Keep this legacy entry point as a thin wrapper around the canonical
    # lazy-trajectory API.
    if ana_period <= 0:
        raise GenesisValidationError(
            f"ana_period must be positive, got {ana_period}"
        )
    if not os.path.exists(dcd_file):
        raise GenesisValidationError(f"DCD file not found: {dcd_file}")
    trj_type = 2 if has_box else 1
    trj_type_name = "COOR+BOX" if has_box else "COOR"
    info = crd_convert_info(
        molecule, [dcd_file], trj_format="DCD", trj_type=trj_type_name
    )
    if not info.frame_counts or info.frame_counts[0] <= 0:
        raise GenesisValidationError("No frames found in DCD file")
    dcd_nframe = info.frame_counts[0]
    required_frames = dcd_nframe // ana_period
    if max_frames is not None and max_frames < required_frames:
        raise GenesisValidationError(
            f"max_frames={max_frames} is smaller than the required "
            f"{required_frames} frames"
        )

    lazy_traj = STrajectories.from_lazy(
        dcd_file=dcd_file,
        trj_type=trj_type,
        nframe=dcd_nframe,
        natom=molecule.num_atoms,
        ana_period=1,
    )
    result = rmsd_analysis(
        molecule=molecule,
        trajs=lazy_traj,
        analysis_selection=analysis_selection,
        fitting_selection=fitting_selection,
        fitting_method=fitting_method,
        ana_period=ana_period,
        mass_weighted=mass_weighted,
        ref_coord=ref_coord,
    )
    return RmsdLazyAnalysisResult(
        result.rmsd,
        dcd_nframe,
        molecule.num_atoms,
    )
