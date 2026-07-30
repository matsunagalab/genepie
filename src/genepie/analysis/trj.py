"""Distance, angle and dihedral measurements along a trajectory."""

import ctypes
from collections import namedtuple
import numpy as np
import numpy.typing as npt
from typing import (
    List,
    Optional,
    Tuple,
)

from ..libgenesis import LibGenesis
from ..s_molecule import SMolecule
from ..s_trajectories import STrajectories
from ..exceptions import GenesisValidationError
from .._fortran import fortran_status
from ._common import _prepare_lazy_trajectory


TrjAnalysisResult = namedtuple(
        'TrjAnalysisResult',
        ['distance',
         'angle',
         'torsion',
         'com_distance',
         'com_angle',
         'com_torsion'])


def _flatten_com_groups(
        groups: List[Tuple[List[int], ...]],
        n_per_measurement: int
) -> Tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.int32]]:
    """
    Flatten COM group data into arrays suitable for Fortran.

    Converts a list of group tuples (each containing atom index lists) into:
    - flat_atoms: all atom indices concatenated
    - offsets: starting position of each group in flat_atoms
    - indices: group indices for each measurement

    Args:
        groups: List of tuples, where each tuple contains atom index lists.
                For distances: [([atoms1], [atoms2]), ...]
                For angles: [([atoms1], [atoms2], [atoms3]), ...]
                For torsions: [([atoms1], [atoms2], [atoms3], [atoms4]), ...]
        n_per_measurement: Number of groups per measurement (2 for dist, 3 for angle, 4 for torsion)

    Returns:
        flat_atoms: 1D array of all atom indices (1-indexed for Fortran)
        offsets: 1D array of group offsets (0-indexed, length = n_groups + 1)
        indices: 1D array of group indices for each measurement (0-indexed)
    """
    all_atoms: List[int] = []
    offsets: List[int] = [0]
    group_indices: List[int] = []

    group_counter = 0
    for group_tuple in groups:
        for atom_list in group_tuple:
            all_atoms.extend(atom_list)
            offsets.append(len(all_atoms))
            group_indices.append(group_counter)
            group_counter += 1

    return (
        np.array(all_atoms, dtype=np.int32),
        np.array(offsets, dtype=np.int32),
        np.array(group_indices, dtype=np.int32)
    )


def _trj_analysis_lazy(
        trajs: STrajectories,
        distance_pairs: Optional[npt.NDArray[np.int32]] = None,
        angle_triplets: Optional[npt.NDArray[np.int32]] = None,
        torsion_quadruplets: Optional[npt.NDArray[np.int32]] = None,
        ana_period: int = 1,
        ) -> TrjAnalysisResult:
    """Private implementation: trj_analysis with lazy DCD loading.

    Called by trj_analysis() when trajs.is_lazy is True. Only atom-based
    distance/angle/torsion measurements are supported in lazy mode.
    """
    lib = LibGenesis().lib

    dcd_filename_bytes, source_selection, effective_period, n_frame = \
        _prepare_lazy_trajectory(trajs, ana_period)
    trj_type = trajs.lazy_trj_type
    n_atoms = trajs.natom

    # Prepare distance list
    n_dist = 0
    dist_list_ptr = ctypes.c_void_p()
    dist_f = None
    if distance_pairs is not None and len(distance_pairs) > 0:
        n_dist = distance_pairs.shape[0]
        dist_f = np.asfortranarray(distance_pairs.T, dtype=np.int32)
        dist_list_ptr = dist_f.ctypes.data_as(ctypes.c_void_p)

    # Prepare angle list
    n_angl = 0
    angl_list_ptr = ctypes.c_void_p()
    angl_f = None
    if angle_triplets is not None and len(angle_triplets) > 0:
        n_angl = angle_triplets.shape[0]
        angl_f = np.asfortranarray(angle_triplets.T, dtype=np.int32)
        angl_list_ptr = angl_f.ctypes.data_as(ctypes.c_void_p)

    # Prepare torsion list
    n_tors = 0
    tors_list_ptr = ctypes.c_void_p()
    tors_f = None
    if torsion_quadruplets is not None and len(torsion_quadruplets) > 0:
        n_tors = torsion_quadruplets.shape[0]
        tors_f = np.asfortranarray(torsion_quadruplets.T, dtype=np.int32)
        tors_list_ptr = tors_f.ctypes.data_as(ctypes.c_void_p)

    if n_dist == 0 and n_angl == 0 and n_tors == 0:
        raise GenesisValidationError(
            "trj_analysis requires at least one distance/angle/torsion "
            "measurement"
        )

    # Pre-allocate result arrays (zerocopy, sized to the analyzed frame count)
    result_distance = np.zeros((n_dist, n_frame), dtype=np.float64, order='F') if n_dist > 0 else None
    result_angle = np.zeros((n_angl, n_frame), dtype=np.float64, order='F') if n_angl > 0 else None
    result_torsion = np.zeros((n_tors, n_frame), dtype=np.float64, order='F') if n_tors > 0 else None

    dist_ptr = result_distance.ctypes.data_as(ctypes.c_void_p) if result_distance is not None else ctypes.c_void_p()
    angl_ptr = result_angle.ctypes.data_as(ctypes.c_void_p) if result_angle is not None else ctypes.c_void_p()
    tors_ptr = result_torsion.ctypes.data_as(ctypes.c_void_p) if result_torsion is not None else ctypes.c_void_p()

    nstru_out = ctypes.c_int()
    dcd_nframe_out = ctypes.c_int()
    dcd_natom_out = ctypes.c_int()

    with fortran_status() as (status, msg, msglen):
        lib.trj_analysis_lazy_c(
            dcd_filename_bytes,
            ctypes.c_int(len(dcd_filename_bytes)),
            ctypes.c_int(trj_type),
            ctypes.c_int(trajs.lazy_dcd_natom),
            source_selection.ctypes.data_as(ctypes.c_void_p),
            ctypes.c_int(n_atoms),
            dist_list_ptr,
            ctypes.c_int(n_dist),
            angl_list_ptr,
            ctypes.c_int(n_angl),
            tors_list_ptr,
            ctypes.c_int(n_tors),
            ctypes.c_int(n_atoms),
            ctypes.c_int(effective_period),
            ctypes.c_int(n_frame),
            dist_ptr,
            ctypes.c_int(n_dist * n_frame),
            angl_ptr,
            ctypes.c_int(n_angl * n_frame),
            tors_ptr,
            ctypes.c_int(n_tors * n_frame),
            ctypes.byref(nstru_out),
            ctypes.byref(dcd_nframe_out),
            ctypes.byref(dcd_natom_out),
            ctypes.byref(status),
            msg,
            ctypes.c_int(msglen),
        )

    n_actual = nstru_out.value

    # Transpose to Python convention and slice to analyzed frames
    final_distance = result_distance[:, :n_actual].T.copy() if result_distance is not None else None
    final_angle = result_angle[:, :n_actual].T.copy() if result_angle is not None else None
    final_torsion = result_torsion[:, :n_actual].T.copy() if result_torsion is not None else None

    return TrjAnalysisResult(
        final_distance, final_angle, final_torsion,
        None, None, None
    )


def trj_analysis(
        trajs: STrajectories,
        distance_pairs: Optional[npt.NDArray[np.int32]] = None,
        angle_triplets: Optional[npt.NDArray[np.int32]] = None,
        torsion_quadruplets: Optional[npt.NDArray[np.int32]] = None,
        cdis_groups: Optional[List[Tuple[List[int], List[int]]]] = None,
        cang_groups: Optional[List[Tuple[List[int], List[int], List[int]]]] = None,
        ctor_groups: Optional[List[Tuple[List[int], List[int], List[int], List[int]]]] = None,
        molecule: Optional[SMolecule] = None,
        ana_period: int = 1,
        ) -> TrjAnalysisResult:
    """
    Executes trajectory analysis for distances, angles, and torsions.

    This function calculates distances, angles, and dihedral angles from
    trajectory data. It supports both atom-based and COM-based measurements.

    Lazy trajectories (from ``crd_convert(..., lazy=True)``) are supported for
    atom-based distance/angle/torsion measurements; frames are read from the
    DCD file on demand. COM-based measurements still require an in-memory
    trajectory.

    Args:
        trajs: STrajectories object containing trajectory data (memory or lazy)
        distance_pairs: 2D array of shape (n_pairs, 2) with atom index pairs
                        (1-indexed as in Fortran convention)
        angle_triplets: 2D array of shape (n_triplets, 3) with atom indices
        torsion_quadruplets: 2D array of shape (n_quadruplets, 4) with atom indices
        cdis_groups: List of tuples for COM distance, each tuple contains two
                     lists of atom indices: [([atoms1], [atoms2]), ...]
                     (1-indexed as in Fortran convention)
        cang_groups: List of tuples for COM angles, each tuple contains three
                     lists of atom indices
        ctor_groups: List of tuples for COM torsions, each tuple contains four
                     lists of atom indices
        molecule: SMolecule object (required only for COM calculations)
        ana_period: Analysis period (default: 1)

    Returns:
        TrjAnalysisResult containing:
        - distance: 2D array of shape (n_frames, n_pairs) or None
        - angle: 2D array of shape (n_frames, n_triplets) or None
        - torsion: 2D array of shape (n_frames, n_quadruplets) or None
        - cdis: 2D array of shape (n_frames, n_cdis) or None
        - cang: 2D array of shape (n_frames, n_cang) or None
        - ctor: 2D array of shape (n_frames, n_ctor) or None

    Example:
        >>> # Atom-based distance
        >>> dist_pairs = np.array([[1, 2], [3, 4]], dtype=np.int32)
        >>> result = trj_analysis(trajs, distance_pairs=dist_pairs)
        >>> print(result.distance)

        >>> # COM-based distance (requires molecule)
        >>> cdis_groups = [([1, 2, 3], [4, 5, 6])]
        >>> result = trj_analysis(trajs, cdis_groups=cdis_groups, molecule=mol)
        >>> print(result.cdis)
    """
    lib = LibGenesis().lib

    # Check if COM groups are provided
    has_com = (cdis_groups is not None and len(cdis_groups) > 0) or \
              (cang_groups is not None and len(cang_groups) > 0) or \
              (ctor_groups is not None and len(ctor_groups) > 0)

    if has_com and molecule is None:
        raise ValueError("molecule is required for COM-based measurements")

    # Handle lazy trajectories (atom-based measurements only for now).
    if trajs.is_lazy:
        if has_com:
            raise GenesisValidationError(
                "COM-based trj_analysis is not supported for lazy "
                "trajectories yet; load the trajectory into memory instead"
            )
        return _trj_analysis_lazy(
            trajs,
            distance_pairs=distance_pairs,
            angle_triplets=angle_triplets,
            torsion_quadruplets=torsion_quadruplets,
            ana_period=ana_period,
        )

    n_frame = int(trajs.nframe / ana_period)

    # Prepare distance list
    n_dist = 0
    dist_list_ptr = ctypes.c_void_p()
    dist_f = None
    if distance_pairs is not None and len(distance_pairs) > 0:
        n_dist = distance_pairs.shape[0]
        dist_f = np.asfortranarray(distance_pairs.T, dtype=np.int32)
        dist_list_ptr = dist_f.ctypes.data_as(ctypes.c_void_p)

    # Prepare angle list
    n_angl = 0
    angl_list_ptr = ctypes.c_void_p()
    angl_f = None
    if angle_triplets is not None and len(angle_triplets) > 0:
        n_angl = angle_triplets.shape[0]
        angl_f = np.asfortranarray(angle_triplets.T, dtype=np.int32)
        angl_list_ptr = angl_f.ctypes.data_as(ctypes.c_void_p)

    # Prepare torsion list
    n_tors = 0
    tors_list_ptr = ctypes.c_void_p()
    tors_f = None
    if torsion_quadruplets is not None and len(torsion_quadruplets) > 0:
        n_tors = torsion_quadruplets.shape[0]
        tors_f = np.asfortranarray(torsion_quadruplets.T, dtype=np.int32)
        tors_list_ptr = tors_f.ctypes.data_as(ctypes.c_void_p)

    # Pre-allocate atom-based result arrays
    result_distance = np.zeros((n_dist, n_frame), dtype=np.float64, order='F') if n_dist > 0 else None
    result_angle = np.zeros((n_angl, n_frame), dtype=np.float64, order='F') if n_angl > 0 else None
    result_torsion = np.zeros((n_tors, n_frame), dtype=np.float64, order='F') if n_tors > 0 else None

    dist_ptr = result_distance.ctypes.data_as(ctypes.c_void_p) if result_distance is not None else ctypes.c_void_p()
    angl_ptr = result_angle.ctypes.data_as(ctypes.c_void_p) if result_angle is not None else ctypes.c_void_p()
    tors_ptr = result_torsion.ctypes.data_as(ctypes.c_void_p) if result_torsion is not None else ctypes.c_void_p()

    nstru_out = ctypes.c_int()

    if has_com:
        # COM-based analysis
        mass = np.ascontiguousarray(molecule.mass, dtype=np.float64)
        mass_ptr = mass.ctypes.data_as(ctypes.c_void_p)
        n_atoms = len(mass)

        # Prepare COM distance groups
        n_cdis = 0
        cdis_atoms = np.array([], dtype=np.int32)
        cdis_offsets = np.array([0], dtype=np.int32)
        cdis_pairs_arr = np.array([], dtype=np.int32)
        if cdis_groups is not None and len(cdis_groups) > 0:
            n_cdis = len(cdis_groups)
            cdis_atoms, cdis_offsets, cdis_pairs_arr = _flatten_com_groups(cdis_groups, 2)

        # Prepare COM angle groups
        n_cang = 0
        cang_atoms = np.array([], dtype=np.int32)
        cang_offsets = np.array([0], dtype=np.int32)
        cang_triplets_arr = np.array([], dtype=np.int32)
        if cang_groups is not None and len(cang_groups) > 0:
            n_cang = len(cang_groups)
            cang_atoms, cang_offsets, cang_triplets_arr = _flatten_com_groups(cang_groups, 3)

        # Prepare COM torsion groups
        n_ctor = 0
        ctor_atoms = np.array([], dtype=np.int32)
        ctor_offsets = np.array([0], dtype=np.int32)
        ctor_quads = np.array([], dtype=np.int32)
        if ctor_groups is not None and len(ctor_groups) > 0:
            n_ctor = len(ctor_groups)
            ctor_atoms, ctor_offsets, ctor_quads = _flatten_com_groups(ctor_groups, 4)

        # Pre-allocate COM result arrays
        result_cdis = np.zeros((n_cdis, n_frame), dtype=np.float64, order='F') if n_cdis > 0 else None
        result_cang = np.zeros((n_cang, n_frame), dtype=np.float64, order='F') if n_cang > 0 else None
        result_ctor = np.zeros((n_ctor, n_frame), dtype=np.float64, order='F') if n_ctor > 0 else None

        cdis_result_ptr = result_cdis.ctypes.data_as(ctypes.c_void_p) if result_cdis is not None else ctypes.c_void_p()
        cang_result_ptr = result_cang.ctypes.data_as(ctypes.c_void_p) if result_cang is not None else ctypes.c_void_p()
        ctor_result_ptr = result_ctor.ctypes.data_as(ctypes.c_void_p) if result_ctor is not None else ctypes.c_void_p()

        with fortran_status() as (status, msg, msglen):
            lib.trj_analysis_com_c(
                mass_ptr,
                ctypes.c_int(n_atoms),
                ctypes.byref(trajs.get_c_obj()),
                ctypes.c_int(ana_period),
                # Atom-based measurements
                dist_list_ptr,
                ctypes.c_int(n_dist),
                angl_list_ptr,
                ctypes.c_int(n_angl),
                tors_list_ptr,
                ctypes.c_int(n_tors),
                # COM distance
                cdis_atoms.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(len(cdis_atoms)),
                cdis_offsets.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(len(cdis_offsets)),
                cdis_pairs_arr.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(n_cdis),
                # COM angle
                cang_atoms.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(len(cang_atoms)),
                cang_offsets.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(len(cang_offsets)),
                cang_triplets_arr.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(n_cang),
                # COM torsion
                ctor_atoms.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(len(ctor_atoms)),
                ctor_offsets.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(len(ctor_offsets)),
                ctor_quads.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(n_ctor),
                # Pre-allocated output arrays
                dist_ptr,
                ctypes.c_int(n_dist * n_frame),
                angl_ptr,
                ctypes.c_int(n_angl * n_frame),
                tors_ptr,
                ctypes.c_int(n_tors * n_frame),
                cdis_result_ptr,
                ctypes.c_int(n_cdis * n_frame),
                cang_result_ptr,
                ctypes.c_int(n_cang * n_frame),
                ctor_result_ptr,
                ctypes.c_int(n_ctor * n_frame),
                # Output
                ctypes.byref(nstru_out),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen),
            )


        actual_frames = nstru_out.value

        # Transpose to Python convention and slice
        final_distance = result_distance.T[:actual_frames].copy() if result_distance is not None else None
        final_angle = result_angle.T[:actual_frames].copy() if result_angle is not None else None
        final_torsion = result_torsion.T[:actual_frames].copy() if result_torsion is not None else None
        final_cdis = result_cdis.T[:actual_frames].copy() if result_cdis is not None else None
        final_cang = result_cang.T[:actual_frames].copy() if result_cang is not None else None
        final_ctor = result_ctor.T[:actual_frames].copy() if result_ctor is not None else None

        return TrjAnalysisResult(
            final_distance, final_angle, final_torsion,
            final_cdis, final_cang, final_ctor
        )

    else:
        # Simple atom-based analysis (no COM)
        with fortran_status() as (status, msg, msglen):
            lib.trj_analysis_c(
                ctypes.byref(trajs.get_c_obj()),
                ctypes.c_int(ana_period),
                dist_list_ptr,
                ctypes.c_int(n_dist),
                angl_list_ptr,
                ctypes.c_int(n_angl),
                tors_list_ptr,
                ctypes.c_int(n_tors),
                dist_ptr,
                ctypes.c_int(n_dist * n_frame),
                angl_ptr,
                ctypes.c_int(n_angl * n_frame),
                tors_ptr,
                ctypes.c_int(n_tors * n_frame),
                ctypes.byref(nstru_out),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen),
            )


        n_actual = nstru_out.value

        # Transpose to Python convention and slice
        final_distance = result_distance[:, :n_actual].T.copy() if result_distance is not None else None
        final_angle = result_angle[:, :n_actual].T.copy() if result_angle is not None else None
        final_torsion = result_torsion[:, :n_actual].T.copy() if result_torsion is not None else None

        return TrjAnalysisResult(
            final_distance, final_angle, final_torsion,
            None, None, None
        )
