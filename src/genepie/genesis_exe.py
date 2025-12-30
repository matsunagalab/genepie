import ctypes
from collections import namedtuple
import io
import os
import tempfile
from typing import Iterable, List, NamedTuple, Optional, Tuple
import numpy as np
import numpy.typing as npt
from .libgenesis import LibGenesis
from .s_molecule import SMolecule
from .s_trajectories import STrajectories, STrajectoriesArray
from . import ctrl_files
from . import c2py_util
from . import py2c_util
from functools import lru_cache
from .exceptions import (
    GenesisFortranError,
    GenesisValidationError,
    raise_fortran_error,
)
from .output_capture import suppress_stdout_capture_stderr
from .validation import validate_positive, validate_non_negative, validate_trajectory_dimensions

_DEFAULT_MSG_LEN = 2048


def selection(molecule: SMolecule, selection_str: str) -> npt.NDArray[np.int32]:
    """Select atoms using GENESIS selection syntax.

    Args:
        molecule: SMolecule object containing molecular structure
        selection_str: GENESIS selection expression (e.g., "an:CA", "rn:ALA", "ri:1-10")

    Returns:
        numpy array of 1-indexed atom indices (Fortran convention)

    Raises:
        GenesisFortranError: If selection fails or returns no atoms

    Examples:
        >>> indices = selection(molecule, "an:CA")  # Select all CA atoms
        >>> indices = selection(molecule, "rn:ALA and an:CA")  # CA atoms in ALA residues
        >>> indices = selection(molecule, "ri:1-10")  # Atoms in residues 1-10
    """
    lib = LibGenesis().lib
    mol_c = molecule.to_SMoleculeC()

    # Prepare arguments
    sel_bytes = selection_str.encode('utf-8')
    indices_ptr = ctypes.c_void_p()
    n_indices = ctypes.c_int()
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    try:
        with suppress_stdout_capture_stderr() as captured:
            lib.selection_c(
                ctypes.byref(mol_c),
                sel_bytes,
                len(sel_bytes),
                ctypes.byref(indices_ptr),
                ctypes.byref(n_indices),
                ctypes.byref(status),
                msg,
                msglen
            )

        if status.value != 0:
            error_msg = msg.value.decode('utf-8', errors='replace').strip()
            raise_fortran_error(
                status.value,
                error_msg,
                stderr_output=captured.stderr
            )

        # Convert C pointer to numpy array (copy to avoid memory issues)
        if n_indices.value > 0:
            arr_ptr = ctypes.cast(indices_ptr, ctypes.POINTER(ctypes.c_int))
            indices = np.ctypeslib.as_array(arr_ptr, shape=(n_indices.value,)).copy()
        else:
            indices = np.array([], dtype=np.int32)

        return indices.astype(np.int32)

    finally:
        # Deallocate Fortran-allocated memory
        lib.deallocate_selection_c()
        # Deallocate molecule C structure
        lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


def crd_convert(
        molecule: SMolecule,
        traj_params: Optional[
            Iterable[ctrl_files.TrajectoryParameters]] = None,
        trj_format: Optional[str] = None,
        trj_type: Optional[str] = None,
        trj_natom: Optional[int] = None,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        fitting_method: Optional[str] = None,
        fitting_atom: Optional[int] = None,
        zrot_ngrid: Optional[int] = None,
        zrot_grid_size: Optional[float] = None,
        mass_weight: Optional[bool] = None,
        check_only: Optional[bool] = None,
        allow_backup: Optional[bool] = None,
        centering: Optional[bool] = None,
        centering_atom: Optional[int] = None,
        center_coord: Optional[tuple[float, float, float]] = None,
        pbc_correct: Optional[str] = None,
        rename_res: Optional[Iterable[str]] = None,
        ) -> tuple[STrajectoriesArray, SMolecule]:
    """
    Executes crd_convert.

    Args:
        molecule: SMolecule object containing molecular structure
        traj_params: List of trajectory parameters
        trj_format: Trajectory format (e.g., "DCD")
        trj_type: Trajectory type (e.g., "COOR")
        trj_natom: Number of atoms in trajectory
        selection_group: List of atom selection groups
        selection_mole_name: List of molecule names for selection
        fitting_method: Fitting method (e.g., "TR+ROT")
        fitting_atom: Fitting atom selection
        zrot_ngrid: Z-rotation grid size
        zrot_grid_size: Z-rotation grid size
        mass_weight: Whether to use mass weighting
        check_only: Whether to only check parameters
        allow_backup: Whether to allow backup files
        centering: Whether to center coordinates
        centering_atom: Centering atom selection
        center_coord: Center coordinates
        pbc_correct: PBC correction method
        rename_res: List of residue names to rename

    Returns:
        Tuple of (STrajectoriesArray, SMolecule) where the SMolecule contains
        only the atoms selected by selection_group, or the original molecule
        if no selection_group is specified.

    Notes:
        When md_step is None in TrajectoryParameters, frame count is
        auto-detected from the DCD header. This only works for DCD format.
    """
    # Warn if auto-detection is used with non-DCD format
    if traj_params is not None:
        import warnings
        for traj in traj_params:
            if traj.md_step is None and trj_format is not None:
                fmt_upper = trj_format.upper()
                if fmt_upper != "DCD":
                    warnings.warn(
                        f"Auto-detection of frame count (md_step=None) is only "
                        f"supported for DCD format, not {fmt_upper}. "
                        f"Please specify md_step explicitly.",
                        UserWarning
                    )

    buf = ctypes.c_void_p(None)
    num_trajs_c = ctypes.c_int(0)
    selected_atom_indices_c = ctypes.c_void_p(None)
    num_selected_atoms_c = ctypes.c_int(0)
    mol_c = molecule.to_SMoleculeC()

    try:
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_output(
                ctrl,
                trjfile="dummy.trj",
                )
        ctrl_files.write_trajectory_info(
                ctrl, traj_params, trj_format, trj_type, trj_natom,)
        # Use "all" as default selection when no selection is specified
        # This ensures that group1 is available for fitting and output
        selection_group_to_use = selection_group if selection_group else ["all"]
        ctrl_files.write_ctrl_selection(
                ctrl, selection_group_to_use, selection_mole_name)
        # Only specify fitting_atom if selection_group is provided
        # When no selection is specified, don't specify fitting_atom
        fitting_atom_to_use = fitting_atom if selection_group else None
        ctrl_files.write_ctrl_fitting(
                ctrl, fitting_method, fitting_atom_to_use,
                zrot_ngrid, zrot_grid_size, mass_weight)
        ctrl.write(b"[OPTION]\n")
        ctrl_files.write_kwargs(
                ctrl,
                check_only=check_only,
                allow_backup=allow_backup,
                centering=centering,
                centering_atom=centering_atom,
                center_coord=center_coord,
                pbc_correct=pbc_correct,
                rename_res=ctrl_files.NumberingData(rename_res),
                )

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)

        msgbuf, MSG_LEN = make_msgbuf()
        status = ctypes.c_int(0)

        with suppress_stdout_capture_stderr() as captured:
            LibGenesis().lib.crd_convert_c(
                    ctypes.byref(mol_c),
                    ctrl_bytes,
                    ctypes.c_int(ctrl_len),
                    ctypes.byref(buf),
                    ctypes.byref(num_trajs_c),
                    ctypes.byref(selected_atom_indices_c),
                    ctypes.byref(num_selected_atoms_c),
                    ctypes.byref(status),
                    msgbuf,
                    ctypes.c_int(MSG_LEN),
                    )
        if status.value != 0:
            error_msg = msgbuf.value.decode("utf-8", "replace")
            raise_fortran_error(
                error_msg,
                code=status.value,
                stderr_output=captured.stderr
            )

    finally:
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))
    
    # Create subset molecule based on selected atom indices
    if num_selected_atoms_c.value > 0 and selected_atom_indices_c:
        atom_indices = c2py_util.conv_int_ndarray(
            selected_atom_indices_c, num_selected_atoms_c.value)
        # Convert from 1-based (FORTRAN) to 0-based (Python) indexing
        subset_mol = molecule.subset_atoms(atom_indices - 1)
    else:
        # No selection or no atoms selected, return original molecule
        subset_mol = molecule
    
    return (STrajectoriesArray(buf, num_trajs_c.value), subset_mol)


TrjAnalysisResult = namedtuple(
        'TrjAnalysisResult',
        ['distance',
         'angle',
         'torsion',
         'com_distance',
         'com_angle',
         'com_torsion'])


def trj_analysis(molecule: SMolecule, trajs: STrajectories,
                 ana_period: Optional[int] = 1,
                 selection_group: Optional[Iterable[str]] = None,
                 selection_mole_name: Optional[Iterable[str]] = None,
                 check_only: Optional[bool] = None,
                 distance: Optional[Iterable[str]] = None,
                 dist_weight: Optional[Iterable[str]] = None,
                 angle: Optional[Iterable[str]] = None,
                 torsion: Optional[Iterable[str]] = None,
                 com_distance: Optional[Iterable[str]] = None,
                 com_angle: Optional[Iterable[str]] = None,
                 com_torsion: Optional[Iterable[str]] = None,
                 ) -> TrjAnalysisResult:
    """
    Executes trj_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
        selection_group:
        selection_mole_name:
    Returns:
        (distance, angle, torsion, cdis, cang, ctor)
    """
    num_distance = ctypes.c_int(0)
    result_distance_c = ctypes.c_void_p(None)
    num_angle = ctypes.c_int(0)
    result_angle_c = ctypes.c_void_p(None)
    num_torsion = ctypes.c_int(0)
    result_torsion_c = ctypes.c_void_p(None)
    num_cdis = ctypes.c_int(0)
    result_cdis_c = ctypes.c_void_p(None)
    num_cang = ctypes.c_int(0)
    result_cang_c = ctypes.c_void_p(None)
    num_ctor = ctypes.c_int(0)
    result_ctor_c = ctypes.c_void_p(None)
    ana_period_c = ctypes.c_int(ana_period)
    mol_c = ctypes.c_void_p(None)
    try:
        mol_c = molecule.to_SMoleculeC()
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_output(
                ctrl,
                disfile="dummy.dis",
                angfile="dummy.ang",
                torfile="dummy.tor",
                comdisfile="dummy.comdis",
                comangfile="dummy.comangr",
                comtorfile="dummy.comtor",
                qntfile="dummy.qnt",
                )
        ctrl_files.write_ctrl_selection(
                ctrl, selection_group, selection_mole_name)
        ctrl.write(b"[OPTION]\n")
        ctrl_files.write_kwargs(
                ctrl,
                check_only=check_only,
                distance=ctrl_files.NumberingData(distance),
                dist_weight=ctrl_files.NumberingData(dist_weight),
                angle=ctrl_files.NumberingData(angle),
                torsion=ctrl_files.NumberingData(torsion),
                com_distance=ctrl_files.NumberingData(com_distance),
                com_angle=ctrl_files.NumberingData(com_angle),
                com_torsion=ctrl_files.NumberingData(com_torsion),
                )

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)
        with suppress_stdout_capture_stderr() as captured:
            LibGenesis().lib.trj_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    ctrl_bytes,
                    ctypes.c_int(ctrl_len),
                    ctypes.byref(result_distance_c),
                    ctypes.byref(num_distance),
                    ctypes.byref(result_angle_c),
                    ctypes.byref(num_angle),
                    ctypes.byref(result_torsion_c),
                    ctypes.byref(num_torsion),
                    ctypes.byref(result_cdis_c),
                    ctypes.byref(num_cdis),
                    ctypes.byref(result_cang_c),
                    ctypes.byref(num_cang),
                    ctypes.byref(result_ctor_c),
                    ctypes.byref(num_ctor),
                    )
        n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
        result_distance = (c2py_util.conv_double_ndarray(
            result_distance_c, [n_frame_c.value, num_distance.value])
                           if result_distance_c else None)
        result_angle = (c2py_util.conv_double_ndarray(
                result_angle_c, [n_frame_c.value, num_angle.value])
                        if result_angle_c else None)
        result_torsion = (c2py_util.conv_double_ndarray(
                result_torsion_c, [n_frame_c.value, num_torsion.value])
                        if result_torsion_c else None)
        result_cdis = (c2py_util.conv_double_ndarray(
            result_cdis_c, [n_frame_c.value, num_cdis.value])
                           if result_cdis_c else None)
        result_cang = (c2py_util.conv_double_ndarray(
                result_cang_c, [n_frame_c.value, num_cang.value])
                        if result_cang_c else None)
        result_ctor = (c2py_util.conv_double_ndarray(
                result_ctor_c, [n_frame_c.value, num_ctor.value])
                        if result_ctor_c else None)
        return TrjAnalysisResult(
                result_distance, result_angle, result_torsion,
                result_cdis, result_cang, result_ctor)
    finally:
        if result_ctor_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_ctor_c),
                    ctypes.byref(n_frame_c), ctypes.byref(num_ctor))
        if result_cang_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_cang_c),
                    ctypes.byref(n_frame_c), ctypes.byref(num_cang))
        if result_cdis_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_cdis_c),
                    ctypes.byref(n_frame_c), ctypes.byref(num_cdis))
        if result_torsion_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_torsion_c),
                    ctypes.byref(n_frame_c), ctypes.byref(num_torsion))
        if result_angle_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_angle_c),
                    ctypes.byref(n_frame_c), ctypes.byref(num_angle))
        if result_distance_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_distance_c),
                    ctypes.byref(n_frame_c), ctypes.byref(num_distance))
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


TrjAnalysisZerocopyResult = namedtuple(
        'TrjAnalysisZerocopyResult',
        ['distance', 'angle', 'torsion'])


def trj_analysis_zerocopy(
        trajs: STrajectories,
        distance_pairs: Optional[npt.NDArray[np.int32]] = None,
        angle_triplets: Optional[npt.NDArray[np.int32]] = None,
        torsion_quadruplets: Optional[npt.NDArray[np.int32]] = None,
        ana_period: int = 1,
        ) -> TrjAnalysisZerocopyResult:
    """
    Executes trj_analysis with zerocopy interface.

    This simplified version computes distances, angles, and torsions using
    atom indices passed directly from Python. It does not support COM-based
    measurements (use the legacy trj_analysis for those).

    Args:
        trajs: STrajectories object containing trajectory data
        distance_pairs: 2D array of shape (n_pairs, 2) with atom index pairs
                        (1-indexed as in Fortran convention)
        angle_triplets: 2D array of shape (n_triplets, 3) with atom indices
                        (1-indexed as in Fortran convention)
        torsion_quadruplets: 2D array of shape (n_quadruplets, 4) with atom indices
                             (1-indexed as in Fortran convention)
        ana_period: Analysis period (default: 1)

    Returns:
        TrjAnalysisZerocopyResult containing:
        - distance: 2D array of shape (n_frames, n_pairs)
        - angle: 2D array of shape (n_frames, n_triplets)
        - torsion: 2D array of shape (n_frames, n_quadruplets)

    Example:
        >>> # Compute distance between atoms 1-2 and 3-4
        >>> dist_pairs = np.array([[1, 2], [3, 4]], dtype=np.int32)
        >>> result = trj_analysis_zerocopy(trajs, distance_pairs=dist_pairs)
        >>> print(result.distance)
    """
    lib = LibGenesis().lib

    # Prepare distance list
    n_dist = 0
    dist_ptr = ctypes.c_void_p()
    if distance_pairs is not None and len(distance_pairs) > 0:
        n_dist = distance_pairs.shape[0]
        # Transpose to (2, n_dist) in Fortran order
        dist_f = np.asfortranarray(distance_pairs.T, dtype=np.int32)
        dist_ptr = dist_f.ctypes.data_as(ctypes.c_void_p)

    # Prepare angle list
    n_angl = 0
    angl_ptr = ctypes.c_void_p()
    if angle_triplets is not None and len(angle_triplets) > 0:
        n_angl = angle_triplets.shape[0]
        # Transpose to (3, n_angl) in Fortran order
        angl_f = np.asfortranarray(angle_triplets.T, dtype=np.int32)
        angl_ptr = angl_f.ctypes.data_as(ctypes.c_void_p)

    # Prepare torsion list
    n_tors = 0
    tors_ptr = ctypes.c_void_p()
    if torsion_quadruplets is not None and len(torsion_quadruplets) > 0:
        n_tors = torsion_quadruplets.shape[0]
        # Transpose to (4, n_tors) in Fortran order
        tors_f = np.asfortranarray(torsion_quadruplets.T, dtype=np.int32)
        tors_ptr = tors_f.ctypes.data_as(ctypes.c_void_p)

    # Output variables
    result_distance_ptr = ctypes.c_void_p()
    result_angle_ptr = ctypes.c_void_p()
    result_torsion_ptr = ctypes.c_void_p()
    n_frames = ctypes.c_int()
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    try:
        with suppress_stdout_capture_stderr() as captured:
            lib.trj_analysis_zerocopy_c(
                ctypes.byref(trajs.get_c_obj()),
                ctypes.c_int(ana_period),
                dist_ptr,
                ctypes.c_int(n_dist),
                angl_ptr,
                ctypes.c_int(n_angl),
                tors_ptr,
                ctypes.c_int(n_tors),
                ctypes.byref(result_distance_ptr),
                ctypes.byref(result_angle_ptr),
                ctypes.byref(result_torsion_ptr),
                ctypes.byref(n_frames),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen),
            )

        # Check for errors
        if status.value != 0:
            error_msg = msg.value.decode('utf-8', errors='replace').strip()
            stderr_output = captured.stderr if captured else ""
            raise_fortran_error(status.value, error_msg, stderr_output)

        n_frame = n_frames.value

        # Convert results to numpy arrays (handle Fortran column-major order)
        result_distance = None
        if n_dist > 0 and result_distance_ptr:
            dist_flat = c2py_util.conv_double_ndarray(result_distance_ptr, n_dist * n_frame)
            result_distance = dist_flat.reshape((n_dist, n_frame), order='F').T

        result_angle = None
        if n_angl > 0 and result_angle_ptr:
            angl_flat = c2py_util.conv_double_ndarray(result_angle_ptr, n_angl * n_frame)
            result_angle = angl_flat.reshape((n_angl, n_frame), order='F').T

        result_torsion = None
        if n_tors > 0 and result_torsion_ptr:
            tors_flat = c2py_util.conv_double_ndarray(result_torsion_ptr, n_tors * n_frame)
            result_torsion = tors_flat.reshape((n_tors, n_frame), order='F').T

        return TrjAnalysisZerocopyResult(result_distance, result_angle, result_torsion)

    finally:
        # Cleanup - deallocate Fortran allocated results
        lib.deallocate_trj_results_c()


TrjAnalysisZerocopyCOMResult = namedtuple(
        'TrjAnalysisZerocopyCOMResult',
        ['distance', 'angle', 'torsion', 'cdis', 'cang', 'ctor'])


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


def trj_analysis_zerocopy_com(
        molecule: SMolecule,
        trajs: STrajectories,
        distance_pairs: Optional[npt.NDArray[np.int32]] = None,
        angle_triplets: Optional[npt.NDArray[np.int32]] = None,
        torsion_quadruplets: Optional[npt.NDArray[np.int32]] = None,
        cdis_groups: Optional[List[Tuple[List[int], List[int]]]] = None,
        cang_groups: Optional[List[Tuple[List[int], List[int], List[int]]]] = None,
        ctor_groups: Optional[List[Tuple[List[int], List[int], List[int], List[int]]]] = None,
        ana_period: int = 1,
        ) -> TrjAnalysisZerocopyCOMResult:
    """
    Executes trj_analysis with zerocopy interface including COM calculations.

    This version supports both atom-based and COM-based measurements.
    COM calculations use the flat array + offset approach for efficient
    variable-length group handling.

    Args:
        molecule: SMolecule object containing mass information
        trajs: STrajectories object containing trajectory data
        distance_pairs: 2D array of shape (n_pairs, 2) with atom index pairs
                        (1-indexed as in Fortran convention)
        angle_triplets: 2D array of shape (n_triplets, 3) with atom indices
                        (1-indexed as in Fortran convention)
        torsion_quadruplets: 2D array of shape (n_quadruplets, 4) with atom indices
                             (1-indexed as in Fortran convention)
        cdis_groups: List of tuples for COM distance, each tuple contains two
                     lists of atom indices: [([atoms1], [atoms2]), ...]
                     (1-indexed as in Fortran convention)
        cang_groups: List of tuples for COM angles, each tuple contains three
                     lists of atom indices
        ctor_groups: List of tuples for COM torsions, each tuple contains four
                     lists of atom indices
        ana_period: Analysis period (default: 1)

    Returns:
        TrjAnalysisZerocopyCOMResult containing:
        - distance: 2D array of shape (n_frames, n_pairs)
        - angle: 2D array of shape (n_frames, n_triplets)
        - torsion: 2D array of shape (n_frames, n_quadruplets)
        - cdis: 2D array of shape (n_frames, n_cdis)
        - cang: 2D array of shape (n_frames, n_cang)
        - ctor: 2D array of shape (n_frames, n_ctor)

    Example:
        >>> # Compute COM distance between residue 1 and residue 2
        >>> # Residue 1 atoms: indices 1-4, Residue 2 atoms: indices 5-8
        >>> cdis_groups = [([1, 2, 3, 4], [5, 6, 7, 8])]
        >>> result = trj_analysis_zerocopy_com(mol, trajs, cdis_groups=cdis_groups)
        >>> print(result.cdis)
    """
    lib = LibGenesis().lib

    # Get mass array (zerocopy)
    mass = molecule.mass.astype(np.float64, copy=False)
    mass = np.ascontiguousarray(mass)
    mass_ptr = mass.ctypes.data_as(ctypes.c_void_p)
    n_atoms = len(mass)

    # Prepare distance list
    n_dist = 0
    dist_ptr = ctypes.c_void_p()
    if distance_pairs is not None and len(distance_pairs) > 0:
        n_dist = distance_pairs.shape[0]
        dist_f = np.asfortranarray(distance_pairs.T, dtype=np.int32)
        dist_ptr = dist_f.ctypes.data_as(ctypes.c_void_p)

    # Prepare angle list
    n_angl = 0
    angl_ptr = ctypes.c_void_p()
    if angle_triplets is not None and len(angle_triplets) > 0:
        n_angl = angle_triplets.shape[0]
        angl_f = np.asfortranarray(angle_triplets.T, dtype=np.int32)
        angl_ptr = angl_f.ctypes.data_as(ctypes.c_void_p)

    # Prepare torsion list
    n_tors = 0
    tors_ptr = ctypes.c_void_p()
    if torsion_quadruplets is not None and len(torsion_quadruplets) > 0:
        n_tors = torsion_quadruplets.shape[0]
        tors_f = np.asfortranarray(torsion_quadruplets.T, dtype=np.int32)
        tors_ptr = tors_f.ctypes.data_as(ctypes.c_void_p)

    # Prepare COM distance groups
    n_cdis = 0
    cdis_atoms = np.array([], dtype=np.int32)
    cdis_offsets = np.array([0], dtype=np.int32)
    cdis_pairs = np.array([], dtype=np.int32)
    if cdis_groups is not None and len(cdis_groups) > 0:
        n_cdis = len(cdis_groups)
        cdis_atoms, cdis_offsets, cdis_pairs = _flatten_com_groups(cdis_groups, 2)

    # Prepare COM angle groups
    n_cang = 0
    cang_atoms = np.array([], dtype=np.int32)
    cang_offsets = np.array([0], dtype=np.int32)
    cang_triplets = np.array([], dtype=np.int32)
    if cang_groups is not None and len(cang_groups) > 0:
        n_cang = len(cang_groups)
        cang_atoms, cang_offsets, cang_triplets = _flatten_com_groups(cang_groups, 3)

    # Prepare COM torsion groups
    n_ctor = 0
    ctor_atoms = np.array([], dtype=np.int32)
    ctor_offsets = np.array([0], dtype=np.int32)
    ctor_quads = np.array([], dtype=np.int32)
    if ctor_groups is not None and len(ctor_groups) > 0:
        n_ctor = len(ctor_groups)
        ctor_atoms, ctor_offsets, ctor_quads = _flatten_com_groups(ctor_groups, 4)

    # Output variables
    result_distance_ptr = ctypes.c_void_p()
    result_angle_ptr = ctypes.c_void_p()
    result_torsion_ptr = ctypes.c_void_p()
    result_cdis_ptr = ctypes.c_void_p()
    result_cang_ptr = ctypes.c_void_p()
    result_ctor_ptr = ctypes.c_void_p()
    n_frames = ctypes.c_int()
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    try:
        with suppress_stdout_capture_stderr() as captured:
            lib.trj_analysis_zerocopy_com_c(
                mass_ptr,
                ctypes.c_int(n_atoms),
                ctypes.byref(trajs.get_c_obj()),
                ctypes.c_int(ana_period),
                # Atom-based measurements
                dist_ptr,
                ctypes.c_int(n_dist),
                angl_ptr,
                ctypes.c_int(n_angl),
                tors_ptr,
                ctypes.c_int(n_tors),
                # COM distance
                cdis_atoms.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(len(cdis_atoms)),
                cdis_offsets.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(len(cdis_offsets)),
                cdis_pairs.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(n_cdis),
                # COM angle
                cang_atoms.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(len(cang_atoms)),
                cang_offsets.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(len(cang_offsets)),
                cang_triplets.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(n_cang),
                # COM torsion
                ctor_atoms.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(len(ctor_atoms)),
                ctor_offsets.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(len(ctor_offsets)),
                ctor_quads.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(n_ctor),
                # Output pointers
                ctypes.byref(result_distance_ptr),
                ctypes.byref(result_angle_ptr),
                ctypes.byref(result_torsion_ptr),
                ctypes.byref(result_cdis_ptr),
                ctypes.byref(result_cang_ptr),
                ctypes.byref(result_ctor_ptr),
                ctypes.byref(n_frames),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen),
            )

        # Check for errors
        if status.value != 0:
            error_msg = msg.value.decode('utf-8', errors='replace').strip()
            stderr_output = captured.stderr if captured else ""
            raise_fortran_error(status.value, error_msg, stderr_output)

        n_frame = n_frames.value

        # Convert atom-based results
        result_distance = None
        if n_dist > 0 and result_distance_ptr:
            dist_flat = c2py_util.conv_double_ndarray(result_distance_ptr, n_dist * n_frame)
            result_distance = dist_flat.reshape((n_dist, n_frame), order='F').T

        result_angle = None
        if n_angl > 0 and result_angle_ptr:
            angl_flat = c2py_util.conv_double_ndarray(result_angle_ptr, n_angl * n_frame)
            result_angle = angl_flat.reshape((n_angl, n_frame), order='F').T

        result_torsion = None
        if n_tors > 0 and result_torsion_ptr:
            tors_flat = c2py_util.conv_double_ndarray(result_torsion_ptr, n_tors * n_frame)
            result_torsion = tors_flat.reshape((n_tors, n_frame), order='F').T

        # Convert COM-based results
        result_cdis = None
        if n_cdis > 0 and result_cdis_ptr:
            cdis_flat = c2py_util.conv_double_ndarray(result_cdis_ptr, n_cdis * n_frame)
            result_cdis = cdis_flat.reshape((n_cdis, n_frame), order='F').T

        result_cang = None
        if n_cang > 0 and result_cang_ptr:
            cang_flat = c2py_util.conv_double_ndarray(result_cang_ptr, n_cang * n_frame)
            result_cang = cang_flat.reshape((n_cang, n_frame), order='F').T

        result_ctor = None
        if n_ctor > 0 and result_ctor_ptr:
            ctor_flat = c2py_util.conv_double_ndarray(result_ctor_ptr, n_ctor * n_frame)
            result_ctor = ctor_flat.reshape((n_ctor, n_frame), order='F').T

        return TrjAnalysisZerocopyCOMResult(
            result_distance, result_angle, result_torsion,
            result_cdis, result_cang, result_ctor
        )

    finally:
        # Cleanup - deallocate Fortran allocated results
        lib.deallocate_trj_results_c()


RgAnalysisResult = namedtuple(
        'RgAnalysisResult',
        ['rg'])


def rg_analysis(molecule: SMolecule, trajs: STrajectories,
                ana_period: Optional[int] = 1,
                selection_group: Optional[Iterable[str]] = None,
                selection_mole_name: Optional[Iterable[str]] = None,
                fitting_method: Optional[str] = None,
                fitting_atom: Optional[int] = None,
                check_only: Optional[bool] = None,
                analysis_atom: Optional[int] = None,
                mass_weighted: Optional[bool] = None,
                ) -> RgAnalysisResult:
    """
    Executes rg_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
    Returns:
        rg
    """
    ana_period_c = ctypes.c_int(ana_period)
    result_rg_c = ctypes.c_void_p()
    mol_c = molecule.to_SMoleculeC()
    try:
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_output(
                ctrl,
                rgfile="dummy.rg")
        ctrl_files.write_ctrl_selection(
                ctrl, selection_group, selection_mole_name)
        ctrl_files.write_ctrl_fitting(
                ctrl, fitting_method, fitting_atom)
        ctrl.write(b"[OPTION]\n")
        ctrl_files.write_kwargs(
                ctrl,
                check_only=check_only,
                analysis_atom=analysis_atom,
                mass_weighted=mass_weighted,
                )

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)
        with suppress_stdout_capture_stderr() as captured:
            LibGenesis().lib.rg_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    ctrl_bytes,
                    ctypes.c_int(ctrl_len),
                    ctypes.byref(result_rg_c),
                    )
        n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
        result_rg = (c2py_util.conv_double_ndarray(
            result_rg_c, n_frame_c.value)
                    if result_rg_c else None)
        return RgAnalysisResult(
                result_rg)
    finally:
        if result_rg_c:
            LibGenesis().lib.deallocate_double(
                    ctypes.byref(result_rg_c),
                    ctypes.byref(n_frame_c))
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


def rg_analysis_zerocopy(
        molecule: SMolecule,
        trajs: STrajectories,
        analysis_selection: str,
        ana_period: int = 1,
        mass_weighted: bool = True,
        ) -> RgAnalysisResult:
    """
    Executes rg_analysis with true zero-copy interface.

    This function passes the mass array pointer directly from Python NumPy
    to Fortran, eliminating all data copies. The Fortran code creates a
    view of the Python array using C_F_POINTER.

    Args:
        molecule: Molecular structure
        trajs: Trajectories to analyze
        analysis_selection: GENESIS selection string (e.g., "an:CA", "heavy")
        ana_period: Analysis period (default: 1)
        mass_weighted: Use mass weighting for RG calculation (default: True)

    Returns:
        RgAnalysisResult containing the radius of gyration array

    Example:
        >>> result = rg_analysis_zerocopy(mol, trajs, "an:CA")
        >>> print(result.rg)
    """
    lib = LibGenesis().lib

    # Get atom indices using GENESIS selection
    analysis_indices = selection(molecule, analysis_selection)
    n_analysis = len(analysis_indices)

    # Ensure mass array is contiguous and correct dtype
    mass = np.ascontiguousarray(molecule.mass, dtype=np.float64)

    # Get pointer to mass array (zero-copy)
    mass_ptr = mass.ctypes.data_as(ctypes.c_void_p)

    # Output variables
    result_rg_c = ctypes.c_void_p()
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    try:
        with suppress_stdout_capture_stderr() as captured:
            lib.rg_analysis_zerocopy_c(
                mass_ptr,
                ctypes.c_int(molecule.num_atoms),
                ctypes.byref(trajs.get_c_obj()),
                ctypes.c_int(ana_period),
                analysis_indices.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(n_analysis),
                ctypes.c_int(1 if mass_weighted else 0),
                ctypes.byref(result_rg_c),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen),
            )

        # Check for errors
        if status.value != 0:
            error_msg = msg.value.decode('utf-8', errors='replace').strip()
            stderr_output = captured.stderr if captured else ""
            raise_fortran_error(status.value, error_msg, stderr_output)

        # Convert result to numpy array
        n_frame = int(trajs.nframe / ana_period)
        result_rg = (c2py_util.conv_double_ndarray(result_rg_c, n_frame)
                     if result_rg_c else None)

        return RgAnalysisResult(result_rg)

    finally:
        # Cleanup - only deallocate results, not the mass array (Python owns it)
        lib.deallocate_rg_results_c()
        # Note: NO deallocate_s_molecule_c call - we didn't allocate any s_molecule_c!


RmsdAnalysisResult = namedtuple(
        'RmsdAnalysisResult',
        ['rmsd'])


def rmsd_analysis(
        molecule: SMolecule, trajs: STrajectories,
        ana_period: Optional[int] = 1,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        fitting_method: Optional[str] = None,
        fitting_atom: Optional[int] = None,
        check_only: Optional[bool] = None,
        analysis_atom: Optional[int] = None,
        ) -> RmsdAnalysisResult:
    """
    Executes rmsd_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
    Returns:
        rmsd
    """
    mol_c = None
    ana_period_c = ctypes.c_int(ana_period)
    result_rmsd_c = ctypes.c_void_p()
    try:
        mol_c = molecule.to_SMoleculeC()
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_output(
                ctrl,
                rmsfile="dummy.rms")
        ctrl_files.write_ctrl_selection(
                ctrl, selection_group, selection_mole_name)
        ctrl_files.write_ctrl_fitting(
                ctrl, fitting_method, fitting_atom)
        ctrl.write(b"[OPTION]\n")
        ctrl_files.write_kwargs(
                ctrl,
                check_only=check_only,
                analysis_atom=analysis_atom,
                )

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)
        msgbuf, MSG_LEN = make_msgbuf()
        status = ctypes.c_int(0)

        with suppress_stdout_capture_stderr() as captured:
            LibGenesis().lib.ra_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    ctrl_bytes,
                    ctypes.c_int(ctrl_len),
                    ctypes.byref(result_rmsd_c),
                    ctypes.byref(status),
                    msgbuf,
                    ctypes.c_int(MSG_LEN),
                    )

        if status.value != 0:
            error_msg = msgbuf.value.decode("utf-8", "replace")
            raise_fortran_error(
                error_msg,
                code=status.value,
                stderr_output=captured.stderr
            )
        n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
        result_rmsd = (c2py_util.conv_double_ndarray(
            result_rmsd_c, n_frame_c.value)
                    if result_rmsd_c else None)
        return RmsdAnalysisResult(
                result_rmsd)
    finally:
        if result_rmsd_c:
            LibGenesis().lib.deallocate_double(
                    ctypes.byref(result_rmsd_c),
                    ctypes.byref(n_frame_c))
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


def rmsd_analysis_zerocopy(
        molecule: SMolecule,
        trajs: STrajectories,
        analysis_selection: str,
        ana_period: int = 1,
        mass_weighted: bool = False,
        ) -> RmsdAnalysisResult:
    """
    Executes rmsd_analysis with true zero-copy interface.

    This function passes the mass and reference coordinate arrays directly from
    Python NumPy to Fortran, eliminating all data copies. The Fortran code
    creates views of the Python arrays using C_F_POINTER.

    NOTE: This version does NOT perform structural alignment (fitting).
    Use this when:
    - Coordinates are already aligned (e.g., trajectory was centered/aligned)
    - Fitting is done in Python (e.g., using MDAnalysis or MDTraj)
    - You want to measure raw coordinate deviation without superposition

    Args:
        molecule: Molecular structure (provides mass and reference coordinates)
        trajs: Trajectories to analyze
        analysis_selection: GENESIS selection string (e.g., "an:CA", "heavy")
        ana_period: Analysis period (default: 1)
        mass_weighted: Use mass weighting for RMSD calculation (default: False)

    Returns:
        RmsdAnalysisResult containing the RMSD array

    Example:
        >>> # Calculate RMSD for CA atoms (assumes pre-aligned trajectory)
        >>> result = rmsd_analysis_zerocopy(mol, trajs, "an:CA")
        >>> print(result.rmsd)
    """
    lib = LibGenesis().lib

    # Get atom indices using GENESIS selection
    analysis_indices = selection(molecule, analysis_selection)
    n_analysis = len(analysis_indices)

    # Ensure arrays are contiguous and correct dtype
    mass = np.ascontiguousarray(molecule.mass, dtype=np.float64)
    # Reference coordinates: Fortran expects (3, n_atoms)
    ref_coord = np.asfortranarray(molecule.atom_coord.T, dtype=np.float64)

    # Get pointers (zero-copy)
    mass_ptr = mass.ctypes.data_as(ctypes.c_void_p)
    ref_coord_ptr = ref_coord.ctypes.data_as(ctypes.c_void_p)

    # Output variables
    result_rmsd_c = ctypes.c_void_p()
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    try:
        with suppress_stdout_capture_stderr() as captured:
            lib.rmsd_analysis_zerocopy_c(
                mass_ptr,
                ref_coord_ptr,
                ctypes.c_int(molecule.num_atoms),
                ctypes.byref(trajs.get_c_obj()),
                ctypes.c_int(ana_period),
                analysis_indices.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(n_analysis),
                ctypes.c_int(1 if mass_weighted else 0),
                ctypes.byref(result_rmsd_c),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen),
            )

        # Check for errors
        if status.value != 0:
            error_msg = msg.value.decode('utf-8', errors='replace').strip()
            stderr_output = captured.stderr if captured else ""
            raise_fortran_error(status.value, error_msg, stderr_output)

        # Convert result to numpy array
        n_frame = int(trajs.nframe / ana_period)
        result_rmsd = (c2py_util.conv_double_ndarray(result_rmsd_c, n_frame)
                       if result_rmsd_c else None)

        return RmsdAnalysisResult(result_rmsd)

    finally:
        # Cleanup - only deallocate results, not the input arrays (Python owns them)
        lib.deallocate_rmsd_results_c()
        # Note: NO deallocate_s_molecule_c call - we didn't allocate any s_molecule_c!


DrmsAnalysisResult = namedtuple(
        'DrmsAnalysisResult',
        ['drms'])


def drms_analysis(
        molecule: SMolecule, trajs: STrajectories,
        ana_period: Optional[int] = 1,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        fitting_method: Optional[str] = None,
        fitting_atom: Optional[int] = None,
        check_only: Optional[bool] = None,
        contact_groups: Optional[int] = None,
        ignore_hydrogen: Optional[bool] = None,
        two_states: Optional[bool] = None,
        avoid_bonding: Optional[bool] = None,
        exclude_residues: Optional[int] = None,
        minimum_distance: Optional[float] = None,
        maximum_distance: Optional[float] = None,
        pbc_correct: Optional[bool] = None,
        verbose: Optional[bool] = None,
        ) -> DrmsAnalysisResult:
    """
    Executes drms_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
    Returns:
        drms
    """
    if ana_period is None:
        ana_period = 1
    mol_c = None
    ana_period_c = ctypes.c_int(ana_period)
    result_drms_c = ctypes.c_void_p()
    try:
        mol_c = molecule.to_SMoleculeC()
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_output(
                ctrl,
                rmsfile="dummy.rms")
        ctrl_files.write_ctrl_selection(
                ctrl, selection_group, selection_mole_name)
        ctrl_files.write_ctrl_fitting(
                ctrl, fitting_method, fitting_atom)
        ctrl.write(b"[OPTION]\n")
        ctrl_files.write_kwargs(
                ctrl,
                check_only=check_only,
                contact_groups=contact_groups,
                ignore_hydrogen=ignore_hydrogen,
                two_states=two_states,
                avoid_bonding=avoid_bonding,
                exclude_residues=exclude_residues,
                minimum_distance=minimum_distance,
                maximum_distance=maximum_distance,
                pbc_correct=pbc_correct,
                verbose=verbose,
                )

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)
        msgbuf, MSG_LEN = make_msgbuf()
        status = ctypes.c_int(0)

        with suppress_stdout_capture_stderr() as captured:
            LibGenesis().lib.dr_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    ctrl_bytes,
                    ctypes.c_int(ctrl_len),
                    ctypes.byref(result_drms_c),
                    ctypes.byref(status),
                    msgbuf,
                    ctypes.c_int(MSG_LEN),
                    )
        if status.value != 0:
            error_msg = msgbuf.value.decode("utf-8", "replace")
            raise_fortran_error(
                error_msg,
                code=status.value,
                stderr_output=captured.stderr
            )

        n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
        result_drms = (c2py_util.conv_double_ndarray(
            result_drms_c, n_frame_c.value)
                    if result_drms_c else None)
        return DrmsAnalysisResult(
                result_drms)
    finally:
        if result_drms_c:
            LibGenesis().lib.deallocate_double(
                    ctypes.byref(result_drms_c),
                    ctypes.byref(n_frame_c))
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


def drms_analysis_zerocopy(
        trajs: STrajectories,
        contact_list: np.ndarray,
        contact_dist: np.ndarray,
        ana_period: int = 1,
        pbc_correct: bool = False,
        ) -> DrmsAnalysisResult:
    """
    Executes drms_analysis with true zero-copy interface.

    This function passes the contact list and distance arrays directly from
    Python NumPy to Fortran, eliminating all data copies. The Fortran code
    creates views of the Python arrays using C_F_POINTER.

    Args:
        trajs: Trajectories to analyze
        contact_list: Contact atom pairs as (2, n_contact) array with 1-indexed
                      atom indices. Each column is [atom1_idx, atom2_idx].
        contact_dist: Reference distances for each contact pair (n_contact,)
        ana_period: Analysis period (default: 1)
        pbc_correct: Apply PBC correction for distances (default: False)

    Returns:
        DrmsAnalysisResult containing the DRMS array

    Example:
        >>> # Create contact list from atom selections
        >>> contact_list = np.array([[1, 2, 3], [10, 11, 12]], dtype=np.int32)
        >>> contact_dist = np.array([5.0, 6.0, 7.0], dtype=np.float64)
        >>> result = drms_analysis_zerocopy(trajs, contact_list, contact_dist)
        >>> print(result.drms)
    """
    lib = LibGenesis().lib

    # Ensure arrays are contiguous and correct dtype
    # contact_list should be (2, n_contact) with Fortran order for column-major
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

    # Get pointers (zero-copy)
    contact_list_ptr = contact_list_f.ctypes.data_as(ctypes.c_void_p)
    contact_dist_ptr = contact_dist_f.ctypes.data_as(ctypes.c_void_p)

    # Output variables
    result_drms_c = ctypes.c_void_p()
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    try:
        with suppress_stdout_capture_stderr() as captured:
            lib.drms_analysis_zerocopy_c(
                contact_list_ptr,
                contact_dist_ptr,
                ctypes.c_int(n_contact),
                ctypes.byref(trajs.get_c_obj()),
                ctypes.c_int(ana_period),
                ctypes.c_int(1 if pbc_correct else 0),
                ctypes.byref(result_drms_c),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen),
            )

        # Check for errors
        if status.value != 0:
            error_msg = msg.value.decode('utf-8', errors='replace').strip()
            stderr_output = captured.stderr if captured else ""
            raise_fortran_error(status.value, error_msg, stderr_output)

        # Convert result to numpy array
        n_frame = int(trajs.nframe / ana_period)
        result_drms = (c2py_util.conv_double_ndarray(result_drms_c, n_frame)
                       if result_drms_c else None)

        return DrmsAnalysisResult(result_drms)

    finally:
        # Cleanup - only deallocate results, not the input arrays (Python owns them)
        lib.deallocate_drms_results_c()
        # Note: NO deallocate_s_molecule_c call - we didn't allocate any s_molecule_c!


MsdAnalysisResult = namedtuple(
        'MsdAnalysisResult',
        ['msd'])


def msd_analysis(
        molecule: SMolecule, trajs: STrajectories,
        ana_period: Optional[int] = 1,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        selection: Optional[Iterable[str]] = None,
        mode: Optional[Iterable[int]] = None,
        check_only: Optional[bool] = None,
        oversample: Optional[bool] = None,
        delta: Optional[int] = None,
        ) -> MsdAnalysisResult:
    """
    Executes msd_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
    Returns:
        msd
    """
    mol_c = None
    ana_period_c = ctypes.c_int(ana_period)
    num_analysis_mols_c = ctypes.c_int(0)
    num_delta_c = ctypes.c_int(0)
    result_msd_c = ctypes.c_void_p()
    try:
        mol_c = molecule.to_SMoleculeC()
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_output(
                ctrl,
                msdfile="dummy.msd")
        ctrl_files.write_ctrl_selection(
                ctrl, selection_group, selection_mole_name)
        ctrl_files.write_ctrl_molecule_selection(
                ctrl, selection, mode)
        ctrl.write(b"[OPTION]\n")
        ctrl_files.write_kwargs(
                ctrl,
                check_only=check_only,
                oversample=oversample,
                delta=delta,
                )

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)
        with suppress_stdout_capture_stderr() as captured:
            LibGenesis().lib.ma_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    ctrl_bytes,
                    ctypes.c_int(ctrl_len),
                    ctypes.byref(result_msd_c),
                    ctypes.byref(num_analysis_mols_c),
                    ctypes.byref(num_delta_c),
                    )
        result_msd = c2py_util.conv_double_ndarray(
            result_msd_c, [num_delta_c.value, num_analysis_mols_c.value])
        return MsdAnalysisResult(
                result_msd)
    finally:
        if result_msd_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_msd_c),
                    ctypes.byref(num_delta_c),
                    ctypes.byref(num_analysis_mols_c))
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


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
    mol_c = None
    ana_period_c = ctypes.c_int(ana_period)

    try:
        mol_c = molecule.to_SMoleculeC()
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

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)
        msgbuf, MSG_LEN = make_msgbuf()
        status = ctypes.c_int(0)

        with suppress_stdout_capture_stderr() as captured:
            LibGenesis().lib.hb_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    ctrl_bytes,
                    ctypes.c_int(ctrl_len),
                    ctypes.byref(result),
                    ctypes.byref(status),
                    msgbuf,
                    ctypes.c_int(MSG_LEN),
                    )

        if status.value != 0:
            error_msg = msgbuf.value.decode("utf-8", "replace")
            raise_fortran_error(
                error_msg,
                code=status.value,
                stderr_output=captured.stderr
            )

        s = c2py_util.conv_string(result)
        return s

    finally:
        if result:
            try:
                LibGenesis().lib.deallocate_c_string(ctypes.byref(result))
            except Exception:
                pass


def diffusion_analysis(msd_data: npt.NDArray[np.float64],
                       time_step: Optional[int] = None,
                       start: Optional[str] = None,
                       ) -> npt.NDArray[np.float64]:
    """
    Executes diffusion_analysis.

    Args:
    Returns:
        diffusion
    """
    result = ctypes.c_void_p(None)
    mol_c = None

    c_msd = None
    c_out = ctypes.c_void_p(0)
    if msd_data.ndim != 2:
        raise GenesisValidationError(f"msd_data must be 2D, got {msd_data.ndim}D")
    d0 = ctypes.c_int(msd_data.shape[0])
    d1 = ctypes.c_int(msd_data.shape[1])

    try:
        c_msd = ctypes.c_void_p(
                LibGenesis().lib.allocate_c_double_array2(
                    ctypes.byref(d0),
                    ctypes.byref(d1)))
        py2c_util.write_double_ndarray(msd_data, c_msd)
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_output(
                ctrl,
                outfile="dummy.out")
        ctrl.write(b"[OPTION]\n")
        ctrl_files.write_kwargs(
                ctrl,
                time_step=time_step,
                start=start,
                )

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)
        msgbuf, MSG_LEN = make_msgbuf()
        status = ctypes.c_int(0)
        with suppress_stdout_capture_stderr() as captured:
            LibGenesis().lib.diffusion_analysis_c(
                ctypes.byref(c_msd),
                ctypes.byref(d0),
                ctypes.byref(d1),
                ctrl_bytes,
                ctypes.c_int(ctrl_len),
                ctypes.byref(c_out),
                ctypes.byref(status),
                msgbuf,
                ctypes.c_int(MSG_LEN),
                )
        if status.value != 0:
            error_msg = msgbuf.value.decode("utf-8", "replace")
            raise_fortran_error(
                error_msg,
                code=status.value,
                stderr_output=captured.stderr
            )

        return c2py_util.conv_double_ndarray(
                c_out, (msd_data.shape[0], msd_data.shape[1] * 2 - 1))
    finally:
        if c_msd:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(c_msd), ctypes.byref(d0), ctypes.byref(d1))
        if c_out:
            ci = ctypes.c_int(msd_data.shape[0] * (msd_data.shape[1] * 2 - 1))
            LibGenesis().lib.deallocate_double(
                    ctypes.byref(c_out), ctypes.byref(ci))


DiffusionAnalysisResult = namedtuple(
        'DiffusionAnalysisResult',
        ['out_data', 'diffusion_coefficients'])


def diffusion_analysis_zerocopy(
        msd_data: npt.NDArray[np.float64],
        time_step: float = 1.0,
        distance_unit: float = 1.0,
        ndofs: int = 3,
        start_step: int = 1,
        stop_step: Optional[int] = None,
        ) -> DiffusionAnalysisResult:
    """
    Executes diffusion_analysis with true zero-copy interface.

    This function passes the MSD data array directly from Python NumPy to
    Fortran, eliminating unnecessary data copies. Parameters are passed directly
    instead of through a control file.

    Args:
        msd_data: 2D numpy array of shape (ndata, ncols) where:
                  - column 0: time steps (integers or floats)
                  - columns 1+: MSD values for each set
        time_step: Time per step in ps (default: 1.0)
        distance_unit: Distance unit factor (default: 1.0 for Angstrom)
        ndofs: Degrees of freedom (default: 3 for 3D diffusion)
        start_step: Start step for linear fitting (1-indexed, default: 1)
        stop_step: Stop step for linear fitting (1-indexed, default: ndata)

    Returns:
        DiffusionAnalysisResult containing:
        - out_data: 2D array with time, MSD, and fitted values
        - diffusion_coefficients: 1D array of diffusion coefficients (cm^2/s)

    Example:
        >>> # MSD data with time column and one MSD set
        >>> msd = np.array([[0, 0.0], [1, 0.1], [2, 0.3], ...])
        >>> result = diffusion_analysis_zerocopy(msd, time_step=0.01)
        >>> print(f"D = {result.diffusion_coefficients[0]:.2e} cm^2/s")
    """
    lib = LibGenesis().lib

    # Validate input
    if msd_data.ndim != 2:
        raise GenesisValidationError(f"msd_data must be 2D, got {msd_data.ndim}D")
    if msd_data.shape[1] < 2:
        raise GenesisValidationError("msd_data must have at least 2 columns")

    ndata = msd_data.shape[0]
    ncols = msd_data.shape[1]
    n_sets = ncols - 1

    # Set default stop_step
    if stop_step is None:
        stop_step = ndata

    # Ensure array is Fortran order (column-major)
    msd_f = np.asfortranarray(msd_data.T, dtype=np.float64)

    # Get pointer (zero-copy)
    msd_ptr = msd_f.ctypes.data_as(ctypes.c_void_p)

    # Output variables
    out_data_ptr = ctypes.c_void_p()
    diff_coeff_ptr = ctypes.c_void_p()
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    try:
        with suppress_stdout_capture_stderr() as captured:
            lib.diffusion_analysis_zerocopy_c(
                msd_ptr,
                ctypes.c_int(ndata),
                ctypes.c_int(ncols),
                ctypes.c_double(time_step),
                ctypes.c_double(distance_unit),
                ctypes.c_int(ndofs),
                ctypes.c_int(start_step),
                ctypes.c_int(stop_step),
                ctypes.byref(out_data_ptr),
                ctypes.byref(diff_coeff_ptr),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen),
            )

        # Check for errors
        if status.value != 0:
            error_msg = msg.value.decode('utf-8', errors='replace').strip()
            stderr_output = captured.stderr if captured else ""
            raise_fortran_error(status.value, error_msg, stderr_output)

        # Convert results to numpy arrays
        out_ncols = 2 * n_sets + 1  # time + (msd + fit) * n_sets

        # Read flat array and reshape with Fortran order
        out_data_flat = c2py_util.conv_double_ndarray(
            out_data_ptr, out_ncols * ndata) if out_data_ptr else None
        diff_coeff = c2py_util.conv_double_ndarray(
            diff_coeff_ptr, n_sets) if diff_coeff_ptr else None

        # Reshape with Fortran order (column-major) and transpose to (ndata, ncols)
        if out_data_flat is not None:
            out_data = out_data_flat.reshape((out_ncols, ndata), order='F').T
        else:
            out_data = None

        return DiffusionAnalysisResult(out_data, diff_coeff)

    finally:
        # Cleanup - deallocate Fortran allocated results
        lib.deallocate_diffusion_results_c()
        # Note: msd_f is a Python array, no need to deallocate


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
    mol_c = None
    pdb_ave_c = ctypes.c_void_p(None)
    ana_period_c = ctypes.c_int(ana_period)

    try:
        mol_c = molecule.to_SMoleculeC()
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

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)
        msgbuf, MSG_LEN = make_msgbuf()
        status = ctypes.c_int(0)

        with suppress_stdout_capture_stderr() as captured:
           LibGenesis().lib.aa_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    ctrl_bytes,
                    ctypes.c_int(ctrl_len),
                    ctypes.byref(pdb_ave_c),
                    ctypes.byref(status),
                    msgbuf,
                    ctypes.c_int(MSG_LEN),
                    )
        if status.value != 0:
            error_msg = msgbuf.value.decode("utf-8", "replace")
            raise_fortran_error(
                error_msg,
                code=status.value,
                stderr_output=captured.stderr
            )

        if pdb_ave_c:
            pdb_ave = c2py_util.conv_string(pdb_ave_c)
            LibGenesis().lib.deallocate_c_string(ctypes.byref(pdb_ave_c))
        else:
            pdb_ave = None
        return AvecrdAnalysisResult(pdb_ave)
    finally:
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


def wham_analysis(
        psffile: Optional[str] = None,
        prmtopfile: Optional[str] = None,
        ambcrdfile: Optional[str] = None,
        grotopfile: Optional[str] = None,
        grocrdfile: Optional[str] = None,
        pdbfile: Optional[str] = None,
        dcdfile: Optional[str] = None,
        cvfile: Optional[str] = None,
        check_only: Optional[bool] = None,
        allow_backup: Optional[bool] = None,
        dimension: Optional[int] = None,
        nblocks: Optional[int] = None,
        temperature: Optional[float] = None,
        tolerance: Optional[float] = None,
        rest_function: Optional[Iterable[int]] = None,
        grids: Optional[Iterable[tuple[float, float, int]]] = None,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        function: Optional[Iterable[str]] = None,
        select_index: Iterable[Iterable[int]] = None,
        constant: Iterable[Iterable[float]] = None,
        reference: Iterable[Iterable[float]] = None,
        is_periodic: Iterable[bool] = None,
        box_size: Iterable[float] = None,
        ):
    """
    Executes wham_analysis.

    Args:
        ctrl_path:

    Returns:
        pmf
    """
    # === INPUT VALIDATION ===
    from .file_validators import validate_file_exists, validate_file_pattern
    from .param_validators import validate_positive, validate_range

    # Validate input files
    validate_file_exists(psffile, "psffile", required=False)
    validate_file_exists(prmtopfile, "prmtopfile", required=False)
    validate_file_exists(grotopfile, "grotopfile", required=False)
    validate_file_pattern(dcdfile, "dcdfile", required=False)
    validate_file_pattern(cvfile, "cvfile", required=False)

    # Validate parameters
    if dimension is not None:
        validate_range(dimension, 1, 2, "dimension")
    if temperature is not None:
        validate_positive(temperature, "temperature")
    if nblocks is not None:
        validate_positive(nblocks, "nblocks")
    # === END VALIDATION ===

    result_pmf_c = ctypes.c_void_p(None)
    n_bins = ctypes.c_int(0)
    n_bin_x = ctypes.c_int(0)
    try:
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_input(
                ctrl,
                psffile=psffile,
                prmtopfile=prmtopfile,
                ambcrdfile=ambcrdfile,
                grotopfile=grotopfile,
                grocrdfile=grocrdfile,
                pdbfile=pdbfile,
                dcdfile=dcdfile,
                cvfile=cvfile,
                )
        ctrl_files.write_ctrl_output(
                ctrl,
                pmffile="dummy.pmf")
        ctrl.write(b'[WHAM]\n')
        ctrl_files.write_kwargs(
                ctrl,
                check_only=check_only,
                allow_backup=allow_backup,
                dimension=dimension,
                nblocks=nblocks,
                temperature=temperature,
                tolerance=tolerance,
                rest_function=ctrl_files.NumberingData(rest_function),
                grids=ctrl_files.NumberingData(grids),
                )
        ctrl_files.write_ctrl_selection(
                ctrl, selection_group, selection_mole_name)
        ctrl.write(b'[RESTRAINTS]\n')
        ctrl_files.write_kwargs(
                ctrl,
                function=ctrl_files.NumberingData(function),
                select_index=ctrl_files.NumberingData(select_index),
                constant=ctrl_files.NumberingData(constant),
                reference=ctrl_files.NumberingData(reference),
                is_periodic=ctrl_files.NumberingData(is_periodic),
                box_size=ctrl_files.NumberingData(box_size),
                )

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)
        msgbuf, MSG_LEN = make_msgbuf()
        status = ctypes.c_int(0)
        with suppress_stdout_capture_stderr() as captured:
            LibGenesis().lib.wa_analysis_c(
                    ctrl_bytes,
                    ctypes.c_int(ctrl_len),
                    ctypes.byref(result_pmf_c),
                    ctypes.byref(n_bins),
                    ctypes.byref(n_bin_x),
                    ctypes.byref(status),
                    msgbuf,
                    ctypes.c_int(MSG_LEN),
                    )
        if status.value != 0:
            error_msg = msgbuf.value.decode("utf-8", "replace")
            raise_fortran_error(
                error_msg,
                code=status.value,
                stderr_output=captured.stderr
            )

        result_pmf = c2py_util.conv_double_ndarray(
                result_pmf_c, [n_bins.value, n_bin_x.value])

        return result_pmf
    finally:
        if result_pmf_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_pmf_c),
                    ctypes.byref(n_bins), ctypes.byref(n_bin_x))


def mbar_analysis(
        psffile: Optional[str] = None,
        prmtopfile: Optional[str] = None,
        ambcrdfile: Optional[str] = None,
        grotopfile: Optional[str] = None,
        grocrdfile: Optional[str] = None,
        pdbfile: Optional[str] = None,
        dcdfile: Optional[str] = None,
        cvfile: Optional[str] = None,
        check_only: Optional[bool] = None,
        allow_backup: Optional[bool] = None,
        nreplica: Optional[int] = None,
        input_type: Optional[str] = None,
        dimension: Optional[int] = None,
        temperature: Optional[float] = None,
        target_temperature: Optional[float] = None,
        nblocks: Optional[int] = None,
        tolerance: Optional[float] = None,
        rest_function: Optional[Iterable[int]] = None,
        grids: Optional[Iterable[tuple[float, float, int]]] = None,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        constant: Iterable[Iterable[float]] = None,
        reference: Iterable[Iterable[float]] = None,
        is_periodic: Iterable[bool] = None,
        box_size: Iterable[float] = None,
        ):
    """
    Executes mbar_analysis.

    Args:
    Returns:
        fene
    """
    # === INPUT VALIDATION ===
    from .file_validators import validate_file_exists, validate_file_pattern
    from .param_validators import validate_positive, validate_range

    # Validate input files
    validate_file_exists(psffile, "psffile", required=False)
    validate_file_exists(prmtopfile, "prmtopfile", required=False)
    validate_file_exists(grotopfile, "grotopfile", required=False)
    validate_file_pattern(dcdfile, "dcdfile", required=False)
    validate_file_pattern(cvfile, "cvfile", required=False)

    # Validate parameters
    if dimension is not None:
        validate_range(dimension, 1, 2, "dimension")
    if temperature is not None:
        validate_positive(temperature, "temperature")
    if nreplica is not None:
        validate_positive(nreplica, "nreplica")
    if nblocks is not None:
        validate_positive(nblocks, "nblocks")
    # === END VALIDATION ===

    result_fene_c = ctypes.c_void_p(None)
    n_replica = ctypes.c_int(0)
    n_blocks = ctypes.c_int(0)
    try:
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_input(
                ctrl,
                psffile=psffile,
                prmtopfile=prmtopfile,
                ambcrdfile=ambcrdfile,
                grotopfile=grotopfile,
                grocrdfile=grocrdfile,
                pdbfile=pdbfile,
                dcdfile=dcdfile,
                cvfile=cvfile,
                )
        ctrl_files.write_ctrl_output(
                ctrl,
                fenefile="fene.dat")
        ctrl.write(b'[MBAR]\n')
        ctrl_files.write_kwargs(
                ctrl,
                check_only=check_only,
                allow_backup=allow_backup,
                nreplica=nreplica,
                input_type=input_type,
                dimension=dimension,
                temperature=temperature,
                target_temperature=target_temperature,
                nblocks=nblocks,
                tolerance=tolerance,
                rest_function=ctrl_files.NumberingData(rest_function),
                grids=ctrl_files.NumberingData(grids),
                )
        ctrl_files.write_ctrl_selection(
                ctrl, selection_group, selection_mole_name)
        ctrl.write(b'[RESTRAINTS]\n')
        ctrl_files.write_kwargs(
                ctrl,
                constant=ctrl_files.NumberingData(constant),
                reference=ctrl_files.NumberingData(reference),
                is_periodic=ctrl_files.NumberingData(is_periodic),
                box_size=ctrl_files.NumberingData(box_size),
                )

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)
        msgbuf, MSG_LEN = make_msgbuf()
        status = ctypes.c_int(0)

        with suppress_stdout_capture_stderr() as captured:
            LibGenesis().lib.mbar_analysis_c(
                    ctrl_bytes,
                    ctypes.c_int(ctrl_len),
                    ctypes.byref(result_fene_c),
                    ctypes.byref(n_replica),
                    ctypes.byref(n_blocks),
                    ctypes.byref(status),
                    msgbuf,
                    ctypes.c_int(MSG_LEN),
                    )
        if status.value != 0:
            error_msg = msgbuf.value.decode("utf-8", "replace")
            raise_fortran_error(
                error_msg,
                code=status.value,
                stderr_output=captured.stderr
            )

        result_fene = c2py_util.conv_double_ndarray(
                result_fene_c, [n_replica.value, n_blocks.value])
        return result_fene
    finally:
        if result_fene_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_fene_c),
                    ctypes.byref(n_replica), ctypes.byref(n_blocks))


class KmeansClusteringResult(NamedTuple):
    mols_from_pdb: SMolecule
    cluster_idxs: npt.NDArray[np.int64]


def kmeans_clustering(
        molecule: SMolecule, trajs: STrajectories,
        ana_period: Optional[int] = 1,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        fitting_method: Optional[str] = None,
        fitting_atom: Optional[int] = None,
        zrot_ngrid: Optional[int] = None,
        zrot_grid_size: Optional[float] = None,
        mass_weight: Optional[bool] = None,
        check_only: Optional[bool] = None,
        allow_backup: Optional[bool] = None,
        analysis_atom: Optional[int] = None,
        num_clusters: Optional[int] = None,
        max_iteration: Optional[int] = None,
        stop_threshold: Optional[float] = None,
        num_iterations: Optional[int] = None,
        trjout_atom: Optional[int] = None,
        trjout_format: Optional[str] = None,
        trjout_type: Optional[str] = None,
        iseed: Optional[int] = None,
        ) -> KmeansClusteringResult:
    """
    Executes kmeans_clustering.

    Args:
        molecule:
        trajs:
        ana_period:
    Returns:
        (pdb string, cluster indices)
    """
    mol_c = None
    pdb_c = ctypes.c_void_p()
    cluster_idxs_c = ctypes.c_void_p()
    ana_period_c = ctypes.c_int(ana_period)
    cluster_size = ctypes.c_int(0)
    try:
        mol_c = molecule.to_SMoleculeC()
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_output(
                ctrl,
                indexfile="dummy.idx",
                pdbfile="dummy_{}.pdb",
                trjfile="dummy{}.trj")
        ctrl_files.write_ctrl_selection(
                ctrl, selection_group, selection_mole_name)
        ctrl_files.write_ctrl_fitting(
                ctrl, fitting_method, fitting_atom,
                zrot_ngrid, zrot_grid_size, mass_weight)
        ctrl.write(b"[OPTION]\n")
        ctrl_files.write_kwargs(
                ctrl,
                check_only=check_only,
                allow_backup=allow_backup,
                analysis_atom=analysis_atom,
                num_clusters=num_clusters,
                max_iteration=max_iteration,
                stop_threshold=stop_threshold,
                num_iterations=num_iterations,
                trjout_atom=trjout_atom,
                trjout_format=trjout_format,
                trjout_type=trjout_type,
                iseed=iseed,
                )

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)
        msgbuf, MSG_LEN = make_msgbuf()
        status = ctypes.c_int(0)
        with suppress_stdout_capture_stderr() as captured:
            LibGenesis().lib.kc_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    ctrl_bytes,
                    ctypes.c_int(ctrl_len),
                    ctypes.byref(pdb_c),
                    ctypes.byref(cluster_idxs_c),
                    ctypes.byref(cluster_size),
                    ctypes.byref(status),
                    msgbuf,
                    ctypes.c_int(MSG_LEN),
                    )

        if status.value != 0:
            error_msg = msgbuf.value.decode("utf-8", "replace")
            raise_fortran_error(
                error_msg,
                code=status.value,
                stderr_output=captured.stderr
            )

        if pdb_c:
            pdb_str = c2py_util.conv_string(pdb_c)
            LibGenesis().lib.deallocate_c_string(ctypes.byref(pdb_c))
            pdb_mols = []
            for pdb_block in extract_model_blocks(pdb_str):
                with tempfile.NamedTemporaryFile(
                        dir=os.getcwd(), delete=True) as pdb_file:
                    pdb_file.write(pdb_block.encode())
                    pdb_file.seek(0)
                    pdb_mols.append(SMolecule.from_file(pdb=pdb_file.name))
        else:
            pdb_mols = None
        cluster_idxs = (c2py_util.conv_int_ndarray(
                cluster_idxs_c, cluster_size.value)
                        if cluster_idxs_c else None)
        return KmeansClusteringResult(pdb_mols, cluster_idxs)
    finally:
        if cluster_idxs_c:
            LibGenesis().lib.deallocate_int(
                    ctypes.byref(cluster_idxs_c), ctypes.byref(cluster_size))
        if pdb_c:
            LibGenesis().lib.deallocate_c_string(ctypes.byref(pdb_c))
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


def extract_model_blocks(pdb_string):
    """
    与えられたPDB形式の文字列から、MODEL行からENDMDL行の間のスライスを順に返す。
    文字列のコピーは行わない。

    Parameters:
        pdb_string (str): PDB形式の文字列。

    Yields:
        str: 各MODELブロック（元の文字列のスライス）。
    """
    start = None  # MODEL行の開始位置
    end = None    # ENDMDL行の終了位置

    i = 0
    while i < len(pdb_string):
        # MODEL行の開始位置を探す
        if pdb_string.startswith("MODEL", i):
            start = i
            # MODEL行の終わりまで進める
            while i < len(pdb_string) and pdb_string[i] != '\n':
                i += 1
            i += 1  # 改行をスキップ

        # ENDMDL行の終了位置を探す
        elif pdb_string.startswith("ENDMDL", i):
            end = i
            while i < len(pdb_string) and pdb_string[i] != '\n':
                i += 1
            i += 1  # 改行をスキップ
            yield pdb_string[start:i]  # MODEL行からENDMDL行までのスライスを返す
            start = None
            end = None

        else:
            i += 1  # 次の文字へ進む

@lru_cache(maxsize=1)
def get_msg_len() -> int:
    lib = LibGenesis().lib
    fn = getattr(lib, "get_default_msg_len", None) or getattr(lib, "genesis_get_default_msg_len", None)
    if fn is None:
        return _DEFAULT_MSG_LEN
    fn.restype = ctypes.c_int
    try:
        n = int(fn())
    except Exception:
        return _DEFAULT_MSG_LEN

    return n if n > 0 else _DEFAULT_MSG_LEN

def make_msgbuf():
    n = get_msg_len()
    buf = (ctypes.c_char * n)()
    return buf, n


# ============================================================================
# ATDYN MD/Minimization Functions
# ============================================================================

class AtdynMDResult(NamedTuple):
    """Result from atdyn MD simulation."""
    energies: npt.NDArray[np.float64]  # Shape: (nterms, nframes)
    final_coords: npt.NDArray[np.float64]  # Shape: (3, natom)
    energy_labels: tuple  # ('total', 'bond', 'angle', ...)


class AtdynMinResult(NamedTuple):
    """Result from atdyn energy minimization."""
    energies: npt.NDArray[np.float64]  # Shape: (nterms, nsteps)
    final_coords: npt.NDArray[np.float64]  # Shape: (3, natom)
    converged: bool
    final_gradient: float
    energy_labels: tuple  # ('total', 'bond', 'angle', ...)


_ENERGY_LABELS = ('total', 'bond', 'angle', 'urey_bradley', 'dihedral',
                  'improper', 'electrostatic', 'van_der_waals')


def run_atdyn_md(
        # Input files
        psffile: Optional[str] = None,
        prmtopfile: Optional[str] = None,
        ambcrdfile: Optional[str] = None,
        pdbfile: Optional[str] = None,
        rstfile: Optional[str] = None,
        topfile: Optional[str] = None,
        parfile: Optional[str] = None,
        strfile: Optional[str] = None,
        grotopfile: Optional[str] = None,
        grocrdfile: Optional[str] = None,
        # Energy parameters
        forcefield: Optional[str] = None,
        electrostatic: Optional[str] = None,
        switchdist: Optional[float] = None,
        cutoffdist: Optional[float] = None,
        pairlistdist: Optional[float] = None,
        vdw_force_switch: Optional[bool] = None,
        implicit_solvent: Optional[str] = None,
        output_style: Optional[str] = None,
        pme_alpha: Optional[float] = None,
        pme_ngrid_x: Optional[int] = None,
        pme_ngrid_y: Optional[int] = None,
        pme_ngrid_z: Optional[int] = None,
        pme_nspline: Optional[int] = None,
        dispersion_corr: Optional[str] = None,
        # Dynamics parameters
        integrator: Optional[str] = None,
        nsteps: Optional[int] = None,
        timestep: Optional[float] = None,
        eneout_period: Optional[int] = None,
        crdout_period: Optional[int] = None,
        rstout_period: Optional[int] = None,
        nbupdate_period: Optional[int] = None,
        iseed: Optional[int] = None,
        verbose: Optional[bool] = None,
        # Boundary parameters
        boundary_type: Optional[str] = None,
        box_size_x: Optional[float] = None,
        box_size_y: Optional[float] = None,
        box_size_z: Optional[float] = None,
        # Ensemble parameters
        ensemble: Optional[str] = None,
        tpcontrol: Optional[str] = None,
        temperature: Optional[float] = None,
        pressure: Optional[float] = None,
        gamma_t: Optional[float] = None,
        # Constraints parameters
        rigid_bond: Optional[bool] = None,
        shake_iteration: Optional[int] = None,
        shake_tolerance: Optional[float] = None,
        water_model: Optional[str] = None,
        # Output files (optional)
        dcdfile: Optional[str] = None,
        ) -> AtdynMDResult:
    """
    Run atdyn molecular dynamics simulation.

    Args:
        psffile: PSF topology file (CHARMM format)
        prmtopfile: AMBER parameter/topology file
        ambcrdfile: AMBER coordinate file
        pdbfile: PDB coordinate file
        rstfile: GENESIS restart file
        topfile: CHARMM topology RTF file
        parfile: CHARMM parameter file
        strfile: CHARMM stream file
        grotopfile: GROMACS topology file
        grocrdfile: GROMACS coordinate file (.gro)
        forcefield: Force field type (CHARMM, AMBER, GROAMBER, etc.)
        electrostatic: Electrostatic method (PME, CUTOFF, etc.)
        switchdist: Switch distance for vdW (angstrom)
        cutoffdist: Cutoff distance (angstrom)
        pairlistdist: Pairlist distance (angstrom)
        vdw_force_switch: Use force switching for vdW
        implicit_solvent: Implicit solvent model (NONE, GBSA, EEF1, etc.)
        output_style: Output style (GENESIS, CHARMM, etc.)
        pme_alpha: PME alpha parameter
        pme_ngrid_x: PME grid points in X
        pme_ngrid_y: PME grid points in Y
        pme_ngrid_z: PME grid points in Z
        pme_nspline: PME spline order
        dispersion_corr: Dispersion correction (NO, epress, etc.)
        integrator: Integrator (VVER, LEAP, etc.)
        nsteps: Number of MD steps
        timestep: Time step (ps)
        eneout_period: Energy output period
        crdout_period: Coordinate output period
        rstout_period: Restart output period
        nbupdate_period: Nonbond list update period
        iseed: Random seed
        verbose: Verbose output
        boundary_type: Boundary type (PBC, NOBC)
        box_size_x: Box size X (angstrom)
        box_size_y: Box size Y (angstrom)
        box_size_z: Box size Z (angstrom)
        ensemble: Ensemble type (NVE, NVT, NPT, etc.)
        tpcontrol: Thermostat/barostat (LANGEVIN, BUSSI, etc.)
        temperature: Temperature (K)
        pressure: Pressure (atm)
        gamma_t: Langevin friction coefficient
        rigid_bond: Use rigid bonds (SHAKE/SETTLE)
        shake_iteration: Maximum SHAKE iterations
        shake_tolerance: SHAKE convergence tolerance
        water_model: Water model for SETTLE (TIP3, WAT, SOL, etc.)
        dcdfile: Output DCD trajectory file

    Returns:
        AtdynMDResult containing:
            - energies: Array of energy terms (nterms x nframes)
            - final_coords: Final coordinates (3 x natom)
            - energy_labels: Tuple of energy term names

    Raises:
        GenesisValidationError: If input validation fails
        GenesisFortranError: If Fortran code returns an error
    """
    # === INPUT VALIDATION ===
    from .file_validators import (
        validate_file_exists,
        validate_topology_combination,
    )
    from .param_validators import (
        validate_enum,
        validate_positive,
        validate_distance_ordering,
        validate_pme_params,
        validate_shake_params,
        validate_ensemble_params,
        validate_pbc_params,
        FORCEFIELDS,
        ELECTROSTATICS,
        INTEGRATORS,
        BOUNDARY_TYPES,
        IMPLICIT_SOLVENTS,
    )

    # Validate topology - at least one format required
    validate_topology_combination(psffile, prmtopfile, grotopfile)

    # Validate input files exist
    validate_file_exists(psffile, "psffile", required=False)
    validate_file_exists(prmtopfile, "prmtopfile", required=False)
    validate_file_exists(ambcrdfile, "ambcrdfile", required=False)
    validate_file_exists(pdbfile, "pdbfile", required=False)
    validate_file_exists(rstfile, "rstfile", required=False)
    validate_file_exists(topfile, "topfile", required=False)
    validate_file_exists(parfile, "parfile", required=False)
    validate_file_exists(strfile, "strfile", required=False)
    validate_file_exists(grotopfile, "grotopfile", required=False)
    validate_file_exists(grocrdfile, "grocrdfile", required=False)

    # Validate enum parameters
    validate_enum(forcefield, FORCEFIELDS, "forcefield")
    validate_enum(electrostatic, ELECTROSTATICS, "electrostatic")
    validate_enum(integrator, INTEGRATORS, "integrator")
    validate_enum(boundary_type, BOUNDARY_TYPES, "boundary_type")
    validate_enum(implicit_solvent, IMPLICIT_SOLVENTS, "implicit_solvent")

    # Validate numeric parameters
    validate_positive(nsteps, "nsteps", allow_none=False)
    validate_positive(timestep, "timestep")
    validate_positive(eneout_period, "eneout_period")
    validate_positive(crdout_period, "crdout_period")
    validate_positive(rstout_period, "rstout_period")
    validate_positive(nbupdate_period, "nbupdate_period")

    # Validate distance ordering
    validate_distance_ordering(switchdist, cutoffdist, pairlistdist)

    # Validate conditional parameters
    validate_pme_params(
        electrostatic, pme_alpha, pme_ngrid_x, pme_ngrid_y, pme_ngrid_z, pme_nspline
    )
    validate_shake_params(rigid_bond, shake_iteration, shake_tolerance)
    validate_ensemble_params(ensemble, temperature, pressure, tpcontrol)
    validate_pbc_params(boundary_type, box_size_x, box_size_y, box_size_z)
    # === END VALIDATION ===

    result_energies_c = ctypes.c_void_p(None)
    result_nframes_c = ctypes.c_int(0)
    result_nterms_c = ctypes.c_int(0)
    result_final_coords_c = ctypes.c_void_p(None)
    result_natom_c = ctypes.c_int(0)
    status_c = ctypes.c_int(0)
    msgbuf, MSG_LEN = make_msgbuf()

    try:
        # Build control string
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_input(
            ctrl,
            psffile=psffile,
            prmtopfile=prmtopfile,
            ambcrdfile=ambcrdfile,
            pdbfile=pdbfile,
            rstfile=rstfile,
            topfile=topfile,
            parfile=parfile,
            strfile=strfile,
            grotopfile=grotopfile,
            grocrdfile=grocrdfile,
        )
        ctrl_files.write_ctrl_output(
            ctrl,
            dcdfile=dcdfile,
        )
        ctrl_files.write_ctrl_energy(
            ctrl,
            forcefield=forcefield,
            electrostatic=electrostatic,
            switchdist=switchdist,
            cutoffdist=cutoffdist,
            pairlistdist=pairlistdist,
            vdw_force_switch=vdw_force_switch,
            implicit_solvent=implicit_solvent,
            output_style=output_style,
            pme_alpha=pme_alpha,
            pme_ngrid_x=pme_ngrid_x,
            pme_ngrid_y=pme_ngrid_y,
            pme_ngrid_z=pme_ngrid_z,
            pme_nspline=pme_nspline,
            dispersion_corr=dispersion_corr,
        )
        ctrl_files.write_ctrl_dynamics(
            ctrl,
            integrator=integrator,
            nsteps=nsteps,
            timestep=timestep,
            eneout_period=eneout_period,
            crdout_period=crdout_period,
            rstout_period=rstout_period,
            nbupdate_period=nbupdate_period,
            iseed=iseed,
            verbose=verbose,
        )
        ctrl_files.write_ctrl_boundary(
            ctrl,
            type=boundary_type,
            box_size_x=box_size_x,
            box_size_y=box_size_y,
            box_size_z=box_size_z,
        )
        ctrl_files.write_ctrl_ensemble(
            ctrl,
            ensemble=ensemble,
            tpcontrol=tpcontrol,
            temperature=temperature,
            pressure=pressure,
            gamma_t=gamma_t,
        )
        ctrl_files.write_ctrl_constraints(
            ctrl,
            rigid_bond=rigid_bond,
            shake_iteration=shake_iteration,
            shake_tolerance=shake_tolerance,
            water_model=water_model,
        )

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)

        with suppress_stdout_capture_stderr() as captured:
            LibGenesis().lib.atdyn_md_c(
                ctrl_bytes,
                ctypes.c_int(ctrl_len),
                ctypes.byref(result_energies_c),
                ctypes.byref(result_nframes_c),
                ctypes.byref(result_nterms_c),
                ctypes.byref(result_final_coords_c),
                ctypes.byref(result_natom_c),
                ctypes.byref(status_c),
                msgbuf,
                ctypes.c_int(MSG_LEN),
            )

        if status_c.value != 0:
            error_msg = msgbuf.value.decode("utf-8", "replace")
            raise_fortran_error(
                error_msg,
                code=status_c.value,
                stderr_output=captured.stderr
            )

        # Convert results to numpy arrays
        nframes = result_nframes_c.value
        nterms = result_nterms_c.value
        natom = result_natom_c.value

        energies = c2py_util.conv_double_ndarray(
            result_energies_c, [nterms, nframes])
        final_coords = c2py_util.conv_double_ndarray(
            result_final_coords_c, [3, natom])

        return AtdynMDResult(
            energies=energies,
            final_coords=final_coords,
            energy_labels=_ENERGY_LABELS,
        )

    finally:
        # Deallocate results
        LibGenesis().lib.deallocate_atdyn_results_c()


def run_atdyn_min(
        # Input files
        psffile: Optional[str] = None,
        prmtopfile: Optional[str] = None,
        ambcrdfile: Optional[str] = None,
        pdbfile: Optional[str] = None,
        rstfile: Optional[str] = None,
        topfile: Optional[str] = None,
        parfile: Optional[str] = None,
        strfile: Optional[str] = None,
        grotopfile: Optional[str] = None,
        grocrdfile: Optional[str] = None,
        # Energy parameters
        forcefield: Optional[str] = None,
        electrostatic: Optional[str] = None,
        switchdist: Optional[float] = None,
        cutoffdist: Optional[float] = None,
        pairlistdist: Optional[float] = None,
        vdw_force_switch: Optional[bool] = None,
        implicit_solvent: Optional[str] = None,
        output_style: Optional[str] = None,
        pme_alpha: Optional[float] = None,
        pme_ngrid_x: Optional[int] = None,
        pme_ngrid_y: Optional[int] = None,
        pme_ngrid_z: Optional[int] = None,
        pme_nspline: Optional[int] = None,
        dispersion_corr: Optional[str] = None,
        # Minimize parameters
        method: Optional[str] = None,
        nsteps: Optional[int] = None,
        eneout_period: Optional[int] = None,
        crdout_period: Optional[int] = None,
        rstout_period: Optional[int] = None,
        nbupdate_period: Optional[int] = None,
        force_scale_init: Optional[float] = None,
        force_scale_max: Optional[float] = None,
        verbose: Optional[bool] = None,
        tol_rmsg: Optional[float] = None,
        tol_maxg: Optional[float] = None,
        # Constraints parameters
        rigid_bond: Optional[bool] = None,
        # Boundary parameters
        boundary_type: Optional[str] = None,
        box_size_x: Optional[float] = None,
        box_size_y: Optional[float] = None,
        box_size_z: Optional[float] = None,
        # Output files (optional)
        dcdfile: Optional[str] = None,
        ) -> AtdynMinResult:
    """
    Run atdyn energy minimization.

    Args:
        psffile: PSF topology file (CHARMM format)
        prmtopfile: AMBER parameter/topology file
        ambcrdfile: AMBER coordinate file
        pdbfile: PDB coordinate file
        rstfile: GENESIS restart file
        topfile: CHARMM topology RTF file
        parfile: CHARMM parameter file
        strfile: CHARMM stream file
        grotopfile: GROMACS topology file
        grocrdfile: GROMACS coordinate file (.gro)
        forcefield: Force field type (CHARMM, AMBER, GROAMBER, etc.)
        electrostatic: Electrostatic method (PME, CUTOFF, etc.)
        switchdist: Switch distance for vdW (angstrom)
        cutoffdist: Cutoff distance (angstrom)
        pairlistdist: Pairlist distance (angstrom)
        vdw_force_switch: Use force switching for vdW
        implicit_solvent: Implicit solvent model (NONE, GBSA, EEF1, etc.)
        output_style: Output style (GENESIS, CHARMM, etc.)
        pme_alpha: PME alpha parameter
        pme_ngrid_x: PME grid points in X
        pme_ngrid_y: PME grid points in Y
        pme_ngrid_z: PME grid points in Z
        pme_nspline: PME spline order
        dispersion_corr: Dispersion correction (NO, epress, etc.)
        method: Minimization method (SD, LBFGS)
        nsteps: Maximum number of minimization steps
        eneout_period: Energy output period
        crdout_period: Coordinate output period
        rstout_period: Restart output period
        nbupdate_period: Nonbond list update period
        force_scale_init: Initial force scale (SD)
        force_scale_max: Maximum force scale (SD)
        verbose: Verbose output
        tol_rmsg: RMS gradient tolerance for convergence
        tol_maxg: Max gradient tolerance for convergence
        boundary_type: Boundary type (PBC, NOBC)
        box_size_x: Box size X (angstrom)
        box_size_y: Box size Y (angstrom)
        box_size_z: Box size Z (angstrom)
        dcdfile: Output DCD trajectory file

    Returns:
        AtdynMinResult containing:
            - energies: Array of energy terms (nterms x nsteps)
            - final_coords: Final coordinates (3 x natom)
            - converged: Whether minimization converged
            - final_gradient: Final RMS gradient
            - energy_labels: Tuple of energy term names

    Raises:
        GenesisValidationError: If input validation fails
        GenesisFortranError: If Fortran code returns an error
    """
    # === INPUT VALIDATION ===
    from .file_validators import (
        validate_file_exists,
        validate_topology_combination,
    )
    from .param_validators import (
        validate_enum,
        validate_positive,
        validate_distance_ordering,
        validate_pme_params,
        validate_pbc_params,
        FORCEFIELDS,
        ELECTROSTATICS,
        MINIMIZERS,
        BOUNDARY_TYPES,
        IMPLICIT_SOLVENTS,
    )

    # Validate topology - at least one format required
    validate_topology_combination(psffile, prmtopfile, grotopfile)

    # Validate input files exist
    validate_file_exists(psffile, "psffile", required=False)
    validate_file_exists(prmtopfile, "prmtopfile", required=False)
    validate_file_exists(ambcrdfile, "ambcrdfile", required=False)
    validate_file_exists(pdbfile, "pdbfile", required=False)
    validate_file_exists(rstfile, "rstfile", required=False)
    validate_file_exists(topfile, "topfile", required=False)
    validate_file_exists(parfile, "parfile", required=False)
    validate_file_exists(strfile, "strfile", required=False)
    validate_file_exists(grotopfile, "grotopfile", required=False)
    validate_file_exists(grocrdfile, "grocrdfile", required=False)

    # Validate enum parameters
    validate_enum(forcefield, FORCEFIELDS, "forcefield")
    validate_enum(electrostatic, ELECTROSTATICS, "electrostatic")
    validate_enum(method, MINIMIZERS, "method")
    validate_enum(boundary_type, BOUNDARY_TYPES, "boundary_type")
    validate_enum(implicit_solvent, IMPLICIT_SOLVENTS, "implicit_solvent")

    # Validate numeric parameters
    validate_positive(nsteps, "nsteps", allow_none=False)
    validate_positive(eneout_period, "eneout_period")
    validate_positive(crdout_period, "crdout_period")
    validate_positive(rstout_period, "rstout_period")
    validate_positive(nbupdate_period, "nbupdate_period")
    validate_positive(force_scale_init, "force_scale_init")
    validate_positive(force_scale_max, "force_scale_max")
    validate_positive(tol_rmsg, "tol_rmsg")
    validate_positive(tol_maxg, "tol_maxg")

    # Validate distance ordering
    validate_distance_ordering(switchdist, cutoffdist, pairlistdist)

    # Validate conditional parameters
    validate_pme_params(
        electrostatic, pme_alpha, pme_ngrid_x, pme_ngrid_y, pme_ngrid_z, pme_nspline
    )
    validate_pbc_params(boundary_type, box_size_x, box_size_y, box_size_z)
    # === END VALIDATION ===

    result_energies_c = ctypes.c_void_p(None)
    result_nsteps_c = ctypes.c_int(0)
    result_nterms_c = ctypes.c_int(0)
    result_final_coords_c = ctypes.c_void_p(None)
    result_natom_c = ctypes.c_int(0)
    result_converged_c = ctypes.c_int(0)
    result_final_gradient_c = ctypes.c_double(0.0)
    status_c = ctypes.c_int(0)
    msgbuf, MSG_LEN = make_msgbuf()

    try:
        # Build control string
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_input(
            ctrl,
            psffile=psffile,
            prmtopfile=prmtopfile,
            ambcrdfile=ambcrdfile,
            pdbfile=pdbfile,
            rstfile=rstfile,
            topfile=topfile,
            parfile=parfile,
            strfile=strfile,
            grotopfile=grotopfile,
            grocrdfile=grocrdfile,
        )
        ctrl_files.write_ctrl_output(
            ctrl,
            dcdfile=dcdfile,
        )
        ctrl_files.write_ctrl_energy(
            ctrl,
            forcefield=forcefield,
            electrostatic=electrostatic,
            switchdist=switchdist,
            cutoffdist=cutoffdist,
            pairlistdist=pairlistdist,
            vdw_force_switch=vdw_force_switch,
            implicit_solvent=implicit_solvent,
            output_style=output_style,
            pme_alpha=pme_alpha,
            pme_ngrid_x=pme_ngrid_x,
            pme_ngrid_y=pme_ngrid_y,
            pme_ngrid_z=pme_ngrid_z,
            pme_nspline=pme_nspline,
            dispersion_corr=dispersion_corr,
        )
        ctrl_files.write_ctrl_minimize(
            ctrl,
            method=method,
            nsteps=nsteps,
            eneout_period=eneout_period,
            crdout_period=crdout_period,
            rstout_period=rstout_period,
            nbupdate_period=nbupdate_period,
            force_scale_init=force_scale_init,
            force_scale_max=force_scale_max,
            verbose=verbose,
            tol_rmsg=tol_rmsg,
            tol_maxg=tol_maxg,
        )
        ctrl_files.write_ctrl_constraints(
            ctrl,
            rigid_bond=rigid_bond,
        )
        ctrl_files.write_ctrl_boundary(
            ctrl,
            type=boundary_type,
            box_size_x=box_size_x,
            box_size_y=box_size_y,
            box_size_z=box_size_z,
        )

        ctrl_bytes = ctrl.getvalue()
        ctrl_len = len(ctrl_bytes)

        with suppress_stdout_capture_stderr() as captured:
            LibGenesis().lib.atdyn_min_c(
                ctrl_bytes,
                ctypes.c_int(ctrl_len),
                ctypes.byref(result_energies_c),
                ctypes.byref(result_nsteps_c),
                ctypes.byref(result_nterms_c),
                ctypes.byref(result_final_coords_c),
                ctypes.byref(result_natom_c),
                ctypes.byref(result_converged_c),
                ctypes.byref(result_final_gradient_c),
                ctypes.byref(status_c),
                msgbuf,
                ctypes.c_int(MSG_LEN),
            )

        if status_c.value != 0:
            error_msg = msgbuf.value.decode("utf-8", "replace")
            raise_fortran_error(
                error_msg,
                code=status_c.value,
                stderr_output=captured.stderr
            )

        # Convert results to numpy arrays
        nsteps_out = result_nsteps_c.value
        nterms = result_nterms_c.value
        natom = result_natom_c.value

        energies = c2py_util.conv_double_ndarray(
            result_energies_c, [nterms, nsteps_out])
        final_coords = c2py_util.conv_double_ndarray(
            result_final_coords_c, [3, natom])

        return AtdynMinResult(
            energies=energies,
            final_coords=final_coords,
            converged=bool(result_converged_c.value),
            final_gradient=result_final_gradient_c.value,
            energy_labels=_ENERGY_LABELS,
        )

    finally:
        # Deallocate results
        LibGenesis().lib.deallocate_atdyn_results_c()


# =============================================================================
# Subprocess Isolation API
# =============================================================================

def _run_isolated_subprocess(
    func_name: str,
    result_fields: list,
    timeout: Optional[float],
    task_description: str,
    **kwargs
) -> dict:
    """
    Common helper for running atdyn functions in isolated subprocess.

    Args:
        func_name: Name of the genesis_exe function to call
        result_fields: List of field names to extract from result
        timeout: Maximum time in seconds to wait
        task_description: Description for error messages (e.g., "MD simulation")
        **kwargs: Arguments to pass to the function

    Returns:
        Dictionary with the result fields
    """
    import subprocess
    import sys
    import pickle
    import base64

    kwargs_bytes = base64.b64encode(pickle.dumps(kwargs)).decode('ascii')
    result_fields_str = repr(result_fields)

    script = f'''
import sys
import pickle
import base64

try:
    from genepie import genesis_exe
except ImportError:
    sys.path.insert(0, "{os.path.dirname(os.path.dirname(__file__))}")
    from genepie import genesis_exe

kwargs = pickle.loads(base64.b64decode("{kwargs_bytes}"))

try:
    result = genesis_exe.{func_name}(**kwargs)
    output = {{"success": True}}
    for field in {result_fields_str}:
        val = getattr(result, field)
        output[field] = val.tolist() if hasattr(val, 'tolist') else val
except Exception as e:
    import traceback
    output = {{
        "success": False,
        "error": str(e),
        "error_type": type(e).__name__,
        "traceback": traceback.format_exc(),
    }}

sys.stdout.buffer.write(base64.b64encode(pickle.dumps(output)))
'''

    try:
        proc = subprocess.run(
            [sys.executable, '-c', script],
            capture_output=True,
            timeout=timeout,
            env={**os.environ, 'OMP_NUM_THREADS': os.environ.get('OMP_NUM_THREADS', '1')},
        )
    except subprocess.TimeoutExpired as e:
        raise TimeoutError(f"atdyn {task_description} timed out after {timeout} seconds") from e

    if proc.returncode != 0:
        stderr_text = proc.stderr.decode('utf-8', errors='replace')
        raise RuntimeError(
            f"atdyn subprocess failed with code {proc.returncode}:\n{stderr_text}"
        )

    try:
        output = pickle.loads(base64.b64decode(proc.stdout))
    except Exception as e:
        raise RuntimeError(
            f"Failed to decode subprocess output: {e}\n"
            f"stdout: {proc.stdout[:500]}\n"
            f"stderr: {proc.stderr.decode('utf-8', errors='replace')}"
        )

    if not output["success"]:
        from .exceptions import GenesisFortranError, GenesisValidationError
        error_type = output.get("error_type", "")
        error_msg = output.get("error", "Unknown error")

        if "GenesisFortran" in error_type:
            raise GenesisFortranError(error_msg)
        elif error_type == "GenesisValidationError":
            raise GenesisValidationError(error_msg)
        else:
            raise RuntimeError(f"{error_type}: {error_msg}\n{output.get('traceback', '')}")

    return output


def run_atdyn_md_isolated(
    timeout: Optional[float] = None,
    **kwargs
) -> AtdynMDResult:
    """
    Run atdyn MD simulation in an isolated subprocess.

    This function runs the MD simulation in a separate Python subprocess,
    providing maximum isolation from Fortran errors or crashes. Use this
    when running multiple sequential simulations or when reliability is
    more important than performance.

    Args:
        timeout: Maximum time in seconds to wait for completion (None = no limit)
        **kwargs: All arguments passed to run_atdyn_md()

    Returns:
        AtdynMDResult containing energies, final coordinates, and labels

    Raises:
        GenesisValidationError: If input validation fails
        GenesisFortranError: If Fortran code returns an error
        TimeoutError: If simulation exceeds timeout
        RuntimeError: If subprocess fails unexpectedly

    Example:
        >>> # Run multiple simulations safely
        >>> for i in range(10):
        ...     result = run_atdyn_md_isolated(
        ...         prmtopfile="system.prmtop",
        ...         ambcrdfile=f"frame_{i}.rst",
        ...         nsteps=1000,
        ...         timeout=300,  # 5 minute timeout
        ...     )
    """
    output = _run_isolated_subprocess(
        func_name="run_atdyn_md",
        result_fields=["energies", "final_coords", "energy_labels"],
        timeout=timeout,
        task_description="MD simulation",
        **kwargs
    )
    return AtdynMDResult(
        energies=np.array(output["energies"]),
        final_coords=np.array(output["final_coords"]),
        energy_labels=output["energy_labels"],
    )


def run_atdyn_min_isolated(
    timeout: Optional[float] = None,
    **kwargs
) -> AtdynMinResult:
    """
    Run atdyn energy minimization in an isolated subprocess.

    This function runs the minimization in a separate Python subprocess,
    providing maximum isolation from Fortran errors or crashes.

    Args:
        timeout: Maximum time in seconds to wait for completion (None = no limit)
        **kwargs: All arguments passed to run_atdyn_min()

    Returns:
        AtdynMinResult containing energies, final coordinates, convergence info

    Raises:
        GenesisValidationError: If input validation fails
        GenesisFortranError: If Fortran code returns an error
        TimeoutError: If minimization exceeds timeout
        RuntimeError: If subprocess fails unexpectedly
    """
    output = _run_isolated_subprocess(
        func_name="run_atdyn_min",
        result_fields=["energies", "final_coords", "converged", "final_gradient", "energy_labels"],
        timeout=timeout,
        task_description="minimization",
        **kwargs
    )
    return AtdynMinResult(
        energies=np.array(output["energies"]),
        final_coords=np.array(output["final_coords"]),
        converged=output["converged"],
        final_gradient=output["final_gradient"],
        energy_labels=output["energy_labels"],
    )
