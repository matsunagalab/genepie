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
from .s_trajectories import STrajectories
from . import ctrl_files
from . import c2py_util
from functools import lru_cache
from .exceptions import (
    GenesisFortranError,
    GenesisValidationError,
    raise_fortran_error,
)
from .output_capture import suppress_stdout_capture_stderr
from .validation import validate_positive, validate_non_negative, validate_trajectory_dimensions

_DEFAULT_MSG_LEN = 2048


# Trajectory format constants (from trajectory_str.fpp)
class TrjFormat:
    PDB = 1
    AMBER = 2
    DCD = 3
    GROMACS = 4
    CHARMM_RST = 5
    NAMD_RST = 6


# Trajectory type constants
class TrjType:
    COOR = 1
    COOR_BOX = 2


# Fitting method constants (from fitting_str.fpp)
class FittingMethod:
    NO = 1
    TR_ROT = 2
    TR = 3
    TR_ZROT = 4
    XYTR = 5
    XYTR_ZROT = 6


# PBC correction mode constants (from pbc_correct.fpp)
class PBCCMode:
    NO = 1
    MOLECULE = 2


# Map string names to constants
_TRJ_FORMAT_MAP = {
    "PDB": TrjFormat.PDB,
    "AMBER": TrjFormat.AMBER,
    "DCD": TrjFormat.DCD,
    "GROMACS": TrjFormat.GROMACS,
    "CHARMM_RST": TrjFormat.CHARMM_RST,
    "NAMD_RST": TrjFormat.NAMD_RST,
}

_TRJ_TYPE_MAP = {
    "COOR": TrjType.COOR,
    "COOR+BOX": TrjType.COOR_BOX,
}

_FITTING_METHOD_MAP = {
    "NO": FittingMethod.NO,
    "TR+ROT": FittingMethod.TR_ROT,
    "TR": FittingMethod.TR,
    "TR+ZROT": FittingMethod.TR_ZROT,
    "XYTR": FittingMethod.XYTR,
    "XYTR+ZROT": FittingMethod.XYTR_ZROT,
}

_PBCC_MODE_MAP = {
    "NO": PBCCMode.NO,
    "MOLECULE": PBCCMode.MOLECULE,
}


CrdConvertInfo = namedtuple('CrdConvertInfo', [
    'frame_counts',           # List[int] - frame counts per trajectory file
    'selected_atom_indices',  # np.ndarray - selected atom indices (1-indexed)
    'num_selected_atoms',     # int - number of selected atoms
])


def _pack_filenames(filenames: List[str]) -> Tuple[bytes, int, int]:
    """Pack list of filenames into fixed-width byte buffer for Fortran.

    Args:
        filenames: List of file path strings

    Returns:
        Tuple of (packed_bytes, n_files, max_filename_len)
    """
    if not filenames:
        return b'', 0, 0
    max_len = max(len(f) for f in filenames)
    # Pad each filename to max_len with null bytes
    packed = b''.join(f.encode('utf-8').ljust(max_len, b'\x00') for f in filenames)
    return packed, len(filenames), max_len


def crd_convert_info(
    molecule: SMolecule,
    trj_files: List[str],
    trj_format: str = "DCD",
    trj_type: str = "COOR+BOX",
) -> CrdConvertInfo:
    """Get trajectory info (frame counts) for zerocopy crd_convert.

    This is Phase 1 of the zerocopy crd_convert pattern. It reads trajectory
    headers to determine frame counts, allowing Python to pre-allocate arrays.

    Args:
        molecule: SMolecule object containing molecular structure
        trj_files: List of trajectory file paths
        trj_format: Trajectory format ("DCD", "AMBER", "PDB", etc.)
        trj_type: Trajectory type ("COOR" or "COOR+BOX")

    Returns:
        CrdConvertInfo namedtuple with:
            - frame_counts: List[int] of frame counts per trajectory file
            - selected_atom_indices: np.ndarray (empty, for API compatibility)
            - num_selected_atoms: 0 (selection is deferred to crd_convert)

    Raises:
        GenesisValidationError: If parameters are invalid
        GenesisFortranError: If Fortran returns an error

    Note:
        Frame count auto-detection only works for DCD format. For other formats,
        frame counts are estimated and may not be exact.
    """
    lib = LibGenesis().lib
    mol_c = molecule.to_SMoleculeC()

    # Validate and convert format/type strings to constants
    fmt_upper = trj_format.upper()
    if fmt_upper not in _TRJ_FORMAT_MAP:
        raise GenesisValidationError(
            f"Invalid trajectory format: {trj_format}. "
            f"Valid formats: {list(_TRJ_FORMAT_MAP.keys())}"
        )
    trj_format_c = _TRJ_FORMAT_MAP[fmt_upper]

    type_upper = trj_type.upper()
    if type_upper not in _TRJ_TYPE_MAP:
        raise GenesisValidationError(
            f"Invalid trajectory type: {trj_type}. "
            f"Valid types: {list(_TRJ_TYPE_MAP.keys())}"
        )
    trj_type_c = _TRJ_TYPE_MAP[type_upper]

    # Pack filenames
    packed_names, n_files, max_len = _pack_filenames(trj_files)

    # Prepare output variables
    frame_counts_ptr = ctypes.c_void_p()
    n_trajs = ctypes.c_int()
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    try:
        with suppress_stdout_capture_stderr() as captured:
            lib.crd_convert_info_c(
                ctypes.byref(mol_c),
                packed_names,
                ctypes.c_int(n_files),
                ctypes.c_int(max_len),
                ctypes.c_int(trj_format_c),
                ctypes.c_int(trj_type_c),
                ctypes.byref(frame_counts_ptr),
                ctypes.byref(n_trajs),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen)
            )

        if status.value != 0:
            error_msg = msg.value.decode('utf-8', errors='replace').strip()
            raise_fortran_error(
                status.value,
                error_msg,
                stderr_output=captured.stderr
            )

        # Convert frame counts to Python list
        if n_trajs.value > 0 and frame_counts_ptr.value:
            arr_ptr = ctypes.cast(
                frame_counts_ptr.value,
                ctypes.POINTER(ctypes.c_int)
            )
            frame_counts = [arr_ptr[i] for i in range(n_trajs.value)]
            # Deallocate Fortran-allocated memory
            lib.deallocate_frame_counts_c(frame_counts_ptr)
        else:
            frame_counts = []

        return CrdConvertInfo(
            frame_counts=frame_counts,
            selected_atom_indices=np.array([], dtype=np.int32),
            num_selected_atoms=0,
        )

    finally:
        lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


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


# Alias for selection function to avoid name collision with parameter
selection_func = selection


def crd_convert(
    molecule: SMolecule,
    trj_files: List[str],
    trj_format: str = "DCD",
    trj_type: str = "COOR+BOX",
    selection: str = "all",
    fitting_selection: Optional[str] = None,
    fitting_method: str = "NO",
    mass_weighted: bool = False,
    centering: bool = False,
    centering_selection: Optional[str] = None,
    center_coord: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    pbc_correct: str = "NO",
    ana_period: int = 1,
    rename_res: Optional[List[str]] = None,
) -> Tuple[List[STrajectories], SMolecule]:
    """Load and convert trajectory files using zerocopy pattern.

    This function reads trajectory files and converts them into STrajectories
    objects that can be used for analysis. It uses zerocopy to avoid unnecessary
    memory copies between Python and Fortran.

    Args:
        molecule: SMolecule object containing molecular structure
        trj_files: List of trajectory file paths
        trj_format: Trajectory format ("DCD", "AMBER", "PDB", "GROMACS", etc.)
        trj_type: Trajectory type ("COOR" or "COOR+BOX")
        selection: Atom selection string (e.g., "all", "an:CA", "rno:1-10")
        fitting_selection: Selection for fitting atoms (None for no fitting)
        fitting_method: Fitting method ("NO", "TR+ROT", "TR", "TR+ZROT", etc.)
        mass_weighted: Use mass weighting for fitting
        centering: Enable coordinate centering
        centering_selection: Selection for centering atoms (default: same as selection)
        center_coord: Target center coordinates as (x, y, z)
        pbc_correct: PBC correction mode ("NO", "MOLECULE")
        ana_period: Analysis period (process every Nth frame)
        rename_res: List of residue name mappings (e.g., ["HSE HIS", "HSD HIS"])
                    Each string should be "FROM TO" format

    Returns:
        Tuple of (List[STrajectories], SMolecule) where:
            - List[STrajectories]: One STrajectories per input trajectory file
            - SMolecule: Subset molecule containing only selected atoms

    Raises:
        GenesisValidationError: If parameters are invalid
        GenesisFortranError: If trajectory reading fails

    Examples:
        >>> trajs, mol = crd_convert(mol, ["traj.dcd"])
        >>> trajs, mol = crd_convert(mol, ["traj.dcd"], selection="an:CA")
        >>> trajs, mol = crd_convert(mol, ["traj.dcd"],
        ...                          fitting_selection="an:CA",
        ...                          fitting_method="TR+ROT")
        >>> trajs, mol = crd_convert(mol, ["traj.dcd"],
        ...                          rename_res=["HSE HIS", "HSD HIS"])
    """
    lib = LibGenesis().lib

    # Apply residue name renaming if specified
    if rename_res:
        for rename_spec in rename_res:
            parts = rename_spec.split()
            if len(parts) != 2:
                raise GenesisValidationError(
                    f"Invalid rename_res format: '{rename_spec}'. "
                    f"Expected 'FROM TO' format (e.g., 'HSE HIS')"
                )
            from_name, to_name = parts
            # Rename matching residues in the molecule
            mask = molecule.residue_name == from_name
            if np.any(mask):
                molecule.residue_name[mask] = to_name

    # Validate and convert format/type strings to constants
    fmt_upper = trj_format.upper()
    if fmt_upper not in _TRJ_FORMAT_MAP:
        raise GenesisValidationError(
            f"Invalid trajectory format: {trj_format}. "
            f"Valid formats: {list(_TRJ_FORMAT_MAP.keys())}"
        )
    trj_format_c = _TRJ_FORMAT_MAP[fmt_upper]

    type_upper = trj_type.upper()
    if type_upper not in _TRJ_TYPE_MAP:
        raise GenesisValidationError(
            f"Invalid trajectory type: {trj_type}. "
            f"Valid types: {list(_TRJ_TYPE_MAP.keys())}"
        )
    trj_type_c = _TRJ_TYPE_MAP[type_upper]

    fit_method_upper = fitting_method.upper()
    if fit_method_upper not in _FITTING_METHOD_MAP:
        raise GenesisValidationError(
            f"Invalid fitting method: {fitting_method}. "
            f"Valid methods: {list(_FITTING_METHOD_MAP.keys())}"
        )
    fitting_method_c = _FITTING_METHOD_MAP[fit_method_upper]

    pbcc_upper = pbc_correct.upper()
    if pbcc_upper not in _PBCC_MODE_MAP:
        raise GenesisValidationError(
            f"Invalid PBC correction mode: {pbc_correct}. "
            f"Valid modes: {list(_PBCC_MODE_MAP.keys())}"
        )
    pbcc_mode_c = _PBCC_MODE_MAP[pbcc_upper]

    # Phase 1: Get trajectory info (frame counts)
    info = crd_convert_info(molecule, trj_files, trj_format, trj_type)
    frame_counts = info.frame_counts

    if not frame_counts or all(fc == 0 for fc in frame_counts):
        raise GenesisValidationError("No frames found in trajectory files")

    # Get selected atom indices
    selected_indices = selection_func(molecule, selection)
    n_selected = len(selected_indices)

    if n_selected == 0:
        raise GenesisValidationError(
            f"Selection '{selection}' returned no atoms"
        )

    # Get fitting atom indices if fitting is requested
    if fitting_selection is not None and fitting_method_c != FittingMethod.NO:
        fitting_indices = selection_func(molecule, fitting_selection)
        n_fitting = len(fitting_indices)
    else:
        fitting_indices = np.array([], dtype=np.int32)
        n_fitting = 0

    # Get centering atom indices if centering is requested
    if centering and centering_selection is not None:
        centering_indices = selection_func(molecule, centering_selection)
        n_centering = len(centering_indices)
    elif centering:
        # Default: use same atoms as selection
        centering_indices = selected_indices.copy()
        n_centering = len(centering_indices)
    else:
        centering_indices = np.array([], dtype=np.int32)
        n_centering = 0

    # Phase 2: Pre-allocate numpy arrays for trajectory data
    n_trajs = len(frame_counts)

    # Calculate actual frame counts after applying ana_period
    actual_frame_counts = [
        (fc + ana_period - 1) // ana_period for fc in frame_counts
    ]

    # Allocate coordinate arrays in Fortran order (column-major)
    # Fortran expects shape (3, n_selected, n_frames) with column-major layout
    coords_list = []
    pbc_box_list = []
    for i, n_frames in enumerate(actual_frame_counts):
        # Allocate as Fortran-contiguous (3, n_selected, n_frames)
        coords = np.zeros((3, n_selected, n_frames), dtype=np.float64, order='F')
        coords_list.append(coords)
        if trj_type_c == TrjType.COOR_BOX:
            # PBC box: (3, 3, n_frames) in Fortran order
            pbc_box = np.zeros((3, 3, n_frames), dtype=np.float64, order='F')
        else:
            pbc_box = np.zeros((0, 0, 0), dtype=np.float64, order='F')
        pbc_box_list.append(pbc_box)

    # Create array of pointers for Fortran
    coords_ptrs = (ctypes.c_void_p * n_trajs)()
    pbc_box_ptrs = (ctypes.c_void_p * n_trajs)()
    for i in range(n_trajs):
        coords_ptrs[i] = coords_list[i].ctypes.data
        pbc_box_ptrs[i] = pbc_box_list[i].ctypes.data if pbc_box_list[i].size > 0 else None

    # Pack filenames and prepare data for Fortran call
    packed_names, n_files, max_len = _pack_filenames(trj_files)

    # Prepare frame counts array
    frame_counts_arr = np.array(frame_counts, dtype=np.int32)

    # Prepare center coordinate
    center_coord_arr = np.array(center_coord, dtype=np.float64)

    # Prepare molecule C structure
    mol_c = molecule.to_SMoleculeC()

    # Prepare output variables
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    try:
        with suppress_stdout_capture_stderr() as captured:
            lib.crd_convert_zerocopy_c(
                ctypes.byref(mol_c),
                packed_names,
                ctypes.c_int(n_files),
                ctypes.c_int(max_len),
                ctypes.c_int(trj_format_c),
                ctypes.c_int(trj_type_c),
                selected_indices.ctypes.data,
                ctypes.c_int(n_selected),
                ctypes.c_int(fitting_method_c),
                fitting_indices.ctypes.data if n_fitting > 0 else None,
                ctypes.c_int(n_fitting),
                ctypes.c_int(1 if mass_weighted else 0),
                ctypes.c_int(1 if centering else 0),
                centering_indices.ctypes.data if n_centering > 0 else None,
                ctypes.c_int(n_centering),
                center_coord_arr.ctypes.data,
                ctypes.c_int(pbcc_mode_c),
                ctypes.c_int(ana_period),
                frame_counts_arr.ctypes.data,
                ctypes.cast(coords_ptrs, ctypes.c_void_p),
                ctypes.cast(pbc_box_ptrs, ctypes.c_void_p),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen)
            )

        if status.value != 0:
            error_msg = msg.value.decode('utf-8', errors='replace').strip()
            raise_fortran_error(
                status.value,
                error_msg,
                stderr_output=captured.stderr
            )

        # Create STrajectories objects from filled arrays
        # Transpose from Fortran order (3, n_selected, n_frames) to
        # Python order (n_frames, n_selected, 3)
        trajectories = []
        for i, n_frames in enumerate(actual_frame_counts):
            # Transpose coords: (3, n_sel, n_frame) -> (n_frame, n_sel, 3)
            coords_py = np.ascontiguousarray(
                np.transpose(coords_list[i], (2, 1, 0))
            )
            if trj_type_c == TrjType.COOR_BOX and pbc_box_list[i].size > 0:
                # Transpose pbc_box: (3, 3, n_frame) -> (n_frame, 3)
                # Extract diagonal (box dimensions)
                pbc_box_py = np.ascontiguousarray(
                    np.transpose(pbc_box_list[i], (2, 0, 1))[:, np.diag_indices(3)[0], np.diag_indices(3)[1]]
                )
            else:
                pbc_box_py = None
            traj = STrajectories.from_numpy(
                coords_py,
                pbc_box=pbc_box_py,
                mem_owner=True
            )
            trajectories.append(traj)

        # Create subset molecule
        # selected_indices is 1-indexed from Fortran, convert to 0-indexed
        subset_mol = molecule.subset_atoms(selected_indices - 1)

        return trajectories, subset_mol

    finally:
        lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


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

    Args:
        trajs: STrajectories object containing trajectory data
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
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

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

        with suppress_stdout_capture_stderr() as captured:
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

        if status.value != 0:
            error_msg = msg.value.decode('utf-8', errors='replace').strip()
            stderr_output = captured.stderr if captured else ""
            raise_fortran_error(status.value, error_msg, stderr_output)

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
        with suppress_stdout_capture_stderr() as captured:
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

        if status.value != 0:
            error_msg = msg.value.decode('utf-8', errors='replace').strip()
            stderr_output = captured.stderr if captured else ""
            raise_fortran_error(status.value, error_msg, stderr_output)

        n_actual = nstru_out.value

        # Transpose to Python convention and slice
        final_distance = result_distance[:, :n_actual].T.copy() if result_distance is not None else None
        final_angle = result_angle[:, :n_actual].T.copy() if result_angle is not None else None
        final_torsion = result_torsion[:, :n_actual].T.copy() if result_torsion is not None else None

        return TrjAnalysisResult(
            final_distance, final_angle, final_torsion,
            None, None, None
        )


RgAnalysisResult = namedtuple(
        'RgAnalysisResult',
        ['rg'])


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
        trajs: Trajectories to analyze
        analysis_selection: GENESIS selection string (e.g., "an:CA", "heavy")
        ana_period: Analysis period (default: 1)
        mass_weighted: Use mass weighting for RG calculation (default: True)

    Returns:
        RgAnalysisResult containing the radius of gyration array

    Example:
        >>> result = rg_analysis(mol, trajs, "an:CA")
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

    # Pre-allocate result array (zero-copy)
    n_frame = int(trajs.nframe / ana_period)
    result_rg = np.zeros(n_frame, dtype=np.float64)
    result_ptr = result_rg.ctypes.data_as(ctypes.c_void_p)

    # Output variables
    nstru_out = ctypes.c_int()
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    with suppress_stdout_capture_stderr() as captured:
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

    # Check for errors
    if status.value != 0:
        error_msg = msg.value.decode('utf-8', errors='replace').strip()
        stderr_output = captured.stderr if captured else ""
        raise_fortran_error(status.value, error_msg, stderr_output)

    # Return the pre-allocated result array (already filled by Fortran)
    return RgAnalysisResult(result_rg[:nstru_out.value])


RmsdAnalysisResult = namedtuple(
        'RmsdAnalysisResult',
        ['rmsd'])


# Fitting method constants
class FittingMethod:
    """Fitting method constants for RMSD analysis with fitting."""
    NO = 1
    TR_ROT = 2
    TR = 3
    TR_ZROT = 4
    XYTR = 5
    XYTR_ZROT = 6


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
    """
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
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    if fitting_selection is None:
        # No fitting - call simple RMSD analysis
        with suppress_stdout_capture_stderr() as captured:
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
        # Fitting method mapping
        method_map = {
            "NO": FittingMethod.NO,
            "TR+ROT": FittingMethod.TR_ROT,
            "TR": FittingMethod.TR,
            "TR+ZROT": FittingMethod.TR_ZROT,
            "XYTR": FittingMethod.XYTR,
            "XYTR+ZROT": FittingMethod.XYTR_ZROT,
        }
        if fitting_method not in method_map:
            raise GenesisValidationError(
                f"Invalid fitting_method: {fitting_method}. "
                f"Valid options: {list(method_map.keys())}"
            )
        method_int = method_map[fitting_method]

        # Get fitting indices
        fitting_indices = selection(molecule, fitting_selection)
        n_fitting = len(fitting_indices)

        with suppress_stdout_capture_stderr() as captured:
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

    # Check for errors
    if status.value != 0:
        error_msg = msg.value.decode('utf-8', errors='replace').strip()
        stderr_output = captured.stderr if captured else ""
        raise_fortran_error(status.value, error_msg, stderr_output)

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
        max_frames: Maximum expected frames (for result array allocation, default: 100000)

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
    lib = LibGenesis().lib

    # Validate DCD file exists
    if not os.path.exists(dcd_file):
        raise GenesisValidationError(f"DCD file not found: {dcd_file}")

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

    # Pre-allocate result array (will be trimmed later)
    result_rmsd = np.zeros(max_frames, dtype=np.float64)
    result_ptr = result_rmsd.ctypes.data_as(ctypes.c_void_p)

    # Convert filename to C string
    dcd_filename_bytes = dcd_file.encode('utf-8')
    filename_len = len(dcd_filename_bytes)

    # Trajectory type: 1=COOR, 2=COOR+BOX (matches TrjTypeCoor=1, TrjTypeCoorBox=2)
    trj_type = 2 if has_box else 1

    # Output variables
    nstru_out = ctypes.c_int()
    dcd_nframe_out = ctypes.c_int()
    dcd_natom_out = ctypes.c_int()
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    # Fitting parameters
    fitting_idx_ptr = ctypes.c_void_p(0)
    n_fitting = 0
    method_int = 0

    if fitting_selection is not None:
        # Fitting method mapping
        method_map = {
            "NO": FittingMethod.NO,
            "TR+ROT": FittingMethod.TR_ROT,
            "TR": FittingMethod.TR,
            "TR+ZROT": FittingMethod.TR_ZROT,
            "XYTR": FittingMethod.XYTR,
            "XYTR+ZROT": FittingMethod.XYTR_ZROT,
        }
        if fitting_method not in method_map:
            raise GenesisValidationError(
                f"Invalid fitting_method: {fitting_method}. "
                f"Valid options: {list(method_map.keys())}"
            )
        method_int = method_map[fitting_method]

        # Get fitting indices
        fitting_indices = selection(molecule, fitting_selection)
        n_fitting = len(fitting_indices)
        fitting_idx_ptr = fitting_indices.ctypes.data_as(ctypes.c_void_p)

    with suppress_stdout_capture_stderr() as captured:
        lib.rmsd_analysis_lazy_c(
            dcd_filename_bytes,
            ctypes.c_int(filename_len),
            ctypes.c_int(trj_type),
            mass_ptr,
            ref_coord_ptr,
            ctypes.c_int(molecule.num_atoms),
            ctypes.c_int(ana_period),
            fitting_idx_ptr,
            ctypes.c_int(n_fitting),
            analysis_indices.ctypes.data_as(ctypes.c_void_p),
            ctypes.c_int(n_analysis),
            ctypes.c_int(method_int),
            ctypes.c_int(1 if mass_weighted else 0),
            result_ptr,
            ctypes.c_int(max_frames),
            ctypes.byref(nstru_out),
            ctypes.byref(dcd_nframe_out),
            ctypes.byref(dcd_natom_out),
            ctypes.byref(status),
            msg,
            ctypes.c_int(msglen),
        )

    # Check for errors
    if status.value != 0:
        error_msg = msg.value.decode('utf-8', errors='replace').strip()
        stderr_output = captured.stderr if captured else ""
        raise_fortran_error(status.value, error_msg, stderr_output)

    # Return the result trimmed to actual size
    return RmsdLazyAnalysisResult(
        result_rmsd[:nstru_out.value],
        dcd_nframe_out.value,
        dcd_natom_out.value
    )


DrmsAnalysisResult = namedtuple(
        'DrmsAnalysisResult',
        ['drms'])


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
        trajs: Trajectories to analyze
        contact_list: Contact atom pairs as (2, n_contact) array with 1-indexed
                      atom indices. Each column is [atom1_idx, atom2_idx].
        contact_dist: Reference distances for each contact pair (n_contact,)
        ana_period: Analysis period (default: 1)
        pbc_correct: Apply PBC correction for distances (default: False)

    Returns:
        DrmsAnalysisResult containing the DRMS array

    Example:
        >>> contact_list = np.array([[1, 2, 3], [10, 11, 12]], dtype=np.int32)
        >>> contact_dist = np.array([5.0, 6.0, 7.0], dtype=np.float64)
        >>> result = drms_analysis(trajs, contact_list, contact_dist)
        >>> print(result.drms)
    """
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
    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    with suppress_stdout_capture_stderr() as captured:
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

    # Check for errors
    if status.value != 0:
        error_msg = msg.value.decode('utf-8', errors='replace').strip()
        stderr_output = captured.stderr if captured else ""
        raise_fortran_error(status.value, error_msg, stderr_output)

    # Return result sliced to actual size
    return DrmsAnalysisResult(result_drms[:nstru_out.value])


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


DiffusionAnalysisResult = namedtuple(
        'DiffusionAnalysisResult',
        ['out_data', 'diffusion_coefficients'])


def diffusion_analysis(
        msd_data: npt.NDArray[np.float64],
        time_step: float = 1.0,
        distance_unit: float = 1.0,
        ndofs: int = 3,
        start_step: int = 1,
        stop_step: Optional[int] = None,
        ) -> DiffusionAnalysisResult:
    """
    Executes diffusion analysis on mean square displacement data.

    This function analyzes MSD data to calculate diffusion coefficients
    using linear fitting. It uses zero-copy memory sharing with pre-allocated
    result arrays for optimal performance.

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
        >>> result = diffusion_analysis(msd, time_step=0.01)
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

    # Pre-allocate output arrays (full zero-copy)
    out_ncols = 2 * n_sets + 1  # time + (msd + fit) * n_sets
    out_data_f = np.zeros((out_ncols, ndata), dtype=np.float64, order='F')
    diff_coeff = np.zeros(n_sets, dtype=np.float64)

    out_data_ptr = out_data_f.ctypes.data_as(ctypes.c_void_p)
    diff_coeff_ptr = diff_coeff.ctypes.data_as(ctypes.c_void_p)

    status = ctypes.c_int()
    msglen = _DEFAULT_MSG_LEN
    msg = ctypes.create_string_buffer(msglen)

    with suppress_stdout_capture_stderr() as captured:
        lib.diffusion_analysis_c(
            msd_ptr,
            ctypes.c_int(ndata),
            ctypes.c_int(ncols),
            ctypes.c_double(time_step),
            ctypes.c_double(distance_unit),
            ctypes.c_int(ndofs),
            ctypes.c_int(start_step),
            ctypes.c_int(stop_step),
            out_data_ptr,
            ctypes.c_int(out_ncols * ndata),
            diff_coeff_ptr,
            ctypes.c_int(n_sets),
            ctypes.byref(status),
            msg,
            ctypes.c_int(msglen),
        )

    # Check for errors
    if status.value != 0:
        error_msg = msg.value.decode('utf-8', errors='replace').strip()
        stderr_output = captured.stderr if captured else ""
        raise_fortran_error(status.value, error_msg, stderr_output)

    # Transpose to (ndata, ncols) for Python convention
    out_data = out_data_f.T

    return DiffusionAnalysisResult(out_data, diff_coeff)


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
