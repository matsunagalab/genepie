"""Trajectory loading, atom selection and coordinate conversion."""

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
from ..constants import (
    FittingMethod,
    PBCCMode,
    TrjFormat,
    TrjType,
)
from ..exceptions import GenesisValidationError
from .._fortran import fortran_status
from ._common import (
    _resolve_enum,
    _pack_filenames,
    _TRJ_FORMAT_MAP,
    _TRJ_TYPE_MAP,
    _FITTING_METHOD_MAP,
    _PBCC_MODE_MAP,
)


CrdConvertInfo = namedtuple('CrdConvertInfo', [
    'frame_counts',           # List[int] - frame counts per trajectory file
    'selected_atom_indices',  # np.ndarray - selected atom indices (1-indexed)
    'num_selected_atoms',     # int - number of selected atoms
])


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
    trj_format_c = _resolve_enum(
        trj_format, _TRJ_FORMAT_MAP, "trajectory format")
    trj_type_c = _resolve_enum(trj_type, _TRJ_TYPE_MAP, "trajectory type")

    # Pack filenames
    packed_names, n_files, max_len = _pack_filenames(trj_files)

    # Prepare output variables
    frame_counts_ptr = ctypes.c_void_p()
    n_trajs = ctypes.c_int()

    try:
        with fortran_status() as (status, msg, msglen):
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

    try:
        with fortran_status() as (status, msg, msglen):
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
    lazy: bool = False,
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
        lazy: If True, create lazy STrajectories without loading data.
              Lazy mode has restrictions: single DCD file, no fitting,
              no centering, no PBC correction. Analysis functions read
              directly from the file. ``selection`` and ``ana_period`` define
              the logical lazy view exactly as in memory mode.

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
        >>> # Lazy mode: no data loading, analysis reads from file
        >>> lazy_trajs, mol = crd_convert(mol, ["traj.dcd"], lazy=True)
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
    trj_format_c = _resolve_enum(
        trj_format, _TRJ_FORMAT_MAP, "trajectory format")
    trj_type_c = _resolve_enum(trj_type, _TRJ_TYPE_MAP, "trajectory type")
    fitting_method_c = _resolve_enum(
        fitting_method, _FITTING_METHOD_MAP, "fitting method")
    pbcc_mode_c = _resolve_enum(
        pbc_correct, _PBCC_MODE_MAP, "PBC correction mode")
    if ana_period <= 0:
        raise GenesisValidationError(
            f"ana_period must be positive, got {ana_period}"
        )

    # Handle lazy mode: create STrajectories without loading data
    if lazy:
        # Validate lazy mode restrictions
        if len(trj_files) != 1:
            raise GenesisValidationError(
                "lazy=True requires exactly one trajectory file, "
                f"got {len(trj_files)}"
            )
        if fitting_method_c != FittingMethod.NO:
            raise GenesisValidationError(
                "lazy=True does not support fitting. "
                "Fitting requires all coordinates in memory."
            )
        if centering:
            raise GenesisValidationError(
                "lazy=True does not support centering. "
                "Centering requires coordinate modification."
            )
        if pbcc_mode_c != PBCCMode.NO:
            raise GenesisValidationError(
                "lazy=True does not support PBC correction. "
                "PBC correction requires coordinate modification."
            )
        if trj_format_c != TrjFormat.DCD:
            raise GenesisValidationError(
                f"lazy=True only supports DCD format, got {trj_format}"
            )

        # Get trajectory info
        info = crd_convert_info(molecule, trj_files, trj_format, trj_type)
        if not info.frame_counts or info.frame_counts[0] == 0:
            raise GenesisValidationError("No frames found in trajectory file")

        nframe = info.frame_counts[0]
        if nframe // ana_period == 0:
            raise GenesisValidationError(
                f"ana_period={ana_period} selects no frames from {nframe} frames"
            )

        # Get selected atom indices (1-indexed for Fortran)
        selected_indices = selection_func(molecule, selection)
        n_selected = len(selected_indices)

        if n_selected == 0:
            raise GenesisValidationError(
                f"Selection '{selection}' returned no atoms"
            )

        # Create lazy STrajectories
        # Note: natom is set to molecule.num_atoms (DCD file atom count)
        # selection_indices tells which atoms to use during analysis
        from ..s_trajectories import TRJ_TYPE_COOR, TRJ_TYPE_COOR_BOX
        lazy_trj_type = TRJ_TYPE_COOR_BOX if trj_type_c == TrjType.COOR_BOX else TRJ_TYPE_COOR

        lazy_traj = STrajectories.from_lazy(
            dcd_file=trj_files[0],
            trj_type=lazy_trj_type,
            nframe=nframe,
            natom=molecule.num_atoms,  # Full atom count from DCD
            selection_indices=selected_indices,
            ana_period=ana_period,
        )

        # Create subset molecule
        subset_mol = molecule.subset_atoms(selected_indices - 1)

        return [lazy_traj], subset_mol

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
    actual_frame_counts = [fc // ana_period for fc in frame_counts]
    if any(fc == 0 for fc in actual_frame_counts):
        raise GenesisValidationError(
            f"ana_period={ana_period} selects no frames from at least one "
            "trajectory"
        )

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

    try:
        with fortran_status() as (status, msg, msglen):
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
