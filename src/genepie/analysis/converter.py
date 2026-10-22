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
        Frame counts are auto-detected for DCD (from the header) and PDB
        (frames delimited by ENDMDL, or END, records; a bare list of ATOM
        records counts as one frame). Other formats are not supported yet.
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
        trj_format: Trajectory format ("DCD" or "PDB"; other formats cannot
                    auto-detect their frame count yet)
        trj_type: Trajectory type ("COOR" or "COOR+BOX"). For PDB input the
                  box is taken from the CRYST1 record.
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
        >>> # A single PDB structure becomes a one-frame trajectory
        >>> trajs, mol = crd_convert(mol, ["model.pdb"], trj_format="PDB",
        ...                          trj_type="COOR")
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
        # PBC box: (3, 3, n_frames) in Fortran order. Fortran writes the box
        # of every frame unconditionally, so the buffer is allocated even for
        # trj_type="COOR"; it is simply discarded below in that case.
        pbc_box = np.zeros((3, 3, n_frames), dtype=np.float64, order='F')
        pbc_box_list.append(pbc_box)

    # Create array of pointers for Fortran
    coords_ptrs = (ctypes.c_void_p * n_trajs)()
    pbc_box_ptrs = (ctypes.c_void_p * n_trajs)()
    for i in range(n_trajs):
        coords_ptrs[i] = coords_list[i].ctypes.data
        pbc_box_ptrs[i] = pbc_box_list[i].ctypes.data

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
            if trj_type_c == TrjType.COOR_BOX:
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


# Fortran rc_option_str_mod: REPLICA = 1, PARAMETER = 2.
_REMD_CONVERT_TYPE = {"REPLICA": 1, "PARAMETER": 2}

RemdConvertResult = namedtuple("RemdConvertResult", [
    "trajectories",  # List[STrajectories], one per parameter or replica id
    "subset",        # SMolecule reduced to the selection
    "param_ids",     # List[int], 1-based ids in output order
])


def remd_convert(
    molecule: SMolecule,
    dcd_files: List[str],
    rem_files: Optional[List[str]] = None,
    convert_type: str = "PARAMETER",
    nsteps: int = 0,
    exchange_period: int = 0,
    crdout_period: int = 0,
    trjout_period: int = 0,
    convert_ids: Optional[List[int]] = None,
    selection: str = "all",
) -> RemdConvertResult:
    """Sort replica DCD files onto parameter or replica ids.

    This is the in-memory form of the ``remd_convert`` program for
    coordinates. The step clock is the same one the program uses: a frame
    is read every ``crdout_period`` MD steps, and the parameter id is
    advanced every ``exchange_period`` steps after that frame is placed.
    ``trjout_period`` defaults to ``crdout_period``.

    Centering, fitting, PBC correction and energy-log sorting stay on the
    ``remd_convert`` command-line program.

    Args:
        molecule: Topology of the atoms stored in the DCD files.
        dcd_files: Replica trajectories, in replica order starting at 1.
        rem_files: One ``.rem`` file per replica. Required for
            ``convert_type="PARAMETER"``.
        convert_type: ``"PARAMETER"`` or ``"REPLICA"``.
        nsteps: ``nsteps`` from the REMD ``[DYNAMICS]`` section.
        exchange_period: ``exchange_period`` from ``[REMD]``.
        crdout_period: ``crdout_period`` from ``[DYNAMICS]``.
        trjout_period: Keep every frame whose MD step is a multiple of this.
            ``0`` means the same value as ``crdout_period``.
        convert_ids: Parameter or replica ids to keep. ``None`` keeps all.
        selection: GENESIS selection string. ``"all"`` keeps every atom.

    Returns:
        RemdConvertResult with one trajectory per requested id, the molecule
        reduced to ``selection``, and those ids in output order.
    """
    lib = LibGenesis().lib
    kind = _resolve_enum(convert_type, _REMD_CONVERT_TYPE, "convert_type")
    if nsteps <= 0 or exchange_period <= 0 or crdout_period <= 0:
        raise GenesisValidationError(
            "nsteps, exchange_period and crdout_period must be positive"
        )
    if trjout_period < 0:
        raise GenesisValidationError(
            f"trjout_period must be >= 0, got {trjout_period}"
        )
    if not dcd_files:
        raise GenesisValidationError("dcd_files is empty")
    nreplica = len(dcd_files)
    if kind == _REMD_CONVERT_TYPE["PARAMETER"]:
        if not rem_files or len(rem_files) != nreplica:
            raise GenesisValidationError(
                "PARAMETER sort needs one rem file per replica trajectory"
            )
        rem_list = list(rem_files)
    else:
        rem_list = ["."] * nreplica

    if convert_ids is None:
        param_ids = list(range(1, nreplica + 1))
        n_ids = 0
        ids_ptr = None
    else:
        param_ids = [int(i) for i in convert_ids]
        if not param_ids:
            raise GenesisValidationError("convert_ids is empty")
        if any(i < 1 or i > nreplica for i in param_ids):
            raise GenesisValidationError(
                f"convert_ids must be in 1..{nreplica}, got {param_ids}"
            )
        if len(set(param_ids)) != len(param_ids):
            raise GenesisValidationError("convert_ids contains a duplicate")
        n_ids = len(param_ids)
        ids_arr = np.array(param_ids, dtype=np.int32)
        ids_ptr = ids_arr.ctypes.data

    selected = np.ascontiguousarray(selection_func(molecule, selection), dtype=np.int32)
    if selected.size == 0:
        raise GenesisValidationError(
            f"Selection '{selection}' returned no atoms"
        )

    width = max(max(len(p) for p in dcd_files), max(len(p) for p in rem_list))
    def _pack(paths: List[str]) -> bytes:
        return b"".join(p.encode("utf-8").ljust(width, b"\0") for p in paths)
    dcd_packed = _pack(list(dcd_files))
    rem_packed = _pack(rem_list)
    n_out = len(param_ids)
    counts = np.zeros(n_out, dtype=np.int32)

    with fortran_status() as (status, msg, msglen):
        lib.remd_convert_counts_c(
            rem_packed,
            ctypes.c_int(nreplica),
            ctypes.c_int(width),
            ctypes.c_int(kind),
            ids_ptr,
            ctypes.c_int(n_ids),
            ctypes.c_int(nsteps),
            ctypes.c_int(exchange_period),
            ctypes.c_int(crdout_period),
            ctypes.c_int(trjout_period),
            counts.ctypes.data,
            ctypes.c_int(n_out),
            ctypes.byref(status),
            msg,
            ctypes.c_int(msglen),
        )

    if int(counts.sum()) == 0:
        raise GenesisValidationError(
            "remd_convert produced no frames. Check nsteps and the periods."
        )

    n_selected = int(selected.size)
    coords_list = []
    box_list = []
    coords_ptrs = (ctypes.c_void_p * n_out)()
    box_ptrs = (ctypes.c_void_p * n_out)()
    for i, nframe in enumerate(int(n) for n in counts):
        coords = np.zeros((3, n_selected, nframe), dtype=np.float64, order="F")
        box = np.zeros((3, 3, nframe), dtype=np.float64, order="F")
        coords_list.append(coords)
        box_list.append(box)
        if nframe > 0:
            coords_ptrs[i] = coords.ctypes.data
            box_ptrs[i] = box.ctypes.data

    with fortran_status() as (status, msg, msglen):
        lib.remd_convert_fill_c(
            dcd_packed,
            rem_packed,
            ctypes.c_int(nreplica),
            ctypes.c_int(width),
            ctypes.c_int(kind),
            ids_ptr,
            ctypes.c_int(n_ids),
            ctypes.c_int(nsteps),
            ctypes.c_int(exchange_period),
            ctypes.c_int(crdout_period),
            ctypes.c_int(trjout_period),
            selected.ctypes.data,
            ctypes.c_int(n_selected),
            ctypes.c_int(molecule.num_atoms),
            counts.ctypes.data,
            ctypes.c_int(n_out),
            ctypes.cast(coords_ptrs, ctypes.c_void_p),
            ctypes.cast(box_ptrs, ctypes.c_void_p),
            ctypes.byref(status),
            msg,
            ctypes.c_int(msglen),
        )

    trajectories = []
    for coords, box in zip(coords_list, box_list):
        coords_py = np.ascontiguousarray(np.transpose(coords, (2, 1, 0)))
        if box.shape[2] == 0:
            pbc_py = np.zeros((0, 3), dtype=np.float64)
        else:
            pbc_py = np.ascontiguousarray(
                np.transpose(box, (2, 0, 1))[:, np.diag_indices(3)[0],
                                             np.diag_indices(3)[1]]
            )
        trajectories.append(STrajectories.from_numpy(
            coords_py, pbc_box=pbc_py, mem_owner=True))

    subset = molecule.subset_atoms(selected.astype(np.int64) - 1)
    return RemdConvertResult(trajectories, subset, param_ids)
