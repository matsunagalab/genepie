import ctypes
import warnings
from typing import Self, Optional
import numpy as np
import numpy.typing as npt
from .s_trajectories_c import STrajectoriesC

# Trajectory type constants (matching Fortran)
TRJ_TYPE_COOR = 1
TRJ_TYPE_COOR_BOX = 2
from .libgenesis import LibGenesis
from .s_molecule import SMolecule
from .exceptions import GenesisValidationError, GenesisMemoryError
from .validation import validate_trajectory_dimensions, validate_non_negative


class STrajectories:
    natom: int
    nframe: int
    coords: npt.NDArray[np.float64]  # shape=(n_frame, n_atom, 3)
    pbc_boxes: npt.NDArray[np.float64]  # shape=(trajs_c.nframe, 3, 3)
    c_obj: STrajectoriesC

    # Lazy loading fields (Python-only, not in C structure)
    is_lazy: bool
    lazy_dcd_file: Optional[str]
    lazy_trj_type: int
    selection_indices: Optional[npt.NDArray[np.int32]]  # 1-indexed atom indices

    def __init__(self, natom: int = 0, nframe: int = 0,
                 trajs_c: STrajectoriesC = None, mem_owner: bool = True):
        # Validate dimensions before allocation
        if trajs_c is None:
            validate_non_negative(natom, "natom")
            validate_non_negative(nframe, "nframe")
            if natom > 0 and nframe > 0:
                validate_trajectory_dimensions(nframe, natom)
        else:
            validate_non_negative(trajs_c.natom, "trajs_c.natom")
            validate_non_negative(trajs_c.nframe, "trajs_c.nframe")
            if trajs_c.natom > 0 and trajs_c.nframe > 0:
                validate_trajectory_dimensions(trajs_c.nframe, trajs_c.natom)

        if not trajs_c:
            trajs_c = STrajectoriesC()
            LibGenesis().lib.init_empty_s_trajectories_c(
                    ctypes.byref(trajs_c),
                    ctypes.byref(ctypes.c_int(natom)),
                    ctypes.byref(ctypes.c_int(nframe)))
        self.c_obj = trajs_c
        self.natom = trajs_c.natom
        self.nframe = trajs_c.nframe

        # Validate coords pointer before array access
        if trajs_c.nframe > 0 and trajs_c.natom > 0:
            if not trajs_c.coords:
                raise GenesisMemoryError(
                    "trajs_c.coords is NULL but nframe > 0 and natom > 0"
                )
            self.coords = np.ctypeslib.as_array(
                    ctypes.cast(trajs_c.coords, ctypes.POINTER(
                        ctypes.c_double * trajs_c.nframe * trajs_c.natom * 3
                        )).contents
                    ).reshape(trajs_c.nframe, trajs_c.natom, 3)
        else:
            self.coords = np.empty((0, 0, 3), dtype=np.float64)

        # Validate pbc_boxes pointer before array access
        if trajs_c.nframe > 0:
            if not trajs_c.pbc_boxes:
                raise GenesisMemoryError(
                    "trajs_c.pbc_boxes is NULL but nframe > 0"
                )
            self.pbc_boxes = np.ctypeslib.as_array(
                    ctypes.cast(trajs_c.pbc_boxes, ctypes.POINTER(
                        ctypes.c_double * trajs_c.nframe * 3 * 3)).contents
                    ).reshape(trajs_c.nframe, 3, 3)
        else:
            self.pbc_boxes = np.empty((0, 3, 3), dtype=np.float64)

        self._mem_owner = mem_owner

        # Initialize lazy loading fields (non-lazy by default)
        self.is_lazy = False
        self.lazy_dcd_file = None
        self.lazy_trj_type = TRJ_TYPE_COOR_BOX
        self.selection_indices = None

    def __del__(self) -> None:
        try:
            self.free()
        except Exception:
            pass

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.free()

    def __deepcopy__(self, memo) -> Self:
        dst_c = STrajectoriesC()
        LibGenesis().lib.deep_copy_s_trajectories_c(
                ctypes.byref(self.c_obj), ctypes.byref(dst_c))
        return STrajectories(trajs_c=dst_c)

    def free(self):
        """Deallocate resources.

        This method handles cleanup of Fortran-allocated memory through the
        shared library. It attempts graceful cleanup but will not raise
        exceptions to avoid issues during garbage collection.
        """
        mem_owner = bool(getattr(self, "_mem_owner", False))
        cobj = getattr(self, "c_obj", None)

        if (not mem_owner) or (cobj is None):
            return
        if hasattr(cobj, "value") and not cobj.value:
            return

        try:
            lib = LibGenesis().lib
        except (OSError, AttributeError) as e:
            # Library not available (e.g., during interpreter shutdown)
            warnings.warn(
                f"Could not access GENESIS library during cleanup: {e}",
                ResourceWarning,
                stacklevel=2
            )
            self._mem_owner = False
            if hasattr(self, "src_c_obj"):
                self.src_c_obj = ctypes.c_void_p()
            return

        try:
            if isinstance(cobj, STrajectoriesC):
                lib.deallocate_s_trajectories_c.argtypes = [ctypes.POINTER(STrajectoriesC)]
                lib.deallocate_s_trajectories_c.restype = None
                lib.deallocate_s_trajectories_c(ctypes.byref(cobj))
            elif isinstance(cobj, ctypes.POINTER(STrajectoriesC)):
                lib.deallocate_s_trajectories_c.argtypes = [ctypes.POINTER(STrajectoriesC)]
                lib.deallocate_s_trajectories_c.restype = None
                lib.deallocate_s_trajectories_c(cobj)
            elif isinstance(cobj, ctypes.c_void_p):
                lib.deallocate_s_trajectories_c.argtypes = [ctypes.c_void_p]
                lib.deallocate_s_trajectories_c.restype = None
                lib.deallocate_s_trajectories_c(cobj)
            else:
                elty = getattr(cobj, "_type_", None)
                if isinstance(cobj, ctypes._Pointer) and isinstance(elty, type) and issubclass(elty, ctypes.Structure):
                    lib.deallocate_s_trajectories_c.argtypes = [ctypes.POINTER(elty)]
                    lib.deallocate_s_trajectories_c.restype = None
                    lib.deallocate_s_trajectories_c(cobj)
                else:
                    raise TypeError(f"Unsupported handle type for free(): {type(cobj)}")
        except (OSError, ctypes.ArgumentError) as e:
            warnings.warn(
                f"Error during trajectory deallocation: {e}",
                ResourceWarning,
                stacklevel=2
            )
        finally:
            self._mem_owner = False
            if hasattr(self, "src_c_obj"):
                self.src_c_obj = ctypes.c_void_p()
            else:
                self.c_obj = None


    def get_c_obj(self) -> STrajectoriesC:
        return self.c_obj

    def from_trajectories_c(trajs_c: STrajectoriesC, mem_owner=True) -> Self:
        return STrajectories(trajs_c=trajs_c, mem_owner=mem_owner)

    @classmethod
    def from_numpy(cls, coords: npt.NDArray[np.float64],
                   pbc_box: npt.NDArray[np.float64] = None,
                   mem_owner: bool = True) -> Self:
        """Create STrajectories from numpy arrays (zerocopy pattern).

        This creates an STrajectories where the coordinate data is owned by
        numpy arrays rather than Fortran. This is used by zerocopy functions
        where Python pre-allocates memory and Fortran writes directly to it.

        Args:
            coords: Coordinate array of shape (nframe, natom, 3)
            pbc_box: Optional PBC box array of shape (nframe, 3)
                     If None, a zero-filled (nframe, 3, 3) array is used
            mem_owner: If True, this object manages the numpy array lifetime

        Returns:
            STrajectories instance with numpy-owned memory

        Note:
            This differs from the standard constructor which creates
            Fortran-allocated memory. With from_numpy, the arrays are
            managed by Python/numpy, and no Fortran deallocation is needed.
        """
        if coords.ndim != 3 or coords.shape[2] != 3:
            raise GenesisValidationError(
                f"coords must have shape (nframe, natom, 3), got {coords.shape}"
            )

        nframe, natom, _ = coords.shape

        # Create instance without Fortran allocation
        obj = object.__new__(cls)
        obj.natom = natom
        obj.nframe = nframe
        obj.coords = coords
        obj._mem_owner = False  # Don't try to deallocate Fortran memory

        # Handle pbc_box
        if pbc_box is not None:
            if pbc_box.ndim == 2 and pbc_box.shape == (nframe, 3):
                # Expand (nframe, 3) to (nframe, 3, 3) diagonal matrix
                pbc_boxes = np.zeros((nframe, 3, 3), dtype=np.float64)
                for i in range(nframe):
                    np.fill_diagonal(pbc_boxes[i], pbc_box[i])
                obj.pbc_boxes = pbc_boxes
            elif pbc_box.ndim == 3 and pbc_box.shape == (nframe, 3, 3):
                obj.pbc_boxes = pbc_box
            else:
                raise GenesisValidationError(
                    f"pbc_box must have shape (nframe, 3) or (nframe, 3, 3), "
                    f"got {pbc_box.shape}"
                )
        else:
            obj.pbc_boxes = np.zeros((nframe, 3, 3), dtype=np.float64)

        # Keep reference to numpy arrays to prevent garbage collection
        if mem_owner:
            obj._numpy_coords = coords
            obj._numpy_pbc_boxes = obj.pbc_boxes

        # Create C structure that wraps the numpy arrays
        # This allows using the trajectory with Fortran analysis functions
        c_obj = STrajectoriesC()
        c_obj.natom = natom
        c_obj.nframe = nframe
        # Point to numpy array data (no copy)
        c_obj.coords = obj.coords.ctypes.data_as(ctypes.c_void_p).value
        c_obj.pbc_boxes = obj.pbc_boxes.ctypes.data_as(ctypes.c_void_p).value
        obj.c_obj = c_obj

        # Initialize lazy loading fields (non-lazy)
        obj.is_lazy = False
        obj.lazy_dcd_file = None
        obj.lazy_trj_type = TRJ_TYPE_COOR_BOX
        obj.selection_indices = None

        return obj

    @classmethod
    def from_lazy(cls, dcd_file: str, trj_type: int, nframe: int, natom: int,
                  selection_indices: Optional[npt.NDArray[np.int32]] = None) -> Self:
        """Create a lazy STrajectories that reads from DCD file on demand.

        This creates an STrajectories that doesn't load coordinate data into
        memory. Instead, it holds a reference to the DCD file, and analysis
        functions will read frames directly from the file.

        Args:
            dcd_file: Path to DCD trajectory file
            trj_type: Trajectory type (TRJ_TYPE_COOR or TRJ_TYPE_COOR_BOX)
            nframe: Number of frames in DCD file (from DCD header)
            natom: Number of atoms in DCD file (from DCD header)
            selection_indices: Optional pre-computed atom indices (1-indexed).
                              If None, all atoms are used.

        Returns:
            STrajectories instance with is_lazy=True

        Note:
            Lazy trajectories have restrictions:
            - No fitting (requires all coordinates in memory)
            - Single DCD file only (no concatenation)
            - Sequential frame access only
        """
        obj = object.__new__(cls)
        obj.natom = natom
        obj.nframe = nframe
        obj.coords = None  # No data loaded
        obj.pbc_boxes = None  # No data loaded
        obj.c_obj = None  # No C structure needed
        obj._mem_owner = False  # Nothing to deallocate

        # Set lazy loading fields
        obj.is_lazy = True
        obj.lazy_dcd_file = dcd_file
        obj.lazy_trj_type = trj_type
        obj.selection_indices = selection_indices

        return obj

    def save_dcd(self, path: str, molecule: SMolecule) -> None:
        """Write the trajectory to a DCD file via mdtraj.

        The coordinates and box vectors are converted from GENESIS angstroms to
        the nanometers mdtraj expects (the same conversion
        :meth:`to_mdtraj_trajectory` performs). The topology is taken from
        ``molecule`` so the DCD can be paired with e.g. ``molecule``'s PDB when
        deposited or reloaded.

        Args:
            path: Destination DCD path.
            molecule: SMolecule providing the topology (atom count must match
                this trajectory's ``natom``).

        Raises:
            ImportError: If mdtraj is not installed (an optional dependency).
            GenesisValidationError: If the trajectory has no in-memory
                coordinates (e.g. a lazy trajectory) or the atom counts differ.
        """
        if not hasattr(self, "to_mdtraj_trajectory"):
            raise ImportError(
                "save_dcd requires mdtraj, which is an optional dependency. "
                "Install it with 'pip install mdtraj' (or 'pip install "
                "genepie[mdtraj]')."
            )
        if getattr(self, "is_lazy", False) or self.coords is None:
            raise GenesisValidationError(
                "save_dcd requires in-memory coordinates; this trajectory is "
                "lazy. Load it eagerly (crd_convert without lazy=True) first."
            )
        if molecule.num_atoms != self.natom:
            raise GenesisValidationError(
                f"save_dcd: molecule has {molecule.num_atoms} atoms but the "
                f"trajectory has {self.natom}. Pass the matching (subset) "
                "molecule."
            )
        self.to_mdtraj_trajectory(molecule).save_dcd(path)


try:
    import mdtraj as md

    def to_mdtraj_trajectory(self, smol: SMolecule) -> md.Trajectory:
        """

        Returns
            MDTraj Trajectory
        -------
        """
        traj = md.Trajectory(xyz=(self.coords / 10),
                             topology=smol.to_mdtraj_topology())
        traj.unitcell_vectors = np.array(self.pbc_boxes / 10)
        return traj

    STrajectories.to_mdtraj_trajectory = to_mdtraj_trajectory

    @staticmethod
    def from_mdtraj_trajectory(src: md.Trajectory) -> tuple[Self, SMolecule]:
        straj = STrajectories(src.n_atoms, src.n_frames)
        straj.coords[:] = src.xyz * 10
        straj.pbc_boxes[:] = src.unitcell_vectors * 10
        return (straj, SMolecule.from_mdtraj_topology(src.topology))

    STrajectories.from_mdtraj_trajectory = from_mdtraj_trajectory

except ImportError:
    pass


try:
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader
    from MDAnalysis.lib.mdamath import triclinic_vectors
    from MDAnalysis.lib.mdamath import triclinic_box

    def to_mdanalysis_universe(self, smol: SMolecule) -> mda.Universe:
        uni = smol.to_mdanalysis_universe()
        self.add_coordinates_to_mdanalysis_universe(uni)
        return uni

    STrajectories.to_mdanalysis_universe = to_mdanalysis_universe

    def add_coordinates_to_mdanalysis_universe(
            self, uni: mda.Universe) -> None:
        uni.load_new(self.coords, format=MemoryReader, order='fac', dt=1.0)
        for sb, ut in zip(self.pbc_boxes,  uni.trajectory):
            ut.dimensions = triclinic_box(sb[0, :], sb[1, :], sb[2, :])

    STrajectories.add_coordinates_to_mdanalysis_universe \
        = add_coordinates_to_mdanalysis_universe

    @staticmethod
    def from_mdanalysis_universe(src: mda.Universe) -> tuple[Self, SMolecule]:
        straj = STrajectories(natom=src.atoms.n_atoms,
                              nframe=len(src.trajectory))
        for i, ts in enumerate(src.trajectory):
            straj.coords[i, :, :] = src.atoms.positions
            straj.pbc_boxes[i, :, :] = triclinic_vectors(src.dimensions)
        mol = SMolecule.from_mdanalysis_universe(src)
        return (straj, mol)

    STrajectories.from_mdanalysis_universe = from_mdanalysis_universe

except ImportError:
    pass


class STrajectoriesArray:
    def __init__(self, traj_c_array: ctypes.c_void_p, len_array: int):
        # Validate inputs
        validate_non_negative(len_array, "len_array")

        if len_array > 0 and not traj_c_array:
            raise GenesisMemoryError(
                "traj_c_array is NULL but len_array > 0"
            )

        self.src_c_obj = ctypes.cast(
                traj_c_array, ctypes.POINTER(STrajectoriesC))
        self.traj_array = []
        for i in range(len_array):
            self.traj_array.append(
                    STrajectories.from_trajectories_c(
                        self.src_c_obj[i], mem_owner=False))

    def __del__(self) -> None:
        try:
            self.free()
        except Exception:
            pass

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.free()

    def __getitem__(self, index) -> STrajectories:
        return self.traj_array[index]

    def __len__(self) -> int:
        return len(self.traj_array)

    def free(self) -> None:
        """Deallocate resources.

        This method handles cleanup of Fortran-allocated trajectory arrays.
        It attempts graceful cleanup but will not raise exceptions to avoid
        issues during garbage collection.
        """
        cobj = getattr(self, "c_obj", None)
        arr = getattr(self, "traj_array", None)

        # Return early if there's nothing to deallocate
        if cobj is None or not bool(cobj):
            return

        n = len(arr) if isinstance(arr, (list, tuple)) else 0
        len_array = ctypes.c_int(n)

        try:
            lib = LibGenesis().lib
        except (OSError, AttributeError) as e:
            # Library not available (e.g., during interpreter shutdown)
            warnings.warn(
                f"Could not access GENESIS library during cleanup: {e}",
                ResourceWarning,
                stacklevel=2
            )
            self.src_c_obj = None
            if arr is not None:
                arr.clear()
            return

        try:
            lib.deallocate_s_trajectories_c_array(
                ctypes.byref(cobj),
                ctypes.byref(len_array))
        except (OSError, ctypes.ArgumentError) as e:
            warnings.warn(
                f"Error during trajectory array deallocation: {e}",
                ResourceWarning,
                stacklevel=2
            )
        finally:
            self.src_c_obj = None
            if arr is not None:
                arr.clear()

    def get_c_obj(self):
        return self.c_obj

    @property
    def c_obj(self):
        return self.src_c_obj

    def __iter__(self):
        return STrajectoriesArrayIterator(self)

    def join(self) -> STrajectories:
        ret = STrajectoriesC()
        LibGenesis().lib.join_s_trajectories_c(
                ctypes.byref(self.src_c_obj),
                ctypes.byref(ctypes.c_int(len(self))),
                ctypes.byref(ret))
        return STrajectories.from_trajectories_c(ret)

    def from_s_trajectories(*trajs: STrajectories) -> Self:
        if len(trajs) <= 0:
            return STrajectoriesArray(ctypes.c_void_p(), 0)
        natom = trajs[0].natom
        for t in trajs:
            if natom != t.natom:
                raise ValueError("different number of atoms in trajectories")
        c_trajs = ctypes.POINTER(STrajectoriesC)
        LibGenesis().lib.allocate_s_trajectories_c_array(
                ctypes.byref(ctypes.c_int(len(trajs))),
                ctypes.byref(c_trajs))
        for i, t in enumerate(trajs):
            c_trajs[i] = t
        return STrajectoriesArray(c_trajs, len(trajs))


class STrajectoriesArrayIterator:
    def __init__(self, src_array: STrajectoriesArray):
        self.i = 0
        self.src_array = src_array

    def __iter__(self):
        return self

    def __next__(self) -> STrajectories:
        if (self.i < len(self.src_array)):
            self.i += 1
            return self.src_array[self.i - 1]
        else:
            raise StopIteration
