import ctypes
from typing import Self
import numpy as np
import numpy.typing as npt
from .s_trajectories_c import STrajectoriesC
from .libgenesis import LibGenesis
from .s_molecule import SMolecule


class STrajectories:
    natom: int
    nframe: int
    coords: npt.NDArray[np.float64]  # shape=(n_frame, n_atom, 3)
    pbc_boxes: npt.NDArray[np.float64]  # shape=(trajs_c.nframe, 3, 3)
    c_obj: STrajectoriesC

    def __init__(self, natom: int = 0, nframe: int = 0,
                 trajs_c: STrajectoriesC = None, mem_owner: bool = True):
        if not trajs_c:
            trajs_c = STrajectoriesC()
            LibGenesis().lib.init_empty_s_trajectories_c(
                    ctypes.byref(trajs_c),
                    ctypes.byref(ctypes.c_int(natom)),
                    ctypes.byref(ctypes.c_int(nframe)))
        self.c_obj = trajs_c
        self.natom = trajs_c.natom
        self.nframe = trajs_c.nframe
        self.coords = np.ctypeslib.as_array(
                ctypes.cast(trajs_c.coords, ctypes.POINTER(
                    ctypes.c_double * trajs_c.nframe * trajs_c.natom * 3
                    )).contents
                ).reshape(trajs_c.nframe, trajs_c.natom, 3)
        self.pbc_boxes = np.ctypeslib.as_array(
                ctypes.cast(trajs_c.pbc_boxes, ctypes.POINTER(
                    ctypes.c_double * trajs_c.nframe * 3 * 3)).contents
                ).reshape(trajs_c.nframe, 3, 3)
        self._mem_owner = mem_owner

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
        """deallocate resources"""
        mem_owner = bool(getattr(self, "_mem_owner", False))
        cobj = getattr(self, "c_obj", None)

        if (not mem_owner) or (cobj is None):
            return
        if hasattr(cobj, "value") and not cobj.value:
            return

        try:
          lib = LibGenesis().lib
        except Exception:
            self._mem_owner = False

            if hasattr(self, "src_c_obj"):
                self.src_c_obj = ctypes.c_void_p()
            else:
                try: 
                    self.c_obj = None
                except Exception:
                    pass
        try:
            STrajectoriesC
        except NameError:
            pass

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
                if isinstance(cobj, ctypes._Pointer) and isinstance(elty,  type) and issubclass(elty, ctypes.Structure):
                    lib.deallocate_s_trajectories_c.argtypes = [ctypes.POINTER(elty)]
                    lib.deallocate_s_trajectories_c.restype = None
                    lib.deallocate_s_trajectories_c(cobj)
                else:
                    raise TypeError(f"Unsupported handle type for free(): {type(cobj)}")


        finally:
            self._mem_owner = False
            try:
                if hasattr(self, "src_c_obj"):
                    self.src_c_obj = ctypes.c_void_p()
                else:
                    self.c_obj = None
            except Exception:
                pass


    def get_c_obj(self) -> STrajectoriesC:
        return self.c_obj

    def from_trajectories_c(trajs_c: STrajectoriesC, mem_owner=True) -> Self:
        return STrajectories(trajs_c=trajs_c, mem_owner=mem_owner)


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
        self.src_c_obj = ctypes.cast(
                traj_c_array, ctypes.POINTER(STrajectoriesC))
        self.traj_array = []
        for i in range(0, len_array):
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
        """deallocate resources"""
        try:
            cobj = getattr(self, "c_obj", None)
            arr = getattr(self, "traj_array", None)

            if not cobj is None or not bool(cobj):
                return

            n = len(arr) if isinstance(arr, (list, tuple)) else 0
            len_array = types.c_int(n)

            try:
                lib = LibGenesis().lib
            except Exception:
                self.src_c_obj = None
                if hasattr(self, "traj_array") and arr is not None:
                    try: arr.clear()
                    except Exception: pass
                return

            try:
                lib.deallocate_s_trajectories_c_array(
                    ctypes.byref(cobj),
                    ctypes.byref(len_array))
            finally:
                self.src_c_obj = None
                if hasattr(self, "traj_array") and arr is not None:
                    try: arr.clear()
                    except Exception: pass
        except Exception:
            try:
                self.src_c_obj = None
            except Exception:
                pass

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
