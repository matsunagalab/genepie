import ctypes
from typing import Self
import numpy as np
import numpy.typing as npt
from s_trajectories_c import STrajectoriesC
from libgenesis import LibGenesis


class STrajectories:
    natom: int
    nframe: int
    coords: npt.NDArray[np.float64]
    pbc_boxes: npt.NDArray[np.float64]
    c_obj: STrajectoriesC

    def __init__(self, natom: int = 0, nframe: int = 0,
                 trajs_c: STrajectoriesC = None, mem_owner = True):
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
            trajs_c.coords, shape=(trajs_c.nframe, trajs_c.natom, 3))
        self.pbc_boxes = np.ctypeslib.as_array(
            trajs_c.pbc_boxes, shape=(trajs_c.nframe, 3, 3))
        self.mem_owner = mem_owner

    def __del__(self) -> None:
        self.free()

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.free()

    def __deepcopy__(self, memo) -> Self:
        dst_c = STrajectoriesC()
        LibGenesis().lib.deep_copy_s_trajectories_c(
                ctypes.byref(self.c_obj), ctypes.byref(dst_c))
        return STrajectories(trajs_c = dst_c)

    def free(self):
        """deallocate resources"""
        if self.mem_owner and self.c_obj:
            LibGenesis().lib.deallocate_s_trajectories_c(
                    ctypes.byref(self.c_obj))
            self.c_obj = None

    def get_c_obj(self) -> STrajectoriesC:
        return self.c_obj

    def from_trajectories_c(trajs_c: STrajectoriesC, mem_owner=True) -> Self:
        return STrajectories(trajs_c=trajs_c, mem_owner=mem_owner)


class STrajectoriesArray:
    def __init__(self, traj_c_array: ctypes.c_void_p, len_array: int):
        self.src_c_obj = ctypes.cast(
                traj_c_array, ctypes.POINTER(STrajectoriesC))
        self.traj_array = []
        for i in range(0, len_array):
            self.traj_array.append(
                    STrajectories.from_trajectories_c(self.src_c_obj[i], mem_owner=False))

    def __del__(self) -> None:
        self.free()

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
        if self.src_c_obj:
            len_array = ctypes.c_int(len(self.traj_array))
            LibGenesis().lib.deallocate_s_trajectories_c_array(
                    ctypes.byref(self.src_c_obj),
                    ctypes.byref(len_array))
            self.src_c_obj = ctypes.POINTER(STrajectoriesC)()
            self.traj_array.clear()
            self.src_c_obj = None


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
