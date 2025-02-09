import ctypes
from typing import Self
import numpy as np
from s_trajectories_c import STrajectoriesC
from libgenesis import LibGenesis


class STrajectories:
    natom: int
    nframe: int
    # coords:
    # pbc_boxes:
    c_obj: STrajectoriesC

    def __init__(self, natom: int = 0, nframe: int = 0,
                 trajs_c: STrajectoriesC = None):
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

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.free()

    def free(self):
        """deallocate resources"""
        LibGenesis().lib.deallocate_s_trajectories_c(
                ctypes.byref(self.c_obj))

    def get_c_obj(self) -> STrajectoriesC:
        return self.c_obj

    def from_trajectories_c(trajs_c: STrajectoriesC) -> Self:
        return STrajectories(trajs_c=trajs_c)


class STrajectoriesArray:
    def __init__(self, traj_c_array: ctypes.c_void_p, len_array: int):
        if traj_c_array:
            self.src_c_obj = traj_c_array
            self.traj_array = []
            self.traj_p = ctypes.cast(traj_c_array, ctypes.POINTER(STrajectoriesC))
            for i in range(0, len_array):
                self.traj_array.append(
                        STrajectories.from_trajectories_c(self.traj_p[i]))

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.free()

    def __getitem__(self, index) -> STrajectories:
        return self.traj_array[index]

    def __len__(self) -> int:
        return len(self.traj_array)

    def free(self):
        """deallocate resources"""
        len_array = ctypes.c_int(len(self.traj_array))
        LibGenesis().lib.deallocate_s_trajectories_c_array(
                ctypes.byref(self.src_c_obj),
                ctypes.byref(len_array))


    def __iter__(self):
        return STrajectoriesArrayIterator(self)

    def get_array(self):
        return self.traj_array

    def join(self) -> STrajectories:
        ret = STrajectoriesC()
        LibGenesis().lib.join_s_trajectories_c(
                ctypes.byref(self.src_c_obj),
                ctypes.byref(ctypes.c_int(len(self))),
                ctypes.byref(ret))
        return STrajectories.from_trajectories_c(ret)


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
