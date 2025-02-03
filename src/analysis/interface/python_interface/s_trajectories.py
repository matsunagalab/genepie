import ctypes
from typing import Self
import numpy as np
from s_trajectories_c import STrajectoriesC
from libgenesis import LibGenesis


class STrajectories:
    natom: int
    nframe: int

    def __init__(self, natom: int = 0, nframe: int = 0,
                 trajs_c: STrajectoriesC = None):
        if not trajs_c:
            trajs_c = STrajectoriesC()
            LibGenesis().lib.init_empty_s_trajectories_c(
                    ctypes.byref(trajs_c),
                ctypes.byref(ctypes.c_int(natom)),
                ctypes.byref(ctypes.c_int(nframe)))
        self._trajs_c = trajs_c
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
                ctypes.byref(self._trajs_c))

    def from_trajectories_c(trajs_c: STrajectoriesC) -> Self:
        return STrajectories(trajs_c=trajs_c)


class STrajectoriesArray:
    def __init__(self, traj_c_array: ctypes.c_void_p, len_array: int):
        if traj_c_array:
            self._src_c_obj = traj_c_array
            self.traj_c_array = []
            traj_p = ctypes.cast(traj_c_array, ctypes.POINTER(STrajectoriesC))
            for i in range(0, len_array):
                self.traj_c_array.append(STrajectories.from_trajectories_c(traj_p[i]))

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.free()

    def free(self):
        """deallocate resources"""
        len_array = ctypes.c_int(len(self.traj_c_array))
        LibGenesis().lib.deallocate_s_trajectories_c_array(
                ctypes.byref(self._src_c_obj),
                ctypes.byref(len_array))

    def get_array(self):
        return self.traj_c_array
