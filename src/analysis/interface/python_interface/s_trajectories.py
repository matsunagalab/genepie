import ctypes
from typing import Self
import numpy as np
from libgenesis import LibGenesis


class STrajectoriesC(ctypes.Structure):
     _fields_ = [("coords", ctypes.c_ptr),
                 ("pbc_boxes", ctypes.c_ptr),
                 ("nframe", ctypes.c_int),
                 ("natom", ctypes.c_int)]


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
