import ctypes


class STrajectoriesC(ctypes.Structure):
     _fields_ = [("coords", ctypes.c_void_p),
                 ("pbc_boxes", ctypes.c_void_p),
                 ("nframe", ctypes.c_int),
                 ("natom", ctypes.c_int)]
