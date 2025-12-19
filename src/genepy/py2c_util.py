import ctypes
import os
import numpy as np
import numpy.typing as npt
from typing import Union


def write_bool_ndarray(src: npt.NDArray[np.bool_], dst: ctypes.c_void_p):
    dst_p = ctypes.cast(dst, ctypes.POINTER(ctypes.c_bool))
    for i, v in enumerate(np.ravel(src)):
        dst_p[i] = ctypes.c_bool(v)


def write_int_ndarray(src: npt.NDArray[np.int64], dst: ctypes.c_void_p):
    dst_p = ctypes.cast(dst, ctypes.POINTER(ctypes.c_int))
    for i, v in enumerate(np.ravel(src)):
        dst_p[i] = ctypes.c_int(v)


def write_double_ndarray(src: npt.NDArray[np.float64], dst: ctypes.c_void_p):
    dst_p = ctypes.cast(dst, ctypes.POINTER(ctypes.c_double))
    for i, v in enumerate(np.ravel(src)):
        dst_p[i] = ctypes.c_double(v)


def write_fixed_length_string_ndarray(
        src: npt.NDArray[np.str_], dst: ctypes.c_void_p) -> None:
    dst_p = ctypes.cast(dst, ctypes.POINTER(ctypes.c_char))
    for i, v in enumerate(np.ravel(src)):
        dst_p[i] = v.encode('ascii')


def write_pystring_ndarray(src: npt.NDArray[np.object_], dst: ctypes.c_void_p,
                           str_size: int = -1) -> None:
    dst_p = ctypes.cast(dst, ctypes.POINTER(ctypes.c_char))
    if (str_size < 0) and (len(src) > 0):
        str_size = len(src[0])
    k = 0
    for i, v in enumerate(np.ravel(src)):
        for j in range(0, str_size):
            if j <= len(v):
                dst_p[k] = v[j].encode('ascii')
            else:
                dst_p[k] = ' '.encode('ascii')
            k = k + 1


def pathlike_to_byte(path: Union[str, bytes, os.PathLike]) -> bytes:
    if (type(path) is str):
        return path.encode()
    elif (type(path) is bytes):
        return path
    else:
        return os.fspath(path).encode()

def pathlike_to_c_char_p(path: Union[str,bytes,os.PathLike]) -> ctypes.c_char_p:
    if not path or str(path).strip() == '':
        return ctypes.c_char_p(None)  # ← NULL pointer を渡す
    if isinstance(path, bytes):
        return ctypes.c_char_p(path)
    return ctypes.c_char_p(os.fspath(path).encode())
