import ctypes
import numpy as np
import numpy.typing as npt


def conv_ndarray(src: ctypes.c_void_p, size: int | list[int],
                 type_clang: type, type_numpy: type) -> npt.NDArray:
    if src:
        if type(int) is int:
            size = tuple(size)
        all_size = np.prod(size)
        ptr = ctypes.cast(src, ctypes.POINTER(type_clang))
        return (np.fromiter((ptr[i] for i in range(0, all_size)),
                            dtype=type_numpy, count=all_size)
                .reshape(size))
    else:
        return np.empty(0, dtype=type_numpy)


def conv_bool_ndarray(src: ctypes.c_void_p, size: int | list[int]) \
        -> npt.NDArray[np.bool_]:
    return conv_ndarray(src, size, ctypes.c_bool, np.bool_)


def conv_int_ndarray(src: ctypes.c_void_p, size: int | list[int]) \
        -> npt.NDArray[np.int64]:
    return conv_ndarray(src, size, ctypes.c_int, np.int64)


def conv_double_ndarray(src: ctypes.c_void_p, size: int | list[int]) \
        -> npt.NDArray[np.float64]:
    return conv_ndarray(src, size, ctypes.c_double, np.float64)


def conv_fixed_length_string_ndarray(
        str_array: ctypes.c_void_p, size: int | list[int]) \
                -> npt.NDArray[np.str_]:
    if str_array:
        if type(int) is int:
            size = tuple(size)
        all_size = np.prod(size)
        ptr = ctypes.cast(str_array, ctypes.POINTER(ctypes.c_char))
        dst = np.empty(all_size, dtype=str)
        for i in range(0, all_size):
            dst[i] = ptr[i].decode('ascii')
        return dst.reshape(size)
    else:
        return np.empty(0, dtype=np.str_)


def conv_pystring_ndarray(
        str_array: ctypes.c_void_p, size: tuple[int, int]) \
                -> npt.NDArray[np.object_]:
    if str_array:
        ptr = ctypes.cast(str_array, ctypes.POINTER(ctypes.c_char))
        dst = np.empty(size[0], dtype=np.object_)
        for i in range(0, size[0]):
            cur_ptr = ctypes.cast(ctypes.addressof(ptr.contents) + i * size[1],
                                  ctypes.POINTER(ctypes.c_char))
            dst[i] = conv_string(cur_ptr, size[1])
        return dst
    else:
        return np.empty(0, dtype=np.object_)


def conv_string(c_char_ptr: ctypes.c_void_p, size: int = -1,
                encoding: str = 'utf-8') -> str:
    if not c_char_ptr:
        raise ValueError("Invalid pointer: c_char_ptr is None")
    byte_str = ctypes.string_at(c_char_ptr, size)
    return byte_str.decode(encoding)
