import ctypes


def conv_bool_array(src: ctypes.c_void_p, size: int):
    if src:
        ptr = ctypes.cast(src, ctypes.POINTER(ctypes.c_bool))
        return list(bool(ptr[i]) for i in range(0, size))
    else:
        return ()


def conv_int_array(src: ctypes.c_void_p, size: int):
    if src:
        ptr = ctypes.cast(src, ctypes.POINTER(ctypes.c_int))
        return list(int(ptr[i]) for i in range(0, size))
    else:
        return ()


def conv_int_array_2d(src: ctypes.c_void_p, size1: int, size2: int):
    if src:
        ptr = ctypes.cast(src, ctypes.POINTER(ctypes.c_int))
        return list(tuple(int(ptr[i * size2 + j]) for j in range(0, size2)) for i in range(0, size1))
    else:
        return ()


def conv_int_array_3d(src: ctypes.c_void_p, size1: int, size2: int, size3: int):
    if src:
        ptr = ctypes.cast(src, ctypes.POINTER(ctypes.c_int))
        return list(
                tuple(
                    tuple(
                        int(ptr[i1 * size2 * size3 + i2 * size3 + i3])
                        for i3 in range(0, size3))
                    for i2 in range(0, size2))
                for i1 in range(0, size1))
    else:
        return ()


def conv_int_array_4d(src: ctypes.c_void_p,
                      size1: int, size2: int, size3: int, size4: int):
    if src:
        ptr = ctypes.cast(src, ctypes.POINTER(ctypes.c_int))
        return list(
                tuple(
                    tuple(
                        tuple(
                            int(ptr[i1 * size2 * size3 * size4
                                    + i2 * size3 * size4 + i3 * size4 + i4])
                            for i4 in range(0, size3))
                        for i3 in range(0, size3))
                    for i2 in range(0, size2))
                for i1 in range(0, size1))
    else:
        return ()


def conv_double_array(src: ctypes.c_void_p, size: int):
    if src:
        ptr = ctypes.cast(src, ctypes.POINTER(ctypes.c_double))
        return list(float(ptr[i]) for i in range(0, size))
    else:
        return ()


def conv_double_array_2d(src: ctypes.c_void_p, size1: int, size2: int):
    if src:
        ptr = ctypes.cast(src, ctypes.POINTER(ctypes.c_double))
        return list(tuple(float(ptr[i * size2 + j]) for j in range(0, size2)) for i in range(0, size1))
    else:
        return ()


def conv_str_array(src: ctypes.c_void_p, size_array: int, size_str: int):
    if src:
        ptr = ctypes.cast(src, ctypes.POINTER(ctypes.c_char))
        return list(conv_str(ptr[i * size_str], size_str) for i in range(0, size_array))
    else:
        return ()


def conv_str(src: ctypes.c_void_p, size: int):
    if src:
        ptr = ctypes.cast(src, ctypes.POINTER(ctypes.c_char))
        return ''.join(str(ptr[i].decode('utf8')) for i in range(0, size))
    return ' ' * size
