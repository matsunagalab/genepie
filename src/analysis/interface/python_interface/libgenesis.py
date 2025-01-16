import ctypes
import os
import threading
from s_molecule import SMoleculeC

class LibGenesis:
    """singleton
    """

    _instance = None
    _lock = threading.Lock()

    def __new__(cls):
        with cls._lock:
            if cls._instance is None:
                cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        lib_name = 'libpython_interface.so'
        lib_dir = os.path.join(os.path.dirname(__file__), '.libs')
        lib_path = os.path.join(lib_dir, lib_name)
        if not os.path.exists(lib_path):
            raise FileNotFoundError(
                    f"Library file {lib_name} not found in {lib_dir}")
        self.lib = ctypes.CDLL(lib_path)
        self.lib.define_molecule_from_pdb.argtypes = [
                ctypes.c_char_p,
                ctypes.POINTER(SMoleculeC)]
        self.lib.define_molecule_from_pdb.restype = None

        self.lib.deallocate_s_molecule_c.argtypes = [
                ctypes.POINTER(SMoleculeC)]
        self.lib.deallocate_s_molecule_c.restype = None

        self.lib.allocate_s_molecule_c.argtypes = [
                ctypes.POINTER(SMoleculeC)]
        self.lib.allocate_s_molecule_c.restype = None

        self.lib.test_conv_c2f.argtypes = [
                ctypes.POINTER(SMoleculeC)]
        self.lib.test_conv_c2f.restype = None
