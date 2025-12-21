import ctypes
import os
import threading
from .libloader import load_genesis_lib
from .s_molecule_c import SMoleculeC
from .s_trajectories_c import STrajectoriesC

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
        if getattr(self, "_initialized", False):
            return
        self.lib = load_genesis_lib()
        self._initialized = True

        self.lib.define_molecule_from_file.argtypes = [
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.POINTER(SMoleculeC),
                ]
        self.lib.define_molecule_from_file.restype = ctypes.c_int

        self.lib.deallocate_s_molecule_c.argtypes = [
                ctypes.POINTER(SMoleculeC)]
        self.lib.deallocate_s_molecule_c.restype = None

        self.lib.allocate_s_molecule_c.argtypes = [
                ctypes.POINTER(SMoleculeC)]
        self.lib.allocate_s_molecule_c.restype = None

        self.lib.deallocate_s_trajectories_c.argtypes = [
                ctypes.POINTER(STrajectoriesC)]
        self.lib.deallocate_s_trajectories_c.restype = None

        self.lib.allocate_s_trajectories_c_array.argtypes = [
                ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
        self.lib.allocate_s_trajectories_c_array.restype = None

        self.lib.deallocate_s_trajectories_c_array.argtypes = [
                ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
        self.lib.deallocate_s_trajectories_c_array.restype = None

        self.lib.crd_convert_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.crd_convert_c.restype = None

        self.lib.trj_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.trj_analysis_c.restype = None

        self.lib.rg_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ]
        self.lib.rg_analysis_c.restype = None

        self.lib.ra_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.ra_analysis_c.restype = None

        self.lib.dr_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.dr_analysis_c.restype = None

        self.lib.ma_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.ma_analysis_c.restype = None

        self.lib.diffusion_analysis_c.argtypes = [
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.diffusion_analysis_c.restype = None

        self.lib.hb_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.hb_analysis_c.restype = None

        self.lib.aa_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.aa_analysis_c.restype = None

        self.lib.wa_analysis_c.argtypes = [
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.wa_analysis_c.restype = None

        self.lib.mbar_analysis_c.argtypes = [
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.mbar_analysis_c.restype = None

        self.lib.kc_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.kc_analysis_c.restype = None

        self.lib.export_pdb_to_string_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.c_void_p,
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.export_pdb_to_string_c.restype = None

        self.lib.allocate_c_int_array.argtypes = [
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.allocate_c_int_array.restype = ctypes.c_void_p

        self.lib.allocate_c_double_array.argtypes = [
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.allocate_c_double_array.restype = ctypes.c_void_p

        self.lib.allocate_c_double_array2.argtypes = [
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.allocate_c_double_array2.restype = ctypes.c_void_p

        self.lib.deallocate_int.argtypes = [
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.deallocate_int.restype = None

        self.lib.deallocate_double.argtypes = [
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.deallocate_double.restype = None

        self.lib.deallocate_double2.argtypes = [
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.deallocate_double2.restype = None

        self.lib.deallocate_c_string.argtypes = [
                ctypes.POINTER(ctypes.c_void_p),
                ]
        self.lib.deallocate_c_string.restype = None

        self.lib.join_s_trajectories_c.argtypes = [
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_void_p),
                ]
        self.lib.join_s_trajectories_c.restype = None

        self.lib.deep_copy_s_trajectories_c.argtypes = [
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_void_p),
                ]
        self.lib.deep_copy_s_trajectories_c.restype = None
