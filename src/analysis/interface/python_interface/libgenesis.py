import ctypes
import os
import threading
from s_molecule_c import SMoleculeC
from s_trajectories_c import STrajectoriesC

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

        self.lib.define_molecule_from_pdb_psf.argtypes = [
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.POINTER(SMoleculeC)]
        self.lib.define_molecule_from_pdb_psf.restype = None

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

        self.lib.test_conv_c2f.argtypes = [
                ctypes.POINTER(SMoleculeC)]
        self.lib.test_conv_c2f.restype = None

        self.lib.crd_convert_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.c_char_p,
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.crd_convert_c.restype = None

        self.lib.trj_analysis_c.argtyes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.trj_analysis_c.restype = None

        self.lib.rg_analysis_c.argtyes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.POINTER(ctypes.c_void_p),
                ]
        self.lib.rg_analysis_c.restype = None

        self.lib.ra_analysis_c.argtyes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.POINTER(ctypes.c_void_p),
                ]
        self.lib.ra_analysis_c.restype = None

        self.lib.dr_analysis_c.argtyes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.POINTER(ctypes.c_void_p),
                ]
        self.lib.dr_analysis_c.restype = None

        self.lib.ma_analysis_c.argtyes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.ma_analysis_c.restype = None

        self.lib.diffusion_analysis_c.argtyes = [
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.POINTER(ctypes.c_void_p),
                ]
        self.lib.diffusion_analysis_c.restype = None

        self.lib.hb_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ]
        self.lib.hb_analysis_c.restype = None

        self.lib.aa_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ]
        self.lib.aa_analysis_c.restype = None

        self.lib.wa_analysis_c.argtyes = [
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ]
        self.lib.wa_analysis_c.restype = None

        self.lib.mbar_analysis_c.argtyes = [
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ]
        self.lib.mbar_analysis_c.restype = None

        self.lib.kc_analysis_c.argtyes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ]
        self.lib.kc_analysis_c.restype = None

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

        
