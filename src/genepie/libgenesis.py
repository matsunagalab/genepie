import ctypes
import os
import threading
from .libloader import load_genesis_lib
from .s_molecule_c import SMoleculeC
from .s_trajectories_c import STrajectoriesC

class LibGenesis:
    """Thread-safe singleton for GENESIS library access.

    Note: While the singleton pattern is thread-safe, the underlying
    Fortran library uses module-level save pointers and is NOT thread-safe.
    All GENESIS function calls should be made from a single thread.
    """

    _instance = None
    _lock = threading.Lock()

    def __new__(cls):
        # Double-check locking pattern for thread safety
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    instance = super().__new__(cls)
                    instance.lib = load_genesis_lib()
                    instance._initialized = True
                    instance._setup_function_signatures()
                    cls._instance = instance
        return cls._instance

    def __init__(self):
        # All initialization done in __new__ to avoid race conditions
        pass

    def _setup_function_signatures(self):

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

        # Atom selection using GENESIS selection syntax
        self.lib.selection_c.argtypes = [
                ctypes.POINTER(SMoleculeC),     # molecule_c
                ctypes.c_char_p,                # selection_str
                ctypes.c_int,                   # selection_len
                ctypes.POINTER(ctypes.c_void_p),  # indices (output)
                ctypes.POINTER(ctypes.c_int),   # n_indices (output)
                ctypes.POINTER(ctypes.c_int),   # status (output)
                ctypes.c_char_p,                # msg (output)
                ctypes.c_int]                   # msglen
        self.lib.selection_c.restype = None

        self.lib.deallocate_selection_c.argtypes = []
        self.lib.deallocate_selection_c.restype = None

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

        # trj_analysis_zerocopy_c(s_trajes_c, ana_period,
        #                         dist_list_ptr, n_dist,
        #                         angl_list_ptr, n_angl,
        #                         tors_list_ptr, n_tors,
        #                         result_distance, result_angle, result_torsion,
        #                         n_frames, status, msg, msglen)
        self.lib.trj_analysis_zerocopy_c.argtypes = [
                ctypes.POINTER(STrajectoriesC),   # s_trajes_c
                ctypes.c_int,                     # ana_period
                ctypes.c_void_p,                  # dist_list_ptr
                ctypes.c_int,                     # n_dist
                ctypes.c_void_p,                  # angl_list_ptr
                ctypes.c_int,                     # n_angl
                ctypes.c_void_p,                  # tors_list_ptr
                ctypes.c_int,                     # n_tors
                ctypes.POINTER(ctypes.c_void_p),  # result_distance
                ctypes.POINTER(ctypes.c_void_p),  # result_angle
                ctypes.POINTER(ctypes.c_void_p),  # result_torsion
                ctypes.POINTER(ctypes.c_int),     # n_frames
                ctypes.POINTER(ctypes.c_int),     # status
                ctypes.c_char_p,                  # msg
                ctypes.c_int,                     # msglen
                ]
        self.lib.trj_analysis_zerocopy_c.restype = None

        # trj_analysis_zerocopy_com_c - trajectory analysis with COM (zerocopy)
        self.lib.trj_analysis_zerocopy_com_c.argtypes = [
                ctypes.c_void_p,                  # mass_ptr
                ctypes.c_int,                     # n_atoms
                ctypes.POINTER(STrajectoriesC),   # s_trajes_c
                ctypes.c_int,                     # ana_period
                # Atom-based measurements
                ctypes.c_void_p,                  # dist_list_ptr
                ctypes.c_int,                     # n_dist
                ctypes.c_void_p,                  # angl_list_ptr
                ctypes.c_int,                     # n_angl
                ctypes.c_void_p,                  # tors_list_ptr
                ctypes.c_int,                     # n_tors
                # COM distance
                ctypes.c_void_p,                  # cdis_atoms_ptr
                ctypes.c_int,                     # n_cdis_atoms
                ctypes.c_void_p,                  # cdis_offsets_ptr
                ctypes.c_int,                     # n_cdis_offsets
                ctypes.c_void_p,                  # cdis_pairs_ptr
                ctypes.c_int,                     # n_cdis
                # COM angle
                ctypes.c_void_p,                  # cang_atoms_ptr
                ctypes.c_int,                     # n_cang_atoms
                ctypes.c_void_p,                  # cang_offsets_ptr
                ctypes.c_int,                     # n_cang_offsets
                ctypes.c_void_p,                  # cang_triplets_ptr
                ctypes.c_int,                     # n_cang
                # COM torsion
                ctypes.c_void_p,                  # ctor_atoms_ptr
                ctypes.c_int,                     # n_ctor_atoms
                ctypes.c_void_p,                  # ctor_offsets_ptr
                ctypes.c_int,                     # n_ctor_offsets
                ctypes.c_void_p,                  # ctor_quads_ptr
                ctypes.c_int,                     # n_ctor
                # Output pointers
                ctypes.POINTER(ctypes.c_void_p),  # result_distance
                ctypes.POINTER(ctypes.c_void_p),  # result_angle
                ctypes.POINTER(ctypes.c_void_p),  # result_torsion
                ctypes.POINTER(ctypes.c_void_p),  # result_cdis
                ctypes.POINTER(ctypes.c_void_p),  # result_cang
                ctypes.POINTER(ctypes.c_void_p),  # result_ctor
                ctypes.POINTER(ctypes.c_int),     # n_frames
                ctypes.POINTER(ctypes.c_int),     # status
                ctypes.c_char_p,                  # msg
                ctypes.c_int,                     # msglen
                ]
        self.lib.trj_analysis_zerocopy_com_c.restype = None

        self.lib.rg_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ]
        self.lib.rg_analysis_c.restype = None

        self.lib.deallocate_rg_results_c.argtypes = []
        self.lib.deallocate_rg_results_c.restype = None

        # RG analysis with true zero-copy (mass array passed directly)
        self.lib.rg_analysis_zerocopy_c.argtypes = [
                ctypes.c_void_p,                  # mass_ptr (pointer to NumPy array)
                ctypes.c_int,                     # n_atoms
                ctypes.POINTER(STrajectoriesC),   # s_trajes_c
                ctypes.c_int,                     # ana_period
                ctypes.c_void_p,                  # analysis_idx (pointer to int array)
                ctypes.c_int,                     # n_analysis
                ctypes.c_int,                     # mass_weighted (0 or 1)
                ctypes.POINTER(ctypes.c_void_p),  # result_rg (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.rg_analysis_zerocopy_c.restype = None

        # RG analysis with full zero-copy (pre-allocated result array)
        self.lib.rg_analysis_zerocopy_full_c.argtypes = [
                ctypes.c_void_p,                  # mass_ptr (pointer to NumPy array)
                ctypes.c_int,                     # n_atoms
                ctypes.POINTER(STrajectoriesC),   # s_trajes_c
                ctypes.c_int,                     # ana_period
                ctypes.c_void_p,                  # analysis_idx (pointer to int array)
                ctypes.c_int,                     # n_analysis
                ctypes.c_int,                     # mass_weighted (0 or 1)
                ctypes.c_void_p,                  # result_ptr (pre-allocated result)
                ctypes.c_int,                     # result_size
                ctypes.POINTER(ctypes.c_int),     # nstru_out (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.rg_analysis_zerocopy_full_c.restype = None

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

        # RMSD zerocopy (no fitting, arrays passed directly)
        self.lib.rmsd_analysis_zerocopy_c.argtypes = [
                ctypes.c_void_p,                  # mass_ptr
                ctypes.c_void_p,                  # ref_coord_ptr
                ctypes.c_int,                     # n_atoms
                ctypes.POINTER(STrajectoriesC),   # s_trajes_c
                ctypes.c_int,                     # ana_period
                ctypes.c_void_p,                  # analysis_idx
                ctypes.c_int,                     # n_analysis
                ctypes.c_int,                     # mass_weighted
                ctypes.POINTER(ctypes.c_void_p),  # result_rmsd (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.rmsd_analysis_zerocopy_c.restype = None

        # RMSD zerocopy with fitting (structural alignment before RMSD)
        self.lib.rmsd_analysis_zerocopy_fitting_c.argtypes = [
                ctypes.c_void_p,                  # mass_ptr
                ctypes.c_void_p,                  # ref_coord_ptr
                ctypes.c_int,                     # n_atoms
                ctypes.POINTER(STrajectoriesC),   # s_trajes_c
                ctypes.c_int,                     # ana_period
                ctypes.c_void_p,                  # fitting_idx_ptr
                ctypes.c_int,                     # n_fitting
                ctypes.c_void_p,                  # analysis_idx_ptr
                ctypes.c_int,                     # n_analysis
                ctypes.c_int,                     # fitting_method
                ctypes.c_int,                     # mass_weighted
                ctypes.POINTER(ctypes.c_void_p),  # result_rmsd (output)
                ctypes.POINTER(ctypes.c_int),     # nframes_out (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.rmsd_analysis_zerocopy_fitting_c.restype = None

        self.lib.deallocate_rmsd_results_c.argtypes = []
        self.lib.deallocate_rmsd_results_c.restype = None

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

        # DRMS zerocopy (contact data passed directly)
        self.lib.drms_analysis_zerocopy_c.argtypes = [
                ctypes.c_void_p,                  # contact_list_ptr
                ctypes.c_void_p,                  # contact_dist_ptr
                ctypes.c_int,                     # n_contact
                ctypes.POINTER(STrajectoriesC),   # s_trajes_c
                ctypes.c_int,                     # ana_period
                ctypes.c_int,                     # pbc_correct
                ctypes.POINTER(ctypes.c_void_p),  # result_dr (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.drms_analysis_zerocopy_c.restype = None

        self.lib.deallocate_drms_results_c.argtypes = []
        self.lib.deallocate_drms_results_c.restype = None

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

        # diffusion_analysis_zerocopy_c(msd_ptr, ndata, ncols,
        #                               time_step, distance_unit, ndofs,
        #                               start_step, stop_step,
        #                               out_data_ptr, diff_coeff_ptr,
        #                               status, msg, msglen)
        self.lib.diffusion_analysis_zerocopy_c.argtypes = [
                ctypes.c_void_p,         # msd_ptr
                ctypes.c_int,            # ndata
                ctypes.c_int,            # ncols
                ctypes.c_double,         # time_step
                ctypes.c_double,         # distance_unit
                ctypes.c_int,            # ndofs
                ctypes.c_int,            # start_step
                ctypes.c_int,            # stop_step
                ctypes.POINTER(ctypes.c_void_p),  # out_data_ptr
                ctypes.POINTER(ctypes.c_void_p),  # diff_coeff_ptr
                ctypes.POINTER(ctypes.c_int),     # status
                ctypes.c_char_p,                  # msg
                ctypes.c_int,                     # msglen
                ]
        self.lib.diffusion_analysis_zerocopy_c.restype = None

        self.lib.deallocate_diffusion_results_c.argtypes = []
        self.lib.deallocate_diffusion_results_c.restype = None

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

        self.lib.deallocate_hb_results_c.argtypes = []
        self.lib.deallocate_hb_results_c.restype = None

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

        # ATDYN MD/Minimization functions
        self.lib.atdyn_md_c.argtypes = [
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),  # result_energies
                ctypes.POINTER(ctypes.c_int),     # result_nframes
                ctypes.POINTER(ctypes.c_int),     # result_nterms
                ctypes.POINTER(ctypes.c_void_p),  # result_final_coords
                ctypes.POINTER(ctypes.c_int),     # result_natom
                ctypes.POINTER(ctypes.c_int),     # status
                ctypes.c_char_p,                  # msg
                ctypes.c_int,                     # msglen
                ]
        self.lib.atdyn_md_c.restype = None

        self.lib.atdyn_min_c.argtypes = [
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),  # result_energies
                ctypes.POINTER(ctypes.c_int),     # result_nsteps
                ctypes.POINTER(ctypes.c_int),     # result_nterms
                ctypes.POINTER(ctypes.c_void_p),  # result_final_coords
                ctypes.POINTER(ctypes.c_int),     # result_natom
                ctypes.POINTER(ctypes.c_int),     # result_converged
                ctypes.POINTER(ctypes.c_double),  # result_final_gradient
                ctypes.POINTER(ctypes.c_int),     # status
                ctypes.c_char_p,                  # msg
                ctypes.c_int,                     # msglen
                ]
        self.lib.atdyn_min_c.restype = None

        self.lib.deallocate_atdyn_results_c.argtypes = []
        self.lib.deallocate_atdyn_results_c.restype = None

        # Reset atdyn global state for multiple sequential runs
        self.lib.reset_atdyn_state_c.argtypes = []
        self.lib.reset_atdyn_state_c.restype = None
