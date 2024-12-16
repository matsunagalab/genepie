import ctypes
import os


class s_molecule_c(ctypes.Structure):
     _fields_ = [("num_deg_freedom", ctypes.c_int),
                 ("num_atoms", ctypes.c_int),
                 ("num_bonds", ctypes.c_int),
                 ("num_enm_bonds", ctypes.c_int),
                 ("num_angles", ctypes.c_int),
                 ("num_dihedrals", ctypes.c_int),
                 ("num_impropers", ctypes.c_int),
                 ("num_cmaps", ctypes.c_int),
                 ("num_residues", ctypes.c_int),
                 ("num_molecules", ctypes.c_int),
                 ("num_segments", ctypes.c_int),
                 ("shift_origin", ctypes.c_bool),
                 ("special_hydrogen", ctypes.c_bool),
                 ("total_charge", ctypes.c_double),
                 ("atom_no", ctypes.c_void_p),
                 ("segment_name", ctypes.c_void_p),
                 ("segment_no", ctypes.c_void_p),
                 ("residue_no", ctypes.c_void_p),
                 ("residue_c_no", ctypes.c_void_p),
                 ("residue_name", ctypes.c_void_p),
                 ("atom_name", ctypes.c_void_p),
                 ("atom_cls_name", ctypes.c_void_p),
                 ("atom_cls_no", ctypes.c_void_p),
                 ("charge", ctypes.c_void_p),
                 ("mass", ctypes.c_void_p),
                 ("inv_mass", ctypes.c_void_p),
                 ("imove", ctypes.c_void_p),
                 ("stokes_radius", ctypes.c_void_p),
                 ("inv_stokes_radius", ctypes.c_void_p),
                 ("chain_id", ctypes.c_void_p),
                 ("atom_coord", ctypes.c_void_p),
                 ("atom_occupancy", ctypes.c_void_p),
                 ("atom_temp_factor", ctypes.c_void_p),
                 ("atom_velocity", ctypes.c_void_p),
                 ("light_atom_name", ctypes.c_void_p),
                 ("light_atom_mass", ctypes.c_void_p),
                 ("molecule_no", ctypes.c_void_p),
                 ("bond_list", ctypes.c_void_p),
                 ("enm_list", ctypes.c_void_p),
                 ("angl_list", ctypes.c_void_p),
                 ("dihe_list", ctypes.c_void_p),
                 ("impr_list", ctypes.c_void_p),
                 ("cmap_list", ctypes.c_void_p),
                 ("molecule_atom_no", ctypes.c_void_p),
                 ("molecule_mass", ctypes.c_void_p),
                 ("molecule_name", ctypes.c_void_p),
                 ("atom_refcoord", ctypes.c_void_p),
                 ("atom_fitcoord", ctypes.c_void_p),
                 ("num_pc_modes", ctypes.c_int),
                 ("pc_mode", ctypes.c_void_p),
                 ("fep_topology", ctypes.c_int),
                 ("num_hbonds_singleA", ctypes.c_int),
                 ("num_hbonds_singleB", ctypes.c_int),
                 ("num_atoms_fep", ctypes.c_void_p),
                 ("num_bonds_fep", ctypes.c_void_p),
                 ("num_angles_fep", ctypes.c_void_p),
                 ("num_dihedrals_fep", ctypes.c_void_p),
                 ("num_impropers_fep", ctypes.c_void_p),
                 ("num_cmaps_fep", ctypes.c_void_p),
                 ("bond_list_fep", ctypes.c_void_p),
                 ("angl_list_fep", ctypes.c_void_p),
                 ("dihe_list_fep", ctypes.c_void_p),
                 ("impr_list_fep", ctypes.c_void_p),
                 ("cmap_list_fep", ctypes.c_void_p),
                 ("id_singleA", ctypes.c_void_p),
                 ("id_singleB", ctypes.c_void_p),
                 ("fepgrp", ctypes.c_void_p),
                 ("fepgrp_bond", ctypes.c_void_p),
                 ("fepgrp_angl", ctypes.c_void_p),
                 ("fepgrp_dihe", ctypes.c_void_p),
                 ("fepgrp_cmap", ctypes.c_void_p),
                 ]

def define_prototypes(lib: ctypes.CDLL):
    # 関数のプロトタイプを定義
    lib.define_molecule_from_pdb.argtypes = [
            ctypes.c_char_p,
            ctypes.POINTER(s_molecule_c),
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_void_p),
            ctypes.POINTER(ctypes.c_void_p)]
    lib.define_molecule_from_pdb.restype = None

    lib.deallocate_s_molecule_c.argtypes = [
            ctypes.POINTER(s_molecule_c)]
    lib.deallocate_s_molecule_c.restype = None

def test():
    # ライブラリをロード
    lib_name = 'libpython_interface.so'
    lib_dir = os.path.join(os.path.dirname(__file__), '.libs')
    lib_path = os.path.join(lib_dir, lib_name)

    if not os.path.exists(lib_path):
        raise FileNotFoundError(f"Library file {lib_name} not found in {lib_dir}")

    lib = ctypes.CDLL(lib_path)

    define_prototypes(lib)

    # 関数を呼び出す
    pdb_filename = b"molecule.pdb"
    mol = s_molecule_c()
    num_atoms = ctypes.c_int()
    atom_names_ptr = ctypes.c_void_p()
    atom_coords_ptr = ctypes.c_void_p()

    lib.define_molecule_from_pdb(
            pdb_filename,
            ctypes.byref(mol),
            ctypes.byref(num_atoms),
            ctypes.byref(atom_names_ptr),
            ctypes.byref(atom_coords_ptr))

    # 結果を処理する
    print("num_atoms = ", num_atoms, type(num_atoms))
    print("num_atoms = ", mol.num_atoms, type(mol.num_atoms))
    atom_no = ctypes.cast(mol.atom_no, ctypes.POINTER(ctypes.c_int))
    for i in range(max(0, mol.num_atoms - 5), mol.num_atoms):
        print(atom_no[i])

    segment_name = ctypes.cast(mol.segment_name, ctypes.POINTER(ctypes.c_char))
    for i in range(max(0, mol.num_atoms - 5), mol.num_atoms):
        print(segment_name[i*4], segment_name[i*4+1],
              segment_name[i*4+2], segment_name[i*4+3])



    # メモリを解放する
    lib.deallocate_s_molecule_c(ctypes.byref(mol))


if __name__ == "__main__":
    test()

