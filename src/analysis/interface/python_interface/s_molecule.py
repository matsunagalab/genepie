import ctypes
import c2py_util


class SMolecule:
    num_deg_freedom: int
    num_atoms: int
    num_bonds: int
    num_enm_bonds: int
    num_angles: int
    num_dihedrals: int
    num_impropers: int
    num_cmaps: int
    num_residues: int
    num_molecules: int
    num_segments: int
    shift_origin: bool
    special_hydrogen: bool
    total_charge: float
    atom_no: list[int]
    segment_name: list[str]
    segment_no: list[int]
    residue_no: list[int]
    residue_c_no: list[int]
    residue_name: list[str]
    atom_name: list[str]
    atom_cls_name: list[str]
    atom_cls_no: list[int]
    charge: list[float]
    mass: list[float]
    inv_mass: list[float]
    imove: list[int]
    stokes_radius: list[float]
    inv_stokes_radius: list[float]
    chain_id: list[str]
    atom_coord: list[list[float]]
    atom_occupancy: list[float]
    atom_temp_factor: list[float]
    atom_velocity: list[list[float]]
    light_atom_name: list[bool]
    light_atom_mass: list[bool]
    molecule_no: list[int]
    bond_list: list[list[int]]
    enm_list: list[list[int]]
    angl_list: list[list[int]]
    dihe_list: list[list[int]]
    impr_list: list[list[int]]
    cmap_list: list[list[int]]
    molecule_atom_no: list[int]
    molecule_mass: list[float]
    molecule_name: list[str]
    atom_refcoord: list[list[float]]
    atom_fitcoord: list[list[float]]
    num_pc_modes: int
    pc_mode: list[float]
    fep_topology: int
    num_hbonds_singleA: int
    num_hbonds_singleB: int
    num_atoms_fep: list[int]
    num_bonds_fep: list[int]
    num_angles_fep: list[int]
    num_dihedrals_fep: list[int]
    num_impropers_fep: list[int]
    num_cmaps_fep: list[int]
    bond_list_fep: list[list[list[int]]]
    angl_list_fep: list[list[list[int]]]
    dihe_list_fep: list[list[list[int]]]
    impr_list_fep: list[list[list[int]]]
    cmap_list_fep: list[list[list[int]]]
    id_singleA: list[int]
    id_singleB: list[int]
    fepgrp: list[int]
    fepgrp_bond: list[list[int]]
    fepgrp_angl: list[list[list[int]]]
    fepgrp_dihe: list[list[list[list[int]]]]
    fepgrp_cmap: list[int]


class SMoleculeC(ctypes.Structure):
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
                 ("nbnd_fep_max", ctypes.c_int),
                 ("nangl_fep_max", ctypes.c_int),
                 ("ndihe_fep_max", ctypes.c_int),
                 ("nimpr_fep_max", ctypes.c_int),
                 ("ncmap_fep_max", ctypes.c_int),
                 ("size_id_singleA", ctypes.c_int),
                 ("size_id_singleB", ctypes.c_int),
                 ("size_fepgrp", ctypes.c_int),
                 ]


def c2py_s_molecule(src: SMoleculeC) -> SMolecule:
    dst = SMolecule()
    dst.num_deg_freedom = src.num_deg_freedom
    dst.num_atoms     = src.num_atoms
    dst.num_bonds     = src.num_bonds
    dst.num_enm_bonds = src.num_enm_bonds
    dst.num_angles    = src.num_angles
    dst.num_dihedrals = src.num_dihedrals
    dst.num_impropers = src.num_impropers
    dst.num_cmaps     = src.num_cmaps
    dst.num_residues  = src.num_residues
    dst.num_molecules = src.num_molecules
    dst.num_segments  = src.num_segments
    dst.shift_origin  = src.shift_origin
    dst.special_hydrogen = src.special_hydrogen
    dst.total_charge  = src.total_charge
    dst.atom_no = c2py_util.conv_int_array(src.atom_no, dst.num_atoms)
    dst.segment_name = c2py_util.conv_str_array(src.segment_name, dst.num_atoms, 4)
    dst.segment_no = c2py_util.conv_int_array(src.segment_no, dst.num_atoms)
    dst.residue_no = c2py_util.conv_int_array(src.residue_no, dst.num_atoms)
    dst.residue_c_no = c2py_util.conv_int_array(src.residue_c_no, dst.num_atoms)
    dst.residue_name = c2py_util.conv_str_array(src.residue_name, dst.num_atoms, 6)
    dst.atom_name = c2py_util.conv_str_array(src.atom_name, dst.num_atoms, 4)
    dst.atom_cls_name = c2py_util.conv_str_array(src.atom_cls_name, dst.num_atoms, 6)
    dst.atom_cls_no = c2py_util.conv_int_array(src.atom_cls_no, dst.num_atoms)
    dst.charge = c2py_util.conv_double_array(src.charge, dst.num_atoms)
    dst.mass = c2py_util.conv_double_array(src.mass, dst.num_atoms)
    dst.inv_mass = c2py_util.conv_double_array(src.inv_mass, dst.num_atoms)
    dst.imove = c2py_util.conv_int_array(src.imove, dst.num_atoms)
    dst.stokes_radius = c2py_util.conv_double_array(src.stokes_radius, dst.num_atoms)
    dst.inv_stokes_radius = c2py_util.conv_double_array(src.inv_stokes_radius, dst.num_atoms)

    dst.chain_id = c2py_util.conv_str_array(src.chain_id, dst.num_atoms, 1)
    dst.atom_coord = c2py_util.conv_double_array_2d(src.atom_coord, src.num_atoms, 3)
    dst.atom_occupancy = c2py_util.conv_double_array(src.atom_occupancy, dst.num_atoms)
    dst.atom_temp_factor = c2py_util.conv_double_array(src.atom_temp_factor, dst.num_atoms)
    dst.atom_velocity = c2py_util.conv_double_array_2d(src.atom_velocity, src.num_atoms, 3)
    dst.light_atom_name = c2py_util.conv_bool_array(src.light_atom_name, src.num_atoms)
    dst.light_atom_mass = c2py_util.conv_bool_array(src.light_atom_mass, src.num_atoms)

    dst.molecule_no = c2py_util.conv_int_array(src.molecule_no, dst.num_atoms)

    dst.bond_list = c2py_util.conv_int_array_2d(src.bond_list, src.num_bonds, 2)
    dst.enm_list = c2py_util.conv_int_array_2d(src.enm_list, src.num_enm_bonds, 2)
    dst.angl_list = c2py_util.conv_int_array_2d(src.angl_list, src.num_angles, 3)
    dst.dihe_list = c2py_util.conv_int_array_2d(src.dihe_list, src.num_dihedrals, 4)
    dst.impr_list = c2py_util.conv_int_array_2d(src.impr_list, src.num_impropers, 4)
    dst.cmap_list = c2py_util.conv_int_array_2d(src.cmap_list, src.num_cmaps, 8)

    dst.molecule_atom_no = c2py_util.conv_int_array(src.molecule_atom_no, dst.num_molecules)
    dst.molecule_mass = c2py_util.conv_double_array(src.molecule_mass, dst.num_molecules)
    dst.molecule_name = c2py_util.conv_str_array(src.molecule_name, dst.num_molecules, 10)
    dst.atom_refcoord = c2py_util.conv_double_array_2d(src.atom_refcoord, src.num_atoms, 3)
    dst.atom_fitcoord = c2py_util.conv_double_array_2d(src.atom_fitcoord, src.num_atoms, 3)

    dst.num_pc_modes = int(src.num_pc_modes)
    dst.pc_mode = c2py_util.conv_double_array(src.pc_mode, dst.num_pc_modes)

    dst.fep_topology = int(src.fep_topology)
    dst.num_hbonds_singleA = int(src.num_hbonds_singleA)
    dst.num_hbonds_singleB = int(src.num_hbonds_singleB)

    dst.num_atoms_fep = int(src.num_atoms_fep)
    dst.num_bonds_fep = int(src.num_bonds_fep)
    dst.num_angles_fep = int(src.num_angles_fep)
    dst.num_dihedrals_fep = int(src.num_dihedrals_fep)
    dst.num_impropers_fep = int(src.num_impropers_fep)
    dst.num_cmaps_fep = int(src.num_cmaps_fep)

    dst.bond_list_fep = c2py_util.conv_int_array_3d(src.bond_list_fep, 5, src.nbnd_fep_max, 2)
    dst.angl_list_fep = c2py_util.conv_int_array_3d(src.angl_list_fep, 5, src.nangl_fep_max, 3)
    dst.dihe_list_fep = c2py_util.conv_int_array_3d(src.dihe_list_fep, 5, src.ndihe_fep_max, 4)
    dst.impr_list_fep = c2py_util.conv_int_array_3d(src.impr_list_fep, 5, src.nimpr_fep_max, 4)
    dst.cmap_list_fep = c2py_util.conv_int_array_3d(src.cmap_list_fep, 5, src.ncmap_fep_max, 8)

    dst.id_singleA = c2py_util.conv_int_array(src.id_singleA, src.size_id_singleA)
    dst.id_singleB = c2py_util.conv_int_array(src.id_singleB, src.size_id_singleB)
    dst.fepgrp = c2py_util.conv_int_array(src.fepgrp, src.size_fepgrp)
    dst.fepgrp_bond = c2py_util.conv_int_array_2d(src.fepgrp_bond, 5, 5)
    dst.fepgrp_angl = c2py_util.conv_int_array_3d(src.fepgrp_angl, 5, 5, 5)
    dst.fepgrp_dihe = c2py_util.conv_int_array_4d(src.fepgrp_dihe, 5, 5, 5, 5)
    dst.fepgrp_cmap = c2py_util.conv_int_array(src.fepgrp_cmap, 5*5*5*5*5*5*5*5)
    return dst
