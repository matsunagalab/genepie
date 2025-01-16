import ctypes
import numpy as np
import numpy.typing as npt
import c2py_util
import py2c_util
from typing import Self
import libgenesis


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
    atom_no: npt.NDArray[np.int64]
    segment_name: npt.NDArray[np.str_]
    segment_no: npt.NDArray[np.int64]
    residue_no: npt.NDArray[np.int64]
    residue_c_no: npt.NDArray[np.int64]
    residue_name: npt.NDArray[np.str_]
    atom_name: npt.NDArray[np.str_]
    atom_cls_name: npt.NDArray[np.str_]
    atom_cls_no: npt.NDArray[np.int64]
    charge: npt.NDArray[np.float64]
    mass: npt.NDArray[np.float64]
    inv_mass: npt.NDArray[np.float64]
    imove: npt.NDArray[np.int64]
    stokes_radius: npt.NDArray[np.float64]
    inv_stokes_radius: npt.NDArray[np.float64]
    chain_id: npt.NDArray[np.str_]
    atom_coord: npt.NDArray[np.float64]
    atom_occupancy: npt.NDArray[np.float64]
    atom_temp_factor: npt.NDArray[np.float64]
    atom_velocity: npt.NDArray[np.float64]
    light_atom_name: npt.NDArray[np.bool_]
    light_atom_mass: npt.NDArray[np.bool_]
    molecule_no: npt.NDArray[np.int64]
    bond_list: npt.NDArray[np.int64]
    enm_list: npt.NDArray[np.int64]
    angl_list: npt.NDArray[np.int64]
    dihe_list: npt.NDArray[np.int64]
    impr_list: npt.NDArray[np.int64]
    cmap_list: npt.NDArray[np.int64]
    molecule_atom_no: npt.NDArray[np.int64]
    molecule_mass: npt.NDArray[np.float64]
    molecule_name: npt.NDArray[np.str_]
    atom_refcoord: npt.NDArray[np.float64]
    atom_fitcoord: npt.NDArray[np.float64]
    num_pc_modes: int
    pc_mode: npt.NDArray[np.float64]
    fep_topology: int
    num_hbonds_singleA: int
    num_hbonds_singleB: int
    num_atoms_fep: npt.NDArray[np.int64]
    num_bonds_fep: npt.NDArray[np.int64]
    num_angles_fep: npt.NDArray[np.int64]
    num_dihedrals_fep: npt.NDArray[np.int64]
    num_impropers_fep: npt.NDArray[np.int64]
    num_cmaps_fep: npt.NDArray[np.int64]
    bond_list_fep: npt.NDArray[np.int64]
    angl_list_fep: npt.NDArray[np.int64]
    dihe_list_fep: npt.NDArray[np.int64]
    impr_list_fep: npt.NDArray[np.int64]
    cmap_list_fep: npt.NDArray[np.int64]
    id_singleA: npt.NDArray[np.int64]
    id_singleB: npt.NDArray[np.int64]
    fepgrp: npt.NDArray[np.int64]
    fepgrp_bond: npt.NDArray[np.int64]
    fepgrp_angl: npt.NDArray[np.int64]
    fepgrp_dihe: npt.NDArray[np.int64]
    fepgrp_cmap: npt.NDArray[np.int64]


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
    dst.atom_no = c2py_util.conv_int_ndarray(src.atom_no, dst.num_atoms)
    dst.segment_name = c2py_util.conv_fixed_length_string_ndarray(src.segment_name, (dst.num_atoms, 4))
    dst.segment_no = c2py_util.conv_int_ndarray(src.segment_no, dst.num_atoms)
    dst.residue_no = c2py_util.conv_int_ndarray(src.residue_no, dst.num_atoms)
    dst.residue_c_no = c2py_util.conv_int_ndarray(src.residue_c_no, dst.num_atoms)
    dst.residue_name = c2py_util.conv_fixed_length_string_ndarray(src.residue_name, (dst.num_atoms, 6))
    dst.atom_name = c2py_util.conv_fixed_length_string_ndarray(src.atom_name, (dst.num_atoms, 4))
    dst.atom_cls_name = c2py_util.conv_fixed_length_string_ndarray(src.atom_cls_name, (dst.num_atoms, 6))
    dst.atom_cls_no = c2py_util.conv_int_ndarray(src.atom_cls_no, dst.num_atoms)
    dst.charge = c2py_util.conv_double_ndarray(src.charge, dst.num_atoms)
    dst.mass = c2py_util.conv_double_ndarray(src.mass, dst.num_atoms)
    dst.inv_mass = c2py_util.conv_double_ndarray(src.inv_mass, dst.num_atoms)
    dst.imove = c2py_util.conv_int_ndarray(src.imove, dst.num_atoms)
    dst.stokes_radius = c2py_util.conv_double_ndarray(src.stokes_radius, dst.num_atoms)
    dst.inv_stokes_radius = c2py_util.conv_double_ndarray(src.inv_stokes_radius, dst.num_atoms)

    dst.chain_id = c2py_util.conv_fixed_length_string_ndarray(src.chain_id, (dst.num_atoms, 1))
    dst.atom_coord = c2py_util.conv_double_ndarray(src.atom_coord, (src.num_atoms, 3))
    dst.atom_occupancy = c2py_util.conv_double_ndarray(src.atom_occupancy, dst.num_atoms)
    dst.atom_temp_factor = c2py_util.conv_double_ndarray(src.atom_temp_factor, dst.num_atoms)
    dst.atom_velocity = c2py_util.conv_double_ndarray(src.atom_velocity, (src.num_atoms, 3))
    dst.light_atom_name = c2py_util.conv_bool_ndarray(src.light_atom_name, src.num_atoms)
    dst.light_atom_mass = c2py_util.conv_bool_ndarray(src.light_atom_mass, src.num_atoms)

    dst.molecule_no = c2py_util.conv_int_ndarray(src.molecule_no, dst.num_atoms)

    dst.bond_list = c2py_util.conv_int_ndarray(src.bond_list, (src.num_bonds, 2))
    dst.enm_list = c2py_util.conv_int_ndarray(src.enm_list, (src.num_enm_bonds, 2))
    dst.angl_list = c2py_util.conv_int_ndarray(src.angl_list, (src.num_angles, 3))
    dst.dihe_list = c2py_util.conv_int_ndarray(src.dihe_list, (src.num_dihedrals, 4))
    dst.impr_list = c2py_util.conv_int_ndarray(src.impr_list, (src.num_impropers, 4))
    dst.cmap_list = c2py_util.conv_int_ndarray(src.cmap_list, (src.num_cmaps, 8))

    dst.molecule_atom_no = c2py_util.conv_int_ndarray(src.molecule_atom_no, dst.num_molecules)
    dst.molecule_mass = c2py_util.conv_double_ndarray(src.molecule_mass, dst.num_molecules)
    dst.molecule_name = c2py_util.conv_fixed_length_string_ndarray(src.molecule_name, (dst.num_molecules, 10))
    dst.atom_refcoord = c2py_util.conv_double_ndarray(src.atom_refcoord, (src.num_atoms, 3))
    dst.atom_fitcoord = c2py_util.conv_double_ndarray(src.atom_fitcoord, (src.num_atoms, 3))

    dst.num_pc_modes = int(src.num_pc_modes)
    dst.pc_mode = c2py_util.conv_double_ndarray(src.pc_mode, dst.num_pc_modes)

    dst.fep_topology = int(src.fep_topology)
    dst.num_hbonds_singleA = int(src.num_hbonds_singleA)
    dst.num_hbonds_singleB = int(src.num_hbonds_singleB)

    dst.num_atoms_fep = int(src.num_atoms_fep)
    dst.num_bonds_fep = int(src.num_bonds_fep)
    dst.num_angles_fep = int(src.num_angles_fep)
    dst.num_dihedrals_fep = int(src.num_dihedrals_fep)
    dst.num_impropers_fep = int(src.num_impropers_fep)
    dst.num_cmaps_fep = int(src.num_cmaps_fep)

    dst.bond_list_fep = c2py_util.conv_int_ndarray(src.bond_list_fep, (5, src.nbnd_fep_max, 2))
    dst.angl_list_fep = c2py_util.conv_int_ndarray(src.angl_list_fep, (5, src.nangl_fep_max, 3))
    dst.dihe_list_fep = c2py_util.conv_int_ndarray(src.dihe_list_fep, (5, src.ndihe_fep_max, 4))
    dst.impr_list_fep = c2py_util.conv_int_ndarray(src.impr_list_fep, (5, src.nimpr_fep_max, 4))
    dst.cmap_list_fep = c2py_util.conv_int_ndarray(src.cmap_list_fep, (5, src.ncmap_fep_max, 8))

    dst.id_singleA = c2py_util.conv_int_ndarray(src.id_singleA, src.size_id_singleA)
    dst.id_singleB = c2py_util.conv_int_ndarray(src.id_singleB, src.size_id_singleB)
    dst.fepgrp = c2py_util.conv_int_ndarray(src.fepgrp, src.size_fepgrp)
    dst.fepgrp_bond = c2py_util.conv_int_ndarray(src.fepgrp_bond, (5, 5))
    dst.fepgrp_angl = c2py_util.conv_int_ndarray(src.fepgrp_angl, (5, 5, 5))
    dst.fepgrp_dihe = c2py_util.conv_int_ndarray(src.fepgrp_dihe, (5, 5, 5, 5))
    dst.fepgrp_cmap = c2py_util.conv_int_ndarray(src.fepgrp_cmap, 5*5*5*5*5*5*5*5)
    return dst


def py2c_s_molecule(src: SMolecule) -> SMoleculeC:
    dst = SMoleculeC()
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
    dst.num_pc_modes = src.num_pc_modes
    dst.fep_topology = src.fep_topology
    dst.num_hbonds_singleA = src.num_hbonds_singleA
    dst.num_hbonds_singleB = src.num_hbonds_singleB
    libgenesis.LibGenesis().lib.allocate_s_molecule_c(dst)

    dst.shift_origin  = src.shift_origin
    dst.special_hydrogen = src.special_hydrogen
    dst.total_charge  = src.total_charge
    py2c_util.write_int_ndarray(src.atom_no, dst.atom_no)
    py2c_util.write_fixed_length_string_ndarray(
            src.segment_name, dst.segment_name)
    py2c_util.write_int_ndarray(src.segment_no, dst.segment_no)
    py2c_util.write_int_ndarray(src.residue_no, dst.residue_no)
    py2c_util.write_int_ndarray(src.residue_c_no, dst.residue_c_no)
    py2c_util.write_fixed_length_string_ndarray(
            src.residue_name, dst.residue_name)
    py2c_util.write_fixed_length_string_ndarray(
            src.atom_name, dst.atom_name)
    py2c_util.write_fixed_length_string_ndarray(
            src.atom_cls_name, dst.atom_cls_name)
    py2c_util.write_int_ndarray(src.atom_cls_no, dst.atom_cls_no)
    py2c_util.write_double_ndarray(src.charge, dst.charge)
    py2c_util.write_double_ndarray(src.mass, dst.mass)
    py2c_util.write_double_ndarray(src.inv_mass, dst.inv_mass)
    py2c_util.write_int_ndarray(src.imove, dst.imove)
    py2c_util.write_double_ndarray(src.stokes_radius, dst.stokes_radius)
    py2c_util.write_double_ndarray(src.inv_stokes_radius, dst.inv_stokes_radius)

    py2c_util.write_fixed_length_string_ndarray(src.chain_id, dst.chain_id)
    py2c_util.write_double_ndarray(src.atom_coord, dst.atom_coord)
    py2c_util.write_double_ndarray(src.atom_occupancy, dst.atom_occupancy)
    py2c_util.write_double_ndarray(src.atom_temp_factor, dst.atom_temp_factor)
    py2c_util.write_double_ndarray(src.atom_velocity, dst.atom_velocity)
    py2c_util.write_bool_ndarray(src.light_atom_name, dst.light_atom_name)
    py2c_util.write_bool_ndarray(src.light_atom_mass, dst.light_atom_mass)
    py2c_util.write_int_ndarray(src.molecule_no, dst.molecule_no)
    py2c_util.write_int_ndarray(src.bond_list, dst.bond_list)
    py2c_util.write_int_ndarray(src.enm_list, dst.enm_list)
    py2c_util.write_int_ndarray(src.angl_list, dst.angl_list)
    py2c_util.write_int_ndarray(src.dihe_list, dst.dihe_list)
    py2c_util.write_int_ndarray(src.impr_list, dst.impr_list)
    py2c_util.write_int_ndarray(src.cmap_list, dst.cmap_list)
    py2c_util.write_int_ndarray(src.molecule_atom_no, dst.molecule_atom_no)
    py2c_util.write_double_ndarray(src.molecule_mass, dst.molecule_mass)
    py2c_util.write_fixed_length_string_ndarray(src.molecule_name, dst.molecule_name)
    py2c_util.write_double_ndarray(src.atom_refcoord, dst.atom_refcoord)
    py2c_util.write_double_ndarray(src.atom_fitcoord, dst.atom_fitcoord)

    py2c_util.write_double_ndarray(src.pc_mode, dst.pc_mode)
    py2c_util.write_int_ndarray(src.num_atoms_fep, dst.num_atoms_fep)
    py2c_util.write_int_ndarray(src.num_bonds_fep, dst.num_bonds_fep)
    py2c_util.write_int_ndarray(src.num_angles_fep, dst.num_angles_fep)
    py2c_util.write_int_ndarray(src.num_dihedrals_fep, dst.num_dihedrals_fep)
    py2c_util.write_int_ndarray(src.num_impropers_fep, dst.num_impropers_fep)
    py2c_util.write_int_ndarray(src.num_cmaps_fep, dst.num_cmaps_fep)
    py2c_util.write_int_ndarray(src.bond_list_fep, dst.bond_list_fep)
    py2c_util.write_int_ndarray(src.angl_list_fep, dst.angl_list_fep)
    py2c_util.write_int_ndarray(src.dihe_list_fep, dst.dihe_list_fep)
    py2c_util.write_int_ndarray(src.impr_list_fep, dst.impr_list_fep)
    py2c_util.write_int_ndarray(src.cmap_list_fep, dst.cmap_list_fep)
    py2c_util.write_int_ndarray(src.id_singleA, dst.id_singleA)
    py2c_util.write_int_ndarray(src.id_singleB, dst.id_singleB)
    py2c_util.write_int_ndarray(src.fepgrp, dst.fepgrp)
    py2c_util.write_int_ndarray(src.fepgrp_bond, dst.fepgrp_bond)
    py2c_util.write_int_ndarray(src.fepgrp_angl, dst.fepgrp_angl)
    py2c_util.write_int_ndarray(src.fepgrp_dihe, dst.fepgrp_dihe)
    py2c_util.write_int_ndarray(src.fepgrp_cmap, dst.fepgrp_cmap)
    return dst
