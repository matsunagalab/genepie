import ctypes
import os
import tempfile
import numpy as np
import numpy.typing as npt
import c2py_util
import py2c_util
from typing import Self, Union
from libgenesis import LibGenesis
from s_molecule_c import SMoleculeC


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
    segment_name: npt.NDArray[np.object_]
    segment_no: npt.NDArray[np.int64]
    residue_no: npt.NDArray[np.int64]
    residue_c_no: npt.NDArray[np.int64]
    residue_name: npt.NDArray[np.object_]
    atom_name: npt.NDArray[np.object_]
    atom_cls_name: npt.NDArray[np.object_]
    atom_cls_no: npt.NDArray[np.int64]
    charge: npt.NDArray[np.float64]
    mass: npt.NDArray[np.float64]
    inv_mass: npt.NDArray[np.float64]
    imove: npt.NDArray[np.int64]
    stokes_radius: npt.NDArray[np.float64]
    inv_stokes_radius: npt.NDArray[np.float64]
    chain_id: npt.NDArray[np.object_]
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
    molecule_name: npt.NDArray[np.object_]
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

    def __init__(self) -> Self:
        self.num_deg_freedom = 0
        self.num_atoms = 0
        self.num_bonds = 0
        self.num_enm_bonds = 0
        self.num_angles = 0
        self.num_dihedrals = 0
        self.num_impropers = 0
        self.num_cmaps = 0
        self.num_residues = 0
        self.num_molecules = 0
        self.num_segments = 0
        self.shift_origin = False
        self.special_hydrogen = False
        self.total_charge = 0.0
        self.atom_no = np.array([])
        self.segment_name = np.array([])
        self.segment_no = np.array([])
        self.residue_no = np.array([])
        self.residue_c_no = np.array([])
        self.residue_name = np.array([])
        self.atom_name = np.array([])
        self.atom_cls_name = np.array([])
        self.atom_cls_no = np.array([])
        self.charge = np.array([])
        self.mass = np.array([])
        self.inv_mass = np.array([])
        self.imove = np.array([])
        self.stokes_radius = np.array([])
        self.inv_stokes_radius = np.array([])
        self.chain_id = np.array([])
        self.atom_coord = np.array([])
        self.atom_occupancy = np.array([])
        self.atom_temp_factor = np.array([])
        self.atom_velocity = np.array([])
        self.light_atom_name = np.array([])
        self.light_atom_mass = np.array([])
        self.molecule_no = np.array([])
        self.bond_list = np.array([])
        self.enm_list = np.array([])
        self.angl_list = np.array([])
        self.dihe_list = np.array([])
        self.impr_list = np.array([])
        self.cmap_list = np.array([])
        self.molecule_atom_no = np.array([])
        self.molecule_mass = np.array([])
        self.molecule_name = np.array([])
        self.atom_refcoord = np.array([])
        self.atom_fitcoord = np.array([])
        self.num_pc_modes = 0
        self.pc_mode = np.array([])
        self.fep_topology = 0
        self.num_hbonds_singleA = 0
        self.num_hbonds_singleB = 0
        self.num_atoms_fep = np.array([])
        self.num_bonds_fep = np.array([])
        self.num_angles_fep = np.array([])
        self.num_dihedrals_fep = np.array([])
        self.num_impropers_fep = np.array([])
        self.num_cmaps_fep = np.array([])
        self.bond_list_fep = np.array([])
        self.angl_list_fep = np.array([])
        self.dihe_list_fep = np.array([])
        self.impr_list_fep = np.array([])
        self.cmap_list_fep = np.array([])
        self.id_singleA = np.array([])
        self.id_singleB = np.array([])
        self.fepgrp = np.array([])
        self.fepgrp_bond = np.array([])
        self.fepgrp_angl = np.array([])
        self.fepgrp_dihe = np.array([])
        self.fepgrp_cmap = np.array([])

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.free()

    def free(self):
        """deallocate resources"""
        pass

    @staticmethod
    def from_SMoleculeC(src: SMoleculeC) -> Self:
        dst = SMolecule()
        dst.num_deg_freedom = src.num_deg_freedom
        dst.num_atoms = src.num_atoms
        dst.num_bonds = src.num_bonds
        dst.num_enm_bonds = src.num_enm_bonds
        dst.num_angles = src.num_angles
        dst.num_dihedrals = src.num_dihedrals
        dst.num_impropers = src.num_impropers
        dst.num_cmaps = src.num_cmaps
        dst.num_residues = src.num_residues
        dst.num_molecules = src.num_molecules
        dst.num_segments = src.num_segments
        dst.shift_origin = src.shift_origin
        dst.special_hydrogen = src.special_hydrogen
        dst.total_charge = src.total_charge
        dst.atom_no = c2py_util.conv_int_ndarray(src.atom_no, dst.num_atoms)
        dst.segment_name = c2py_util.conv_pystring_ndarray(
            src.segment_name, (dst.num_atoms, 4)
        )
        dst.segment_no = c2py_util.conv_int_ndarray(
                src.segment_no, dst.num_atoms)
        dst.residue_no = c2py_util.conv_int_ndarray(
                src.residue_no, dst.num_atoms)
        dst.residue_c_no = c2py_util.conv_int_ndarray(
                src.residue_c_no, dst.num_atoms)
        dst.residue_name = c2py_util.conv_pystring_ndarray(
            src.residue_name, (dst.num_atoms, 6)
        )
        dst.atom_name = c2py_util.conv_pystring_ndarray(
            src.atom_name, (dst.num_atoms, 4)
        )
        dst.atom_cls_name = c2py_util.conv_pystring_ndarray(
            src.atom_cls_name, (dst.num_atoms, 6)
        )
        dst.atom_cls_no = c2py_util.conv_int_ndarray(
                src.atom_cls_no, dst.num_atoms)
        dst.charge = c2py_util.conv_double_ndarray(src.charge, dst.num_atoms)
        dst.mass = c2py_util.conv_double_ndarray(src.mass, dst.num_atoms)
        dst.inv_mass = c2py_util.conv_double_ndarray(
                src.inv_mass, dst.num_atoms)
        dst.imove = c2py_util.conv_int_ndarray(src.imove, dst.num_atoms)
        dst.stokes_radius = c2py_util.conv_double_ndarray(
                src.stokes_radius, dst.num_atoms)
        dst.inv_stokes_radius = c2py_util.conv_double_ndarray(
            src.inv_stokes_radius, dst.num_atoms)

        dst.chain_id = c2py_util.conv_pystring_ndarray(
            src.chain_id, (dst.num_atoms, 1)
        )
        dst.atom_coord = c2py_util.conv_double_ndarray(
                src.atom_coord, (src.num_atoms, 3))
        dst.atom_occupancy = c2py_util.conv_double_ndarray(
            src.atom_occupancy, dst.num_atoms)
        dst.atom_temp_factor = c2py_util.conv_double_ndarray(
            src.atom_temp_factor, dst.num_atoms)
        dst.atom_velocity = c2py_util.conv_double_ndarray(
            src.atom_velocity, (src.num_atoms, 3))
        dst.light_atom_name = c2py_util.conv_bool_ndarray(
            src.light_atom_name, src.num_atoms)
        dst.light_atom_mass = c2py_util.conv_bool_ndarray(
            src.light_atom_mass, src.num_atoms)

        dst.molecule_no = c2py_util.conv_int_ndarray(
                src.molecule_no, dst.num_atoms)
        dst.bond_list = c2py_util.conv_int_ndarray(
                src.bond_list, (src.num_bonds, 2))
        dst.enm_list = c2py_util.conv_int_ndarray(
                src.enm_list, (src.num_enm_bonds, 2))
        dst.angl_list = c2py_util.conv_int_ndarray(
                src.angl_list, (src.num_angles, 3))
        dst.dihe_list = c2py_util.conv_int_ndarray(
                src.dihe_list, (src.num_dihedrals, 4))
        dst.impr_list = c2py_util.conv_int_ndarray(
                src.impr_list, (src.num_impropers, 4))
        dst.cmap_list = c2py_util.conv_int_ndarray(
                src.cmap_list, (src.num_cmaps, 8))

        dst.molecule_atom_no = c2py_util.conv_int_ndarray(
            src.molecule_atom_no, dst.num_molecules)
        dst.molecule_mass = c2py_util.conv_double_ndarray(
            src.molecule_mass, dst.num_molecules)
        dst.molecule_name = c2py_util.conv_pystring_ndarray(
            src.molecule_name, (dst.num_molecules, 10))
        dst.atom_refcoord = c2py_util.conv_double_ndarray(
            src.atom_refcoord, (src.num_atoms, 3))
        dst.atom_fitcoord = c2py_util.conv_double_ndarray(
            src.atom_fitcoord, (src.num_atoms, 3))

        dst.num_pc_modes = int(src.num_pc_modes)
        dst.pc_mode = c2py_util.conv_double_ndarray(
                src.pc_mode, dst.num_pc_modes)

        dst.fep_topology = int(src.fep_topology)
        dst.num_hbonds_singleA = int(src.num_hbonds_singleA)
        dst.num_hbonds_singleB = int(src.num_hbonds_singleB)

        dst.num_atoms_fep = int(src.num_atoms_fep)
        dst.num_bonds_fep = int(src.num_bonds_fep)
        dst.num_angles_fep = int(src.num_angles_fep)
        dst.num_dihedrals_fep = int(src.num_dihedrals_fep)
        dst.num_impropers_fep = int(src.num_impropers_fep)
        dst.num_cmaps_fep = int(src.num_cmaps_fep)

        dst.bond_list_fep = c2py_util.conv_int_ndarray(
            src.bond_list_fep, (5, src.nbnd_fep_max, 2))
        dst.angl_list_fep = c2py_util.conv_int_ndarray(
            src.angl_list_fep, (5, src.nangl_fep_max, 3))
        dst.dihe_list_fep = c2py_util.conv_int_ndarray(
            src.dihe_list_fep, (5, src.ndihe_fep_max, 4))
        dst.impr_list_fep = c2py_util.conv_int_ndarray(
            src.impr_list_fep, (5, src.nimpr_fep_max, 4))
        dst.cmap_list_fep = c2py_util.conv_int_ndarray(
            src.cmap_list_fep, (5, src.ncmap_fep_max, 8))

        dst.id_singleA = c2py_util.conv_int_ndarray(
                src.id_singleA, src.size_id_singleA)
        dst.id_singleB = c2py_util.conv_int_ndarray(
                src.id_singleB, src.size_id_singleB)
        dst.fepgrp = c2py_util.conv_int_ndarray(src.fepgrp, src.size_fepgrp)
        dst.fepgrp_bond = c2py_util.conv_int_ndarray(src.fepgrp_bond, (5, 5))
        dst.fepgrp_angl = c2py_util.conv_int_ndarray(
                src.fepgrp_angl, (5, 5, 5))
        dst.fepgrp_dihe = c2py_util.conv_int_ndarray(
                src.fepgrp_dihe, (5, 5, 5, 5))
        dst.fepgrp_cmap = c2py_util.conv_int_ndarray(
            src.fepgrp_cmap, 5 * 5 * 5 * 5 * 5 * 5 * 5 * 5)
        return dst

    def to_SMoleculeC(self) -> SMoleculeC:
        dst = SMoleculeC()
        dst.num_deg_freedom = self.num_deg_freedom
        dst.num_atoms = self.num_atoms
        dst.num_bonds = self.num_bonds
        dst.num_enm_bonds = self.num_enm_bonds
        dst.num_angles = self.num_angles
        dst.num_dihedrals = self.num_dihedrals
        dst.num_impropers = self.num_impropers
        dst.num_cmaps = self.num_cmaps
        dst.num_residues = self.num_residues
        dst.num_molecules = self.num_molecules
        dst.num_segments = self.num_segments
        dst.num_pc_modes = self.num_pc_modes
        dst.fep_topology = self.fep_topology
        dst.num_hbonds_singleA = self.num_hbonds_singleA
        dst.num_hbonds_singleB = self.num_hbonds_singleB
        LibGenesis().lib.allocate_s_molecule_c(dst)

        dst.shift_origin = self.shift_origin
        dst.special_hydrogen = self.special_hydrogen
        dst.total_charge = self.total_charge
        py2c_util.write_int_ndarray(self.atom_no, dst.atom_no)
        py2c_util.write_pystring_ndarray(self.segment_name, dst.segment_name)
        py2c_util.write_int_ndarray(self.segment_no, dst.segment_no)
        py2c_util.write_int_ndarray(self.residue_no, dst.residue_no)
        py2c_util.write_int_ndarray(self.residue_c_no, dst.residue_c_no)
        py2c_util.write_pystring_ndarray(self.residue_name, dst.residue_name)
        py2c_util.write_pystring_ndarray(self.atom_name, dst.atom_name)
        py2c_util.write_pystring_ndarray(self.atom_cls_name, dst.atom_cls_name)
        py2c_util.write_int_ndarray(self.atom_cls_no, dst.atom_cls_no)
        py2c_util.write_double_ndarray(self.charge, dst.charge)
        py2c_util.write_double_ndarray(self.mass, dst.mass)
        py2c_util.write_double_ndarray(self.inv_mass, dst.inv_mass)
        py2c_util.write_int_ndarray(self.imove, dst.imove)
        py2c_util.write_double_ndarray(self.stokes_radius, dst.stokes_radius)
        py2c_util.write_double_ndarray(
                self.inv_stokes_radius, dst.inv_stokes_radius)

        py2c_util.write_pystring_ndarray(self.chain_id, dst.chain_id)
        py2c_util.write_double_ndarray(self.atom_coord, dst.atom_coord)
        py2c_util.write_double_ndarray(self.atom_occupancy, dst.atom_occupancy)
        py2c_util.write_double_ndarray(
                self.atom_temp_factor, dst.atom_temp_factor)
        py2c_util.write_double_ndarray(self.atom_velocity, dst.atom_velocity)
        py2c_util.write_bool_ndarray(self.light_atom_name, dst.light_atom_name)
        py2c_util.write_bool_ndarray(self.light_atom_mass, dst.light_atom_mass)
        py2c_util.write_int_ndarray(self.molecule_no, dst.molecule_no)
        py2c_util.write_int_ndarray(self.bond_list, dst.bond_list)
        py2c_util.write_int_ndarray(self.enm_list, dst.enm_list)
        py2c_util.write_int_ndarray(self.angl_list, dst.angl_list)
        py2c_util.write_int_ndarray(self.dihe_list, dst.dihe_list)
        py2c_util.write_int_ndarray(self.impr_list, dst.impr_list)
        py2c_util.write_int_ndarray(self.cmap_list, dst.cmap_list)
        py2c_util.write_int_ndarray(
                self.molecule_atom_no, dst.molecule_atom_no)
        py2c_util.write_double_ndarray(self.molecule_mass, dst.molecule_mass)
        py2c_util.write_pystring_ndarray(self.molecule_name, dst.molecule_name)
        py2c_util.write_double_ndarray(self.atom_refcoord, dst.atom_refcoord)
        py2c_util.write_double_ndarray(self.atom_fitcoord, dst.atom_fitcoord)

        py2c_util.write_double_ndarray(self.pc_mode, dst.pc_mode)
        py2c_util.write_int_ndarray(self.num_atoms_fep, dst.num_atoms_fep)
        py2c_util.write_int_ndarray(self.num_bonds_fep, dst.num_bonds_fep)
        py2c_util.write_int_ndarray(self.num_angles_fep, dst.num_angles_fep)
        py2c_util.write_int_ndarray(
                self.num_dihedrals_fep, dst.num_dihedrals_fep)
        py2c_util.write_int_ndarray(
                self.num_impropers_fep, dst.num_impropers_fep)
        py2c_util.write_int_ndarray(self.num_cmaps_fep, dst.num_cmaps_fep)
        py2c_util.write_int_ndarray(self.bond_list_fep, dst.bond_list_fep)
        py2c_util.write_int_ndarray(self.angl_list_fep, dst.angl_list_fep)
        py2c_util.write_int_ndarray(self.dihe_list_fep, dst.dihe_list_fep)
        py2c_util.write_int_ndarray(self.impr_list_fep, dst.impr_list_fep)
        py2c_util.write_int_ndarray(self.cmap_list_fep, dst.cmap_list_fep)
        py2c_util.write_int_ndarray(self.id_singleA, dst.id_singleA)
        py2c_util.write_int_ndarray(self.id_singleB, dst.id_singleB)
        py2c_util.write_int_ndarray(self.fepgrp, dst.fepgrp)
        py2c_util.write_int_ndarray(self.fepgrp_bond, dst.fepgrp_bond)
        py2c_util.write_int_ndarray(self.fepgrp_angl, dst.fepgrp_angl)
        py2c_util.write_int_ndarray(self.fepgrp_dihe, dst.fepgrp_dihe)
        py2c_util.write_int_ndarray(self.fepgrp_cmap, dst.fepgrp_cmap)
        return dst

    @staticmethod
    def from_file(pdb: Union[str, bytes, os.PathLike] = '',
                  top: Union[str, bytes, os.PathLike] = '',
                  gpr: Union[str, bytes, os.PathLike] = '',
                  psf: Union[str, bytes, os.PathLike] = '',
                  ref: Union[str, bytes, os.PathLike] = '',
                  fit: Union[str, bytes, os.PathLike] = '',
                  prmtop: Union[str, bytes, os.PathLike] = '',
                  ambcrd: Union[str, bytes, os.PathLike] = '',
                  ambref: Union[str, bytes, os.PathLike] = '',
                  grotop: Union[str, bytes, os.PathLike] = '',
                  grocrd: Union[str, bytes, os.PathLike] = '',
                  groref: Union[str, bytes, os.PathLike] = '',
                  ) -> Self:
        mol_c = SMoleculeC()
        try:
            LibGenesis().lib.define_molecule_from_file(
                py2c_util.pathlike_to_byte(pdb),
                py2c_util.pathlike_to_byte(top),
                py2c_util.pathlike_to_byte(gpr),
                py2c_util.pathlike_to_byte(psf),
                py2c_util.pathlike_to_byte(ref),
                py2c_util.pathlike_to_byte(fit),
                py2c_util.pathlike_to_byte(prmtop),
                py2c_util.pathlike_to_byte(ambcrd),
                py2c_util.pathlike_to_byte(ambref),
                py2c_util.pathlike_to_byte(grotop),
                py2c_util.pathlike_to_byte(grocrd),
                py2c_util.pathlike_to_byte(groref),
                ctypes.byref(mol_c)
            )
            mol_py = SMolecule.from_SMoleculeC(mol_c)
            return mol_py
        finally:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))

    @staticmethod
    def guess_atom_element(atom_name) -> str:
        """
        Guess the element from an atom name.

        Parameters
        ----------
        atom_name :

        Returns
        -------
        md.element.Element or None
        """
        special_cases = {
            "FE": "Fe",
            "ZN": "Zn",
            "CL": "Cl",
        }
        if atom_name in special_cases:
            return special_cases[atom_name]

        stripped_name = atom_name.strip()
        if stripped_name and stripped_name[0].isalpha():
            try:
                return stripped_name[0]
            except KeyError:
                pass
        return None


try:
    import mdtraj as md

    def to_mdtraj_topology(self) -> md.Topology:
        """

        Returns
            MDTraj Topology
        -------
        """
        topology = md.Topology()
        # create atoms
        cur_chain_id = None
        cur_residue_no = None
        md_chain = None
        md_residue = None
        atom_no_to_idx = dict()
        for i in range(0, self.num_atoms):
            if cur_chain_id != self.chain_id[i]:
                cur_chain_id = self.chain_id[i]
                md_chain = topology.add_chain()
            if cur_residue_no != self.residue_no[i]:
                cur_residue_no = self.residue_no[i]
                md_residue = topology.add_residue(
                        self.residue_name[i], md_chain,
                        segment_id=self.segment_no[i])
            atom_name = "".join(self.atom_name[i])
            element = guess_atom_element_mdtraj(atom_name)
            topology.add_atom(atom_name, element, md_residue)
            atom_no_to_idx[self.atom_no[i]] = i
        # create bonds
        for i in range(0, self.num_bonds):
            atom1 = topology.atom(atom_no_to_idx[self.bond_list[i, 1]])
            atom2 = topology.atom(atom_no_to_idx[self.bond_list[i, 2]])
            topology.add_bond(atom1, atom2)
        return topology

    SMolecule.to_mdtraj_topology = to_mdtraj_topology

    def from_mdtraj_topology(src: md.Topology) -> Self:
        return from_mdtraj_topology_via_pdb(src)

    SMolecule.from_mdtraj_topology = from_mdtraj_topology

    def from_mdtraj_trajectory(src: md.Trajectory) -> Self:
        return from_mdtraj_trajectory_via_pdb(src)

    SMolecule.from_mdtraj_trajectory = from_mdtraj_trajectory

    @staticmethod
    def from_mdtraj_topology_via_mem(src: md.Topology) -> Self:
        mol = SMolecule()
        mol.num_atoms = src.n_atoms
        mol.atom_name = np.empty(src.n_atoms, dtype=np.object_)
        mol.atom_no = np.empty(src.n_atoms, dtype=np.int64)
        mol.charge = np.empty(src.n_atoms, dtype=np.float64)
        mol.residue_no = np.empty(src.n_atoms, dtype=np.int64)
        mol.residue_c_no = np.empty(src.n_atoms, dtype=np.object_)
        mol.residue_name = np.empty(src.n_atoms, dtype=np.object_)
        mol.segment_no = np.empty(src.n_atoms, dtype=np.int64)
        mol.segment_name = np.empty(src.n_atoms, dtype=np.object_)
        mol.chain_id = np.empty(src.n_atoms, dtype=np.object_)
        prev_residue = None
        for i, atom in enumerate(src.atoms):
            mol.atom_name[i] = np.object_(atom.name)
            mol.atom_no[i] = atom.index + 1
            mol.charge[i] = atom.formal_charge
            if prev_residue != atom.residue.index:
                prev_residue = atom.residue.index
                mol.num_residues += 1
            mol.residue_no[i] = atom.residue.index + 1
            mol.residue_c_no[i] = 0
            mol.residue_name[i] = np.object_(atom.residue.name)
            mol.segment_no[i] = 0
            mol.segment_name[i] = atom.segment_id
            mol.chain_id[i] = np.object_(atom.residue.chain.chain_id)
        mol.num_bonds = src.n_bonds
        mol.bond_list = np.empty((src.n_bonds, 2), dtype=np.int64)
        for i, bond in enumerate(src.bonds):
            mol.bond_list[i, 0] = bond.atom1.index + 1
            mol.bond_list[i, 1] = bond.atom2.index + 1
        mol.num_molecules = 1
        return mol

    SMolecule.from_mdtraj_topology = from_mdtraj_topology

    @staticmethod
    def from_mdtraj_trajectory_via_pdb(src: md.Trajectory) -> Self:
        with tempfile.NamedTemporaryFile(
                suffix=".pdb", dir=os.getcwd(), delete=True) as pdb_file:
            src.save_pdb(pdb_file.name)
            pdb_file.seek(0)
            return SMolecule.from_file(pdb=pdb_file.name)

    @staticmethod
    def from_mdtraj_topology_via_pdb(src: md.Topology) -> Self:
        empty_traj = md.Trajectory(xyz=np.zeros((1, src.n_atoms, 3)),
                                   topology=src)
        return from_mdtraj_trajectory_via_pdb(empty_traj)

    @staticmethod
    def guess_atom_element_mdtraj(atom_name) -> md.element.Element:
        """
        Guess the element from an atom name.

        Parameters
        ----------
        atom_name :

        Returns
        -------
        md.element.Element or None
        """
        try:
            el = SMolecule.guess_atom_element(atom_name)
            if el:
                return md.element.Element.getBySymbol(el)
        except KeyError:
            pass
        return None

    SMolecule.guess_atom_element_mdtraj = guess_atom_element_mdtraj

except ImportError:
    pass


try:
    import MDAnalysis as mda
    import MDAnalysis.core.topologyattrs as mdaattr
    from MDAnalysis.core.topology import Topology as mdaTopology

    def to_mdanalysis_topology(self) -> mdaTopology:
        """

        Returns
            MDAnalysis Topology
        -------
        """
        ids = []
        names = []
        chain_ids = []
        temp_factors = []
        elements = []
        occupancies = []
        formal_chages = []
        residue_names = []
        residue_ids = []
        seg_ids = []
        atom_resindex = []
        residue_segindex = []
        prev_resid = None
        prev_segid = None
        for i in range(0, self.num_atoms):
            ids.append(self.atom_no[i])
            names.append(self.atom_name[i])
            chain_ids.append(str(self.chain_id[i]))
            temp_factors.append(self.atom_temp_factor[i])
            elements.append(
                    SMolecule.guess_atom_element(self.atom_name[i]))
            occupancies.append(self.atom_occupancy[i])
            formal_chages.append(self.charge[i])
            if prev_segid != self.segment_name[i]:
                seg_ids.append(self.segment_name[i])
                prev_segid = self.segment_name[i]
            if prev_resid != self.residue_no[i]:
                residue_ids.append(self.residue_no[i])
                residue_names.append(self.residue_name[i])
                prev_resid = self.residue_no[i]
                residue_segindex.append(len(seg_ids) - 1)
            atom_resindex.append(len(residue_ids) - 1)

        attrs = []
        for vals, Attr, dtype in (
                (ids, mdaattr.Atomids, np.int64),
                (names, mdaattr.Atomnames, object),
                (chain_ids, mdaattr.ChainIDs, object),
                (temp_factors, mdaattr.Tempfactors, np.float64),
                (elements, mdaattr.Elements, object),
                (occupancies, mdaattr.Occupancies, np.float64),
                (formal_chages, mdaattr.FormalCharges, np.float64),
                (residue_ids, mdaattr.Resids, np.int64),
                (residue_names, mdaattr.Resnames, object),
                (seg_ids, mdaattr.Segids, object),
                ):
            attrs.append(Attr(np.array(vals, dtype=dtype)))
        top = mdaTopology(
                n_atoms=self.num_atoms,
                n_res=self.num_residues,
                n_seg=self.num_segments,
                attrs=attrs,
                atom_resindex=atom_resindex,
                residue_segindex=residue_segindex,
                )
        return top

    SMolecule.to_mdanalysis_topology = to_mdanalysis_topology

    def to_mdanalysis_universe(self) -> mda.Universe:
        top = self.to_mdanalysis_topology()
        uni = mda.Universe(top)
        return uni

    SMolecule.to_mdanalysis_universe = to_mdanalysis_universe

    @staticmethod
    def from_mdanalysis_universe(uni: mda.Universe) -> Self:
        return from_mdanalysis_universe_via_pdb(uni)

    SMolecule.from_mdanalysis_universe = from_mdanalysis_universe

    @staticmethod
    def from_mdanalysis_universe_via_mem(uni: mda.Universe) -> Self:
        mol = SMolecule()
        mol.num_atoms = len(uni.atoms)

        mol.atom_name = np.empty(mol.num_atoms, dtype=np.object_)
        mol.atom_no = np.empty(mol.num_atoms, dtype=np.int64)
        mol.charge = np.empty(mol.num_atoms, dtype=np.float64)
        mol.residue_no = np.empty(mol.num_atoms, dtype=np.int64)
        mol.residue_name = np.empty(mol.num_atoms, dtype=np.object_)
        mol.chain_id = np.empty(mol.num_atoms, dtype=np.object_)
        for i, atom in enumerate(uni.atoms):
            mol.atom_name[i] = atom.name
            mol.atom_no[i] = atom.index
            try:
                mol.charge[i] = atom.formalcharge
            except mda.exceptions.NoDataError:
                mol.charge[i] = 0.0
            mol.residue_no[i] = atom.residue.resindex
            try:
                mol.residue_name[i] = atom.residue.name
            except mda.exceptions.NoDataError:
                mol.residue_name[i] = ''
            mol.chain_id[i] = atom.chainID
        try:
            mol.num_bonds = len(uni.bonds)
            mol.bond_list = np.empty((mol.num_bonds, 2), dtype=np.int64)
            for i, bond in enumerate(uni.bonds):
                mol.bond_list[i, 0] = bond.atoms[0].index
                mol.bond_list[i, 1] = bond.atoms[1].index
        except mda.exceptions.NoDataError:
            mol.bond_list = np.empty((0, 2), dtype=np.int64)
        mol.num_molecules = 1
        return mol

    @staticmethod
    def from_mdanalysis_universe_via_pdb(uni: mda.Universe) -> Self:
        with tempfile.NamedTemporaryFile(
                suffix=".pdb", dir=os.getcwd(), delete=True) as pdb_file:
            uni.atoms.write(pdb_file.name)
            pdb_file.seek(0)
            return SMolecule.from_file(pdb=pdb_file.name)

    @staticmethod
    def from_mdanalysis_topology(top: mdaTopology) -> Self:
        return from_mdanalysis_universe(mda.Universe(topology=top))

    SMolecule.from_mdanalysis_topology = from_mdanalysis_topology

except ImportError:
    pass
