import ctypes
import os
import sys
import tempfile
from typing import Optional, Union

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

import numpy as np
import numpy.typing as npt
from . import c2py_util
from . import py2c_util
import traceback
from .libgenesis import LibGenesis
from .s_molecule_c import SMoleculeC


class SMolecule:
    """
    Python implementation that holds data equivalent to Genesis' s_molecule type.

    This class holds data equivalent to Genesis' s_molecule type.
    Data is held on the memory area allocated by python.
    """
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
        """
        Initialize an empty instance.

        Returns
        -------
        An instance that has no data.
        """
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

    @staticmethod
    def from_SMoleculeC(src: SMoleculeC) -> Self:
        """
        Create molecular structure from a Genesis's s_molecule object.

        Parameters
        ----------
        src: a Genesis's s_molecule object.

        Returns
        -------
        An instance that has same data as src.

        """
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
        """
        Create a Genesis's s_molecule object for C language interface.

        Returns
        -------
        An instance that can adapt to C language interface input.

        """
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

    def from_file(pdb= None,
                  top= None,
                  gpr= None,
                  psf= None,
                  ref= None,
                  fit= None,
                  prmtop= None,
                  ambcrd= None,
                  ambref= None,
                  grotop= None,
                  grocrd= None,
                  groref= None,
                  ) -> Self:
        """
        Load molecular structure from the specified files.

        This method initializes an object by loading various types of molecular
        structure files (e.g., PDB, topology, coordinate files).

        Parameters
        ----------
        pdb : Path to the PDB file. containing atomic coordinates.
        top : Path to the topology file.
        gpr : Path to the GPR file (GROMACS parameter file).
        psf : Path to the PSF file (CHARMM/NAMD Protein Structure File).
        ref : Path to the reference file for alignment or comparison.
        fit : Path to the fitting file for structural alignment.
        prmtop : Path to the AMBER PRMTOP file (topology file).
        ambcrd : Path to the AMBER coordinate file (CRD format).
        ambref : Path to the AMBER reference file.
        grotop : Path to the GROMACS topology file (TOP format).
        grocrd : Path to the GROMACS coordinate file (GRO format).
        groref : Path to the GROMACS reference file.

        Returns
        -------
        An instance of the class initialized with the provided file paths.

        Notes
        -----
        - At least one of the input file paths must be provided to initialize the object.

        Examples
        --------
        >>> obj = SMolecule.from_file(pdb="example.pdb", top="example.top")
        """
        paths = [pdb, top, gpr, psf, ref, fit, prmtop, ambcrd, ambref, grotop, grocrd, groref]
        labels = ['pdb', 'top', 'gpr', 'psf', 'ref', 'fit', 'prmtop', 'ambcrd', 'ambref', 'grotop', 'grocrd', 'groref']
        args = [py2c_util.pathlike_to_c_char_p(p) for p in paths]

        for i, (label, val) in enumerate(zip(labels, args)):
            print(f"arg[{i}] {label}: {val}, type = {type(val)}")

        mol_c = SMoleculeC()

        try:
            ret = LibGenesis().lib.define_molecule_from_file(
                 *args,
                 ctypes.byref(mol_c)  # 13番目の引数
            )
            print(f"define_molecule_from_file retturned: {ret}")
            if ret != 0:
                raise RuntimeError("define_molecule_from_file failed")

            mol_py = SMolecule.from_SMoleculeC(mol_c)
            return mol_py
        finally:
            if mol_c:
                LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))

    def subset_atoms(self, atom_indices: npt.NDArray[np.int64]) -> Self:
        """
        Create a subset molecule based on the given list of atom indices.

        Parameters
        ----------
        atom_indices : npt.NDArray[np.int64]
            List of atom indices to extract (0-based).

        Returns
        -------
        Self
            A new SMolecule object containing only the specified atoms.

        Notes
        -----
        - atom_indices must be 0-based indices.
        - Only bonds, angles, dihedrals, etc. between the specified atoms are retained.
        - Counts such as number of molecules and residues are updated appropriately.
        """
        if len(atom_indices) == 0:
            return SMolecule()
        
        # Sort and deduplicate indices
        atom_indices = np.unique(atom_indices)
        
        # Check index range
        if np.any(atom_indices < 0) or np.any(atom_indices >= self.num_atoms):
            raise ValueError("atom_indices contains invalid indices")
        
        # Create a new SMolecule object
        new_mol = SMolecule()
        
        # Set number of atoms and other basic information
        new_mol.num_atoms = len(atom_indices)
        new_mol.shift_origin = self.shift_origin
        new_mol.special_hydrogen = self.special_hydrogen
        
        # Extract atom information
        new_mol.atom_no = self.atom_no[atom_indices].copy()
        new_mol.segment_name = self.segment_name[atom_indices].copy()
        new_mol.segment_no = self.segment_no[atom_indices].copy()
        new_mol.residue_no = self.residue_no[atom_indices].copy()
        new_mol.residue_c_no = self.residue_c_no[atom_indices].copy()
        new_mol.residue_name = self.residue_name[atom_indices].copy()
        new_mol.atom_name = self.atom_name[atom_indices].copy()
        new_mol.atom_cls_name = self.atom_cls_name[atom_indices].copy()
        new_mol.atom_cls_no = self.atom_cls_no[atom_indices].copy()
        new_mol.charge = self.charge[atom_indices].copy()
        new_mol.mass = self.mass[atom_indices].copy()
        new_mol.inv_mass = self.inv_mass[atom_indices].copy()
        new_mol.imove = self.imove[atom_indices].copy()
        new_mol.stokes_radius = self.stokes_radius[atom_indices].copy()
        new_mol.inv_stokes_radius = self.inv_stokes_radius[atom_indices].copy()
        new_mol.chain_id = self.chain_id[atom_indices].copy()
        new_mol.atom_coord = self.atom_coord[atom_indices].copy()
        new_mol.atom_occupancy = self.atom_occupancy[atom_indices].copy()
        new_mol.atom_temp_factor = self.atom_temp_factor[atom_indices].copy()
        new_mol.atom_velocity = self.atom_velocity[atom_indices].copy()
        new_mol.light_atom_name = self.light_atom_name[atom_indices].copy()
        new_mol.light_atom_mass = self.light_atom_mass[atom_indices].copy()
        new_mol.molecule_no = self.molecule_no[atom_indices].copy()
        
        # Update reference and fit coordinates if present
        if self.atom_refcoord.size > 0:
            new_mol.atom_refcoord = self.atom_refcoord[atom_indices].copy()
        if self.atom_fitcoord.size > 0:
            new_mol.atom_fitcoord = self.atom_fitcoord[atom_indices].copy()
        
        # Create mapping from old index to new index (original index -> new index)
        old_to_new = {old_idx: new_idx for new_idx, old_idx in enumerate(atom_indices)}
        
        # Update bond information
        if self.num_bonds > 0 and self.bond_list.size > 0:
            valid_bonds = []
            for i in range(self.num_bonds):
                atom1, atom2 = self.bond_list[i, 0] - 1, self.bond_list[i, 1] - 1  # 1-based to 0-based
                if atom1 in old_to_new and atom2 in old_to_new:
                    valid_bonds.append([old_to_new[atom1] + 1, old_to_new[atom2] + 1])  # 0-based to 1-based
            
            if valid_bonds:
                new_mol.num_bonds = len(valid_bonds)
                new_mol.bond_list = np.array(valid_bonds, dtype=np.int64)
            else:
                new_mol.num_bonds = 0
                new_mol.bond_list = np.array([], dtype=np.int64).reshape(0, 2)
        
        # Update ENM bond information
        if self.num_enm_bonds > 0 and self.enm_list.size > 0:
            valid_enm_bonds = []
            for i in range(self.num_enm_bonds):
                atom1, atom2 = self.enm_list[i, 0] - 1, self.enm_list[i, 1] - 1  # 1-based to 0-based
                if atom1 in old_to_new and atom2 in old_to_new:
                    valid_enm_bonds.append([old_to_new[atom1] + 1, old_to_new[atom2] + 1])  # 0-based to 1-based
            
            if valid_enm_bonds:
                new_mol.num_enm_bonds = len(valid_enm_bonds)
                new_mol.enm_list = np.array(valid_enm_bonds, dtype=np.int64)
            else:
                new_mol.num_enm_bonds = 0
                new_mol.enm_list = np.array([], dtype=np.int64).reshape(0, 2)
        
        # Update angle information
        if self.num_angles > 0 and self.angl_list.size > 0:
            valid_angles = []
            for i in range(self.num_angles):
                atom1, atom2, atom3 = self.angl_list[i, 0] - 1, self.angl_list[i, 1] - 1, self.angl_list[i, 2] - 1
                if atom1 in old_to_new and atom2 in old_to_new and atom3 in old_to_new:
                    valid_angles.append([old_to_new[atom1] + 1, old_to_new[atom2] + 1, old_to_new[atom3] + 1])
            
            if valid_angles:
                new_mol.num_angles = len(valid_angles)
                new_mol.angl_list = np.array(valid_angles, dtype=np.int64)
            else:
                new_mol.num_angles = 0
                new_mol.angl_list = np.array([], dtype=np.int64).reshape(0, 3)
        
        # Update dihedral information
        if self.num_dihedrals > 0 and self.dihe_list.size > 0:
            valid_dihedrals = []
            for i in range(self.num_dihedrals):
                atom1, atom2, atom3, atom4 = (self.dihe_list[i, 0] - 1, self.dihe_list[i, 1] - 1, 
                                              self.dihe_list[i, 2] - 1, self.dihe_list[i, 3] - 1)
                if (atom1 in old_to_new and atom2 in old_to_new and 
                    atom3 in old_to_new and atom4 in old_to_new):
                    valid_dihedrals.append([old_to_new[atom1] + 1, old_to_new[atom2] + 1, 
                                            old_to_new[atom3] + 1, old_to_new[atom4] + 1])
            
            if valid_dihedrals:
                new_mol.num_dihedrals = len(valid_dihedrals)
                new_mol.dihe_list = np.array(valid_dihedrals, dtype=np.int64)
            else:
                new_mol.num_dihedrals = 0
                new_mol.dihe_list = np.array([], dtype=np.int64).reshape(0, 4)
        
        # Update improper dihedral information
        if self.num_impropers > 0 and self.impr_list.size > 0:
            valid_impropers = []
            for i in range(self.num_impropers):
                atom1, atom2, atom3, atom4 = (self.impr_list[i, 0] - 1, self.impr_list[i, 1] - 1, 
                                              self.impr_list[i, 2] - 1, self.impr_list[i, 3] - 1)
                if (atom1 in old_to_new and atom2 in old_to_new and 
                    atom3 in old_to_new and atom4 in old_to_new):
                    valid_impropers.append([old_to_new[atom1] + 1, old_to_new[atom2] + 1, 
                                            old_to_new[atom3] + 1, old_to_new[atom4] + 1])
            
            if valid_impropers:
                new_mol.num_impropers = len(valid_impropers)
                new_mol.impr_list = np.array(valid_impropers, dtype=np.int64)
            else:
                new_mol.num_impropers = 0
                new_mol.impr_list = np.array([], dtype=np.int64).reshape(0, 4)
        
        # Update CMAP information
        if self.num_cmaps > 0 and self.cmap_list.size > 0:
            valid_cmaps = []
            for i in range(self.num_cmaps):
                atoms = [self.cmap_list[i, j] - 1 for j in range(8)]  # 1-based to 0-based
                if all(atom in old_to_new for atom in atoms):
                    valid_cmaps.append([old_to_new[atom] + 1 for atom in atoms])  # 0-based to 1-based
            
            if valid_cmaps:
                new_mol.num_cmaps = len(valid_cmaps)
                new_mol.cmap_list = np.array(valid_cmaps, dtype=np.int64)
            else:
                new_mol.num_cmaps = 0
                new_mol.cmap_list = np.array([], dtype=np.int64).reshape(0, 8)
        
        # Update residue and segment information
        unique_residues = np.unique(new_mol.residue_no)
        new_mol.num_residues = len(unique_residues)
        
        unique_segments = np.unique(new_mol.segment_name)
        new_mol.num_segments = len(unique_segments)
        
        # Update molecule information
        unique_molecules = np.unique(new_mol.molecule_no)
        new_mol.num_molecules = len(unique_molecules)
        
        # Update detailed molecule information
        if self.num_molecules > 0 and self.molecule_atom_no.size > 0:
            new_mol.molecule_atom_no = np.array([new_mol.num_atoms], dtype=np.int64)
            new_mol.molecule_mass = np.array([np.sum(new_mol.mass)], dtype=np.float64)
            new_mol.molecule_name = np.array(['SUBSET'], dtype=np.object_)
        
        # Update total charge
        new_mol.total_charge = np.sum(new_mol.charge)
        
        # Update degrees of freedom
        new_mol.num_deg_freedom = new_mol.num_atoms * 3
        
        # Update PC mode information (if applicable)
        new_mol.num_pc_modes = self.num_pc_modes
        if self.pc_mode.size > 0:
            new_mol.pc_mode = self.pc_mode.copy()
        
        # Update FEP information
        new_mol.fep_topology = self.fep_topology
        new_mol.num_hbonds_singleA = self.num_hbonds_singleA
        new_mol.num_hbonds_singleB = self.num_hbonds_singleB
        
        # Update FEP-related arrays
        if hasattr(self.num_atoms_fep, 'copy'):
            new_mol.num_atoms_fep = self.num_atoms_fep.copy()
        else:
            new_mol.num_atoms_fep = self.num_atoms_fep
        
        if hasattr(self.num_bonds_fep, 'copy'):
            new_mol.num_bonds_fep = self.num_bonds_fep.copy()
        else:
            new_mol.num_bonds_fep = self.num_bonds_fep
            
        if hasattr(self.num_angles_fep, 'copy'):
            new_mol.num_angles_fep = self.num_angles_fep.copy()
        else:
            new_mol.num_angles_fep = self.num_angles_fep
            
        if hasattr(self.num_dihedrals_fep, 'copy'):
            new_mol.num_dihedrals_fep = self.num_dihedrals_fep.copy()
        else:
            new_mol.num_dihedrals_fep = self.num_dihedrals_fep
            
        if hasattr(self.num_impropers_fep, 'copy'):
            new_mol.num_impropers_fep = self.num_impropers_fep.copy()
        else:
            new_mol.num_impropers_fep = self.num_impropers_fep
            
        if hasattr(self.num_cmaps_fep, 'copy'):
            new_mol.num_cmaps_fep = self.num_cmaps_fep.copy()
        else:
            new_mol.num_cmaps_fep = self.num_cmaps_fep
        
        # FEP bond/angle/dihedral/improper/cmap arrays are complex, so initialize as empty arrays for now
        if self.bond_list_fep.size > 0:
            new_mol.bond_list_fep = np.array([], dtype=np.int64).reshape(0, 0, 2)
        if self.angl_list_fep.size > 0:
            new_mol.angl_list_fep = np.array([], dtype=np.int64).reshape(0, 0, 3)
        if self.dihe_list_fep.size > 0:
            new_mol.dihe_list_fep = np.array([], dtype=np.int64).reshape(0, 0, 4)
        if self.impr_list_fep.size > 0:
            new_mol.impr_list_fep = np.array([], dtype=np.int64).reshape(0, 0, 4)
        if self.cmap_list_fep.size > 0:
            new_mol.cmap_list_fep = np.array([], dtype=np.int64).reshape(0, 0, 8)
        
        # Other FEP-related arrays
        new_mol.id_singleA = np.array([], dtype=np.int64)
        new_mol.id_singleB = np.array([], dtype=np.int64)
        new_mol.fepgrp = np.array([], dtype=np.int64)
        
        if hasattr(self.fepgrp_bond, 'copy'):
            new_mol.fepgrp_bond = self.fepgrp_bond.copy()
        else:
            new_mol.fepgrp_bond = self.fepgrp_bond
            
        if hasattr(self.fepgrp_angl, 'copy'):
            new_mol.fepgrp_angl = self.fepgrp_angl.copy()
        else:
            new_mol.fepgrp_angl = self.fepgrp_angl
            
        if hasattr(self.fepgrp_dihe, 'copy'):
            new_mol.fepgrp_dihe = self.fepgrp_dihe.copy()
        else:
            new_mol.fepgrp_dihe = self.fepgrp_dihe
            
        if hasattr(self.fepgrp_cmap, 'copy'):
            new_mol.fepgrp_cmap = self.fepgrp_cmap.copy()
        else:
            new_mol.fepgrp_cmap = self.fepgrp_cmap
        
        return new_mol

    def save_psf(self, file_name: str) -> None:
        """
        Save PSF file.
        
        Parameters
        ----------
        file_name : str
            Path to the PSF file to save
        """
        with open(file_name, 'w') as f:
            # PSF header
            f.write("PSF ")
            if self.num_cmaps > 0:
                f.write("CMAP")
            f.write("\n\n")
            
            # Title section
            f.write(f"{1:8d} !NTITLE\n")
            f.write("CREATED by GENESIS Python interface\n")
            f.write("\n")
            
            # Atom section
            f.write(f"{self.num_atoms:8d} !NATOM\n")
            for i in range(self.num_atoms):
                # Atom ID
                f.write(f"{self.atom_no[i]:8d}")
                
                # Segment name (used as chain name)
                if self.segment_name.size > 0 and i < len(self.segment_name):
                    segment_name = str(self.segment_name[i])[:4].ljust(4)
                else:
                    segment_name = "NONE"
                f.write(f" {segment_name}")
                
                # Residue number
                if self.residue_no.size > 0 and i < len(self.residue_no):
                    f.write(f" {self.residue_no[i]:-4d}")
                else:
                    f.write(" 0")
                
                # Residue name
                if self.residue_name.size > 0 and i < len(self.residue_name):
                    residue_name = str(self.residue_name[i])[:4].ljust(4)
                else:
                    residue_name = "NONE"
                f.write(f" {residue_name}")
                
                # Atom name
                if self.atom_name.size > 0 and i < len(self.atom_name):
                    atom_name = str(self.atom_name[i])[:4].ljust(4)
                else:
                    atom_name = "NONE"
                f.write(f" {atom_name}")
                
                # Atom type (using atom class name)
                if self.atom_cls_name.size > 0 and i < len(self.atom_cls_name):
                    atom_type = str(self.atom_cls_name[i])[:4].ljust(4)
                else:
                    atom_type = "NONE"
                f.write(f" {atom_type}")
                
                # Charge
                if self.charge.size > 0 and i < len(self.charge):
                    f.write(f"{self.charge[i]:14.9f}")
                else:
                    f.write("     0.000000000")
                
                # Mass
                if self.mass.size > 0 and i < len(self.mass):
                    f.write(f"{self.mass[i]:14.7f}")
                else:
                    f.write("     0.0000000")
                
                # Reserved field
                f.write("        0")
                f.write("\n")
            
            f.write("\n")
            
            # Bond section
            if self.num_bonds > 0 and self.bond_list.size > 0:
                f.write(f"{self.num_bonds:8d} !NBOND\n")
                for i in range(self.num_bonds):
                    for j in range(2):
                        f.write(f"{self.bond_list[i, j]:8d}")
                    if ((i + 1) % 4 == 0) or (i == self.num_bonds - 1):
                        f.write("\n")
                f.write("\n")
            else:
                f.write("       0 !NBOND\n\n")
            
            # Angle section
            if self.num_angles > 0 and self.angl_list.size > 0:
                f.write(f"{self.num_angles:8d} !NTHETA\n")
                for i in range(self.num_angles):
                    for j in range(3):
                        f.write(f"{self.angl_list[i, j]:8d}")
                    if ((i + 1) % 3 == 0) or (i == self.num_angles - 1):
                        f.write("\n")
                f.write("\n")
            else:
                f.write("       0 !NTHETA\n\n")
            
            # Dihedral section
            if self.num_dihedrals > 0 and self.dihe_list.size > 0:
                f.write(f"{self.num_dihedrals:8d} !NPHI\n")
                for i in range(self.num_dihedrals):
                    for j in range(4):
                        f.write(f"{self.dihe_list[i, j]:8d}")
                    if ((i + 1) % 2 == 0) or (i == self.num_dihedrals - 1):
                        f.write("\n")
                f.write("\n")
            else:
                f.write("       0 !NPHI\n\n")
            
            # Improper dihedral section
            if self.num_impropers > 0 and self.impr_list.size > 0:
                f.write(f"{self.num_impropers:8d} !NIMPHI\n")
                for i in range(self.num_impropers):
                    for j in range(4):
                        f.write(f"{self.impr_list[i, j]:8d}")
                    if ((i + 1) % 2 == 0) or (i == self.num_impropers - 1):
                        f.write("\n")
                f.write("\n")
            else:
                f.write("       0 !NIMPHI\n\n")
            
            # Donor section (empty)
            f.write("       0 !NDON: donors\n\n")
            
            # Acceptor section (empty)
            f.write("       0 !NACC\n\n")
            
            # NNB section (empty)
            f.write("       0 !NNB\n\n")
            
            # CMAP section
            if self.num_cmaps > 0 and self.cmap_list.size > 0:
                f.write(f"{self.num_cmaps:8d} !NCRTERM\n")
                for i in range(self.num_cmaps):
                    for j in range(8):
                        f.write(f"{self.cmap_list[i, j]:8d}")
                    f.write("\n")
                f.write("\n")
            else:
                f.write("       0 !NCRTERM\n\n")

    @staticmethod
    def guess_atom_element(atom_name: str) -> Optional[str]:
        """
        Guess the element from an atom name.

        Parameters
        ----------
        atom_name : atom name

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
                        segment_id=self.segment_name[i])
            atom_name = "".join(self.atom_name[i])
            element = guess_atom_element_mdtraj(atom_name)
            topology.add_atom(atom_name, element, md_residue)
            atom_no_to_idx[self.atom_no[i]] = i
        # create bonds
        for i in range(0, self.num_bonds):
            atom1 = topology.atom(atom_no_to_idx[self.bond_list[i, 0]])
            atom2 = topology.atom(atom_no_to_idx[self.bond_list[i, 1]])
            topology.add_bond(atom1, atom2)
        return topology

    SMolecule.to_mdtraj_topology = to_mdtraj_topology

    @staticmethod
    def from_mdtraj_topology(src: md.Topology) -> Self:
        mol = from_mdtraj_topology_via_pdb(src)
        mol.add_topology_info_from_mdtraj_topology(src)
        return mol

    SMolecule.from_mdtraj_topology = from_mdtraj_topology

    @staticmethod
    def from_mdtraj_trajectory(src: md.Trajectory) -> Self:
        mol = from_mdtraj_trajectory_via_pdb(src)
        mol.add_topology_info_from_mdtraj_topology(src.topology)
        return mol

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
    def guess_atom_element_mdtraj(atom_name) -> Optional[md.element.Element]:
        """
        Guess the MDTraj element from an atom name.

        Parameters
        ----------
        atom_name : atom name

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

    def add_topology_info_from_mdtraj_topology(self, top: md.Topology) -> None:
        # add bonds
        bond_set = {(b[0], b[1]) for b in self.bond_list}
        for bond in top.bonds:
            bond_set.add((bond.atom1.index, bond.atom2.index))
        self.num_bonds = len(bond_set)
        self.bond_list = np.array(list(bond_set), dtype=np.int64)

    SMolecule.add_topology_info_from_mdtraj_topology \
            = add_topology_info_from_mdtraj_topology

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
        # add bonds
        bonds = set()
        for bond in self.bond_list:
            bonds.add((int(bond[0]), int(bond[1])))
        if len(bonds) > 0:
            bonds = tuple(bonds)
            top.add_TopologyAttr(mdaattr.Bonds(bonds))
        return top

    SMolecule.to_mdanalysis_topology = to_mdanalysis_topology

    def to_mdanalysis_universe(self) -> mda.Universe:
        top = self.to_mdanalysis_topology()
        uni = mda.Universe(top)
        return uni

    SMolecule.to_mdanalysis_universe = to_mdanalysis_universe

    @staticmethod
    def from_mdanalysis_universe(uni: mda.Universe) -> Self:
        mol = from_mdanalysis_universe_via_pdb(uni)
        mol.add_topology_info_from_mdanalysis_universe(uni)
        return mol

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

    def add_topology_info_from_mdanalysis_universe(
            self, uni: mda.Universe) -> None:
        # add bonds
        try:
            bond_set = {(b[0], b[1]) for b in self.bond_list}
            for b in uni.bonds:
                bond_set.add((b.atoms[0].index, b.atoms[1].index))
            self.num_bonds = len(bond_set)
            self.bond_list = np.array(list(bond_set), dtype=np.int64)
        except mda.exceptions.NoDataError:
            pass
        try:
            angl_set = {(a[0], a[1], a[2]) for a in self.angl_list}
            for a in uni.angles:
                angl_set.add((a.atoms[0].index,
                              a.atoms[1].index,
                              a.atoms[2].index))
            self.num_angles = len(angl_set)
            self.angl_list = np.array(list(angl_set), dtype=np.int64)
        except mda.exceptions.NoDataError:
            pass

    SMolecule.add_topology_info_from_mdanalysis_universe \
            = add_topology_info_from_mdanalysis_universe

except ImportError:
    pass
