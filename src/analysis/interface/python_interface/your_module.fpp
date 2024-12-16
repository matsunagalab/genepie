module your_module
  use, intrinsic :: iso_c_binding
  use input_mod
  use molecules_mod
  use fileio_pdb_mod
  use constants_mod      ! Likely contains MaxFilename
  use molecules_str_mod  ! Contains s_molecule definition
  implicit none

  private
  public :: define_molecule_from_pdb

  ! Define MaxFilename if it"s not available from constants_mod
  integer, parameter :: MaxFilename = 256  ! Adjust this value as needed

  type, public, bind(C) :: s_molecule_c
    integer(c_int) :: num_deg_freedom

    integer(c_int) :: num_atoms
    integer(c_int) :: num_bonds
    integer(c_int) :: num_enm_bonds
    integer(c_int) :: num_angles
    integer(c_int) :: num_dihedrals
    integer(c_int) :: num_impropers
    integer(c_int) :: num_cmaps

    integer(c_int) :: num_residues
    integer(c_int) :: num_molecules
    integer(c_int) :: num_segments
    ! shift or not shift
    logical(c_bool) :: shift_origin
    logical(c_bool) :: special_hydrogen

    real(c_double)  :: total_charge
    ! atom
    type(c_ptr) :: atom_no ! c_int[num_atoms]
    type(c_ptr) :: segment_name ! c_char[num_atoms][4]
    type(c_ptr) :: segment_no ! c_int[num_atoms]
    type(c_ptr) :: residue_no ! c_int[num_atoms]
    type(c_ptr) :: residue_c_no ! c_int[num_atoms]
    type(c_ptr) :: residue_name ! c_char[num_atoms][6]
    type(c_ptr) :: atom_name ! c_char[num_atoms][4]
    type(c_ptr) :: atom_cls_name ! c_char[num_atoms][6]
    type(c_ptr) :: atom_cls_no ! c_int[num_atoms]
    type(c_ptr) :: charge ! c_double[num_atoms]
    type(c_ptr) :: mass ! c_double[num_atoms]
    type(c_ptr) :: inv_mass ! c_double[num_atoms]
    type(c_ptr) :: imove ! c_int[num_atoms]
    type(c_ptr) :: stokes_radius ! c_double[num_atoms]
    type(c_ptr) :: inv_stokes_radius ! c_double[num_atoms]

    type(c_ptr) :: chain_id ! c_char[num_atoms][1]
    type(c_ptr) :: atom_coord ! c_double[num_atoms][3]
    type(c_ptr) :: atom_occupancy ! c_double[num_atoms]
    type(c_ptr) :: atom_temp_factor ! c_double[num_atoms]
    type(c_ptr) :: atom_velocity ! c_double[num_atoms][3]
    type(c_ptr) :: light_atom_name ! c_bool[num_atoms]
    type(c_ptr) :: light_atom_mass ! c_bool[num_atoms]

    type(c_ptr) :: molecule_no ! c_int[num_atoms]

    ! bond
    type(c_ptr) :: bond_list ! c_int[num_bonds][2]
    ! enm
    type(c_ptr) :: enm_list ! c_int[num_enm_bonds][2]
    ! angle
    type(c_ptr) :: angl_list ! c_int[num_angles][3]
    ! dihedral
    type(c_ptr) :: dihe_list ! c_int[num_dihedrals][4]
    ! improper
    type(c_ptr) :: impr_list ! c_int[num_impropers][4]
    ! cmap
    type(c_ptr) :: cmap_list ! c_int[num_cmaps][8]
    ! molecule
    type(c_ptr) :: molecule_atom_no ! c_int[num_molecules]
    type(c_ptr) :: molecule_mass ! c_double[num_molecules]
    type(c_ptr) :: molecule_name ! c_char[num_molecules][10]
    ! ref coord
    type(c_ptr) :: atom_refcoord ! c_double[][3]
    ! fit coord
    type(c_ptr) :: atom_fitcoord ! c_double[][3]
    ! principal component mode
    integer(c_int) :: num_pc_modes
    type(c_ptr) :: pc_mode ! c_double[]
    ! ! FEP
    integer(c_int)           :: fep_topology
    integer(c_int)           :: num_hbonds_singleA
    integer(c_int)           :: num_hbonds_singleB
    type(c_ptr)              :: num_atoms_fep     ! c_int[5]
    type(c_ptr)              :: num_bonds_fep     ! c_int[5]
    type(c_ptr)              :: num_angles_fep    ! c_int[5]
    type(c_ptr)              :: num_dihedrals_fep ! c_int[5]
    type(c_ptr)              :: num_impropers_fep ! c_int[5]
    type(c_ptr)              :: num_cmaps_fep     ! c_int[5]
    type(c_ptr) :: bond_list_fep ! c_int[][][]
    type(c_ptr) :: angl_list_fep ! c_int[][][]
    type(c_ptr) :: dihe_list_fep ! c_int[][][]
    type(c_ptr) :: impr_list_fep ! c_int[][][]
    type(c_ptr) :: cmap_list_fep ! c_int[][][]
    type(c_ptr) :: id_singleA ! c_int[]
    type(c_ptr) :: id_singleB ! c_int[]
    type(c_ptr) :: fepgrp ! c_int[]
    type(c_ptr) :: fepgrp_bond ! c_int[5][5]
    type(c_ptr) :: fepgrp_angl ! c_int[5][5][5]
    type(c_ptr) :: fepgrp_dihe ! c_int[5][5][5][5]
    type(c_ptr) :: fepgrp_cmap ! c_int[5*5*5*5*5*5*5*5]
  end type s_molecule_c

contains

  subroutine define_molecule_from_pdb( &
      pdb_filename, out_mol, num_atoms, atom_names, atom_coords) &
      bind(C, name="define_molecule_from_pdb")
    use conv_f_c_util
    implicit none
    ! Input parameters
    character(kind=c_char), intent(in) :: pdb_filename(*)
    ! Output parameters
    type(s_molecule_c), intent(out) :: out_mol
    integer(c_int), intent(out) :: num_atoms
    type(c_ptr), intent(out) :: atom_names
    type(c_ptr), intent(out) :: atom_coords
    ! Local variables
    type(s_inp_info) :: inp_info
    type(s_pdb) :: pdb
    type(s_molecule) :: molecule
    character(MaxFilename) :: filename
    integer :: i
    character(4), pointer :: temp_atom_names(:)
    real(c_float), pointer :: temp_atom_coords(:,:)

    ! Convert C string to Fortran string
    call c2f_string(pdb_filename, filename)
    ! Set up input info
    inp_info%pdbfile = trim(filename)
    ! Read PDB file
    call input_files(inp_info, pdb=pdb)
    ! Define molecule
    call define_molecules(molecule, pdb=pdb)
    ! Set output values
    num_atoms = molecule%num_atoms

    ! Allocate temporary arrays
    allocate(temp_atom_names(num_atoms))
    allocate(temp_atom_coords(3, num_atoms))
    ! Copy data to temporary arrays
    do i = 1, num_atoms
      temp_atom_names(i) = molecule%atom_name(i)
      temp_atom_coords(:,i) = molecule%atom_coord(:,i)
    end do
    ! Convert Fortran arrays to C pointers
    atom_names = c_loc(temp_atom_names)
    atom_coords = c_loc(temp_atom_coords)

    call f2c_s_molecule(molecule, out_mol)
  end subroutine define_molecule_from_pdb

  subroutine f2c_s_molecule(f_src, c_dst)
    use conv_f_c_util
    implicit none
    type(s_molecule), intent(in) :: f_src
    type(s_molecule_c), intent(out) :: c_dst
    integer, pointer :: ptr_int(:)
    integer :: i

    c_dst%num_deg_freedom = f_src%num_deg_freedom
    c_dst%num_atoms = f_src%num_atoms
    c_dst%num_bonds = f_src%num_bonds
    c_dst%num_enm_bonds = f_src%num_enm_bonds
    c_dst%num_angles = f_src%num_angles
    c_dst%num_dihedrals = f_src%num_dihedrals
    c_dst%num_impropers = f_src%num_impropers
    c_dst%num_cmaps     = f_src%num_cmaps

    c_dst%num_residues = f_src%num_residues
    c_dst%num_molecules = f_src%num_molecules
    c_dst%num_segments = f_src%num_segments

    c_dst%shift_origin = f_src%shift_origin
    c_dst%special_hydrogen = f_src%special_hydrogen

    c_dst%total_charge = f_src%total_charge

    c_dst%atom_no = f2c_int_array_allocatable(f_src%atom_no)
    c_dst%segment_name  = f2c_string_array_allocatable(f_src%segment_name)
    c_dst%segment_no    = f2c_int_array_allocatable(f_src%segment_no)
    c_dst%residue_no    = f2c_int_array_allocatable(f_src%residue_no)
    c_dst%residue_c_no  = f2c_int_array_allocatable(f_src%residue_c_no)
    c_dst%residue_name  = f2c_string_array_allocatable(f_src%residue_name)
    c_dst%atom_name     = f2c_string_array_allocatable(f_src%atom_name)
    c_dst%atom_cls_name = f2c_string_array_allocatable(f_src%atom_cls_name)
    c_dst%atom_cls_no   = f2c_int_array_allocatable(f_src%atom_cls_no)
    c_dst%charge            = f2c_double_array_allocatable(f_src%charge)
    c_dst%mass              = f2c_double_array_allocatable(f_src%mass)
    c_dst%inv_mass          = f2c_double_array_allocatable(f_src%inv_mass)
    c_dst%imove             = f2c_int_array_allocatable(f_src%imove)
    c_dst%stokes_radius     = f2c_double_array_allocatable(f_src%stokes_radius)
    c_dst%inv_stokes_radius = f2c_double_array_allocatable(f_src%inv_stokes_radius)

    c_dst%chain_id = f2c_string_array_allocatable(f_src%chain_id )
    c_dst%atom_coord = f2c_double_array_2dim_allocatable(f_src%atom_coord)
    c_dst%atom_occupancy = f2c_double_array_allocatable(f_src%atom_occupancy)
    c_dst%atom_temp_factor = f2c_double_array_allocatable(f_src%atom_temp_factor)
    c_dst%atom_velocity = f2c_double_array_2dim_allocatable(f_src%atom_velocity)
    c_dst%light_atom_name = f2c_bool_array_allocatable(f_src%light_atom_name)
    c_dst%light_atom_mass = f2c_bool_array_allocatable(f_src%light_atom_mass)

    c_dst%molecule_no = f2c_int_array_allocatable(f_src%molecule_no)

    c_dst%bond_list = f2c_int_array_2dim_allocatable(f_src%bond_list)
    c_dst%enm_list = f2c_int_array_2dim_allocatable(f_src%enm_list)
    c_dst%angl_list = f2c_int_array_2dim_allocatable(f_src%angl_list)
    c_dst%dihe_list = f2c_int_array_2dim_allocatable(f_src%dihe_list)
    c_dst%impr_list = f2c_int_array_2dim_allocatable(f_src%impr_list)
    c_dst%cmap_list = f2c_int_array_2dim_allocatable(f_src%cmap_list)
    c_dst%molecule_atom_no = f2c_int_array_allocatable(f_src%molecule_atom_no)
    c_dst%molecule_mass = f2c_double_array_allocatable(f_src%molecule_mass)
    c_dst%molecule_name = f2c_string_array_allocatable(f_src%molecule_name)
    c_dst%atom_refcoord = f2c_double_array_2dim_allocatable(f_src%atom_refcoord)
    c_dst%atom_fitcoord = f2c_double_array_2dim_allocatable(f_src%atom_fitcoord)
    c_dst%num_pc_modes = f_src%num_pc_modes
    c_dst%pc_mode = f2c_double_array_allocatable(f_src%pc_mode)

    c_dst%fep_topology       = f_src%fep_topology
    c_dst%num_hbonds_singleA = f_src%num_hbonds_singleA
    c_dst%num_hbonds_singleB = f_src%num_hbonds_singleB

    c_dst%num_atoms_fep     = f2c_int_array(f_src%num_atoms_fep)
    c_dst%num_bonds_fep     = f2c_int_array(f_src%num_bonds_fep)
    c_dst%num_angles_fep    = f2c_int_array(f_src%num_angles_fep)
    c_dst%num_dihedrals_fep = f2c_int_array(f_src%num_dihedrals_fep)
    c_dst%num_impropers_fep = f2c_int_array(f_src%num_impropers_fep)
    c_dst%num_cmaps_fep     = f2c_int_array(f_src%num_cmaps_fep)

    c_dst%bond_list_fep = f2c_int_array_3dim_allocatable(f_src%bond_list_fep)
    c_dst%angl_list_fep = f2c_int_array_3dim_allocatable(f_src%angl_list_fep)
    c_dst%dihe_list_fep = f2c_int_array_3dim_allocatable(f_src%dihe_list_fep)
    c_dst%impr_list_fep = f2c_int_array_3dim_allocatable(f_src%impr_list_fep)
    c_dst%cmap_list_fep = f2c_int_array_3dim_allocatable(f_src%cmap_list_fep)

    c_dst%id_singleA  = f2c_int_array_allocatable(f_src%id_singleA)
    c_dst%id_singleB  = f2c_int_array_allocatable(f_src%id_singleB)
    c_dst%fepgrp      = f2c_int_array_allocatable(f_src%fepgrp)
    c_dst%fepgrp_bond = f2c_int_array_2dim(f_src%fepgrp_bond)
    c_dst%fepgrp_angl = f2c_int_array_3dim(f_src%fepgrp_angl)
    c_dst%fepgrp_dihe = f2c_int_array_4dim(f_src%fepgrp_dihe)
    c_dst%fepgrp_cmap = f2c_int_array(f_src%fepgrp_cmap)
  end subroutine

  subroutine deallocate_s_molecule_c(self) &
      bind(C, name="deallocate_s_molecule_c")
    use, intrinsic :: iso_c_binding
    implicit none
    type(s_molecule_c), intent(inout) :: self
    character(kind=c_char), pointer :: work_string(:,:)
    integer(c_int), pointer :: work_int(:)
    integer(c_int), pointer :: work_int2(:,:)
    real(c_double), pointer :: work_double(:)
    real(c_double), pointer :: work_double2(:,:)

    if (c_associated(self%atom_no)) then
      call C_F_pointer(self%atom_no, work_int, [self%num_atoms])
      deallocate(work_int)
    end if
    if (c_associated(self%segment_name)) then
      call C_F_pointer(self%segment_name, work_int2, [4, self%num_atoms])
      deallocate(work_int2)
    end if
    if (c_associated(self%segment_no)) then
      call C_F_pointer(self%segment_no, work_int, [self%num_atoms])
      deallocate(work_int)
    end if
    if (c_associated(self%residue_no)) then
      call C_F_pointer(self%residue_no, work_int, [self%num_atoms])
      deallocate(work_int)
    end if
    if (c_associated(self%residue_c_no)) then
      call C_F_pointer(self%residue_c_no, work_int, [self%num_atoms])
      deallocate(work_int)
    end if
    if (c_associated(self%residue_name)) then
      call C_F_pointer(self%residue_name, work_string, [6, self%num_atoms])
      deallocate(work_string)
    end if
    if (c_associated(self%atom_name)) then
      call C_F_pointer(self%atom_name, work_string, [4, self%num_atoms])
      deallocate(work_string)
    end if
    if (c_associated(self%atom_cls_name)) then
      call C_F_pointer(self%atom_cls_name, work_string, [6, self%num_atoms])
      deallocate(work_string)
    end if
    if (c_associated(self%atom_cls_no)) then
      call C_F_pointer(self%atom_cls_no, work_int, [self%num_atoms])
      deallocate(work_int)
    end if
    if (c_associated(self%charge)) then
      call C_F_pointer(self%charge, work_double, [self%num_atoms])
      deallocate(work_double)
    end if
    if (c_associated(self%mass)) then
      call C_F_pointer(self%mass, work_double, [self%num_atoms])
      deallocate(work_double)
    end if
    if (c_associated(self%inv_mass)) then
      call C_F_pointer(self%inv_mass, work_double, [self%num_atoms])
      deallocate(work_double)
    end if
    if (c_associated(self%imove)) then
      call C_F_pointer(self%imove, work_int, [self%num_atoms])
      deallocate(work_int)
    end if
    if (c_associated(self%stokes_radius)) then
      call C_F_pointer(self%stokes_radius, work_double, [self%num_atoms])
      deallocate(work_double)
    end if
    if (c_associated(self%inv_stokes_radius)) then
      call C_F_pointer(self%inv_stokes_radius, work_double, [self%num_atoms])
      deallocate(work_double)
    end if
    ! TODO

  !   deallocate(self%chain_id)
  !   deallocate(self%atom_coord)
  !   deallocate(self%atom_occupancy)
  !   deallocate(self%atom_temp_factor)
  !   deallocate(self%atom_velocity)
  !   deallocate(self%light_atom_name)
  !   deallocate(self%light_atom_mass)
  !   deallocate(self%molecule_no)
  !   deallocate(self%bond_list)
  !   deallocate(self%enm_list)
  !   deallocate(self%angl_list)
  !   deallocate(self%dihe_list)
  !   deallocate(self%impr_list)
  !   deallocate(self%cmap_list)
  !   deallocate(self%molecule_atom_no)
  !   deallocate(self%molecule_mass)
  !   deallocate(self%molecule_name)
  !   deallocate(self%atom_refcoord)
  !   deallocate(self%atom_fitcoord)
  !   deallocate(self%pc_mode)
  !   deallocate(self%num_atoms_fep)
  !   deallocate(self%num_bonds_fep)
  !   deallocate(self%num_angles_fep)
  !   deallocate(self%num_dihedrals_fep)
  !   deallocate(self%num_impropers_fep)
  !   deallocate(self%num_cmaps_fep)
  !   deallocate(self%bond_list_fep)
  !   deallocate(self%angl_list_fep)
  !   deallocate(self%dihe_list_fep)
  !   deallocate(self%impr_list_fep)
  !   deallocate(self%cmap_list_fep)
  !   deallocate(self%id_singleA)
  !   deallocate(self%id_singleB)
  !   deallocate(self%fepgrp)
  !   deallocate(self%fepgrp_bond)
  !   deallocate(self%fepgrp_angl)
  !   deallocate(self%fepgrp_dihe)
  !   deallocate(self%fepgrp_cmap)
  end subroutine deallocate_s_molecule_c

end module your_module
