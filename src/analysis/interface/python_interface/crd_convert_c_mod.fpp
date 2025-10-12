!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  cc_main
!! @brief   convert MD trajectory format
!! @authors Norio Takase (NT), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module crd_convert_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use crd_convert_convert
  use conv_f_c_util

  use cc_control_mod
  use cc_option_str_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  implicit none

contains
  subroutine crd_convert_c( &
          molecule, ctrl_path, s_trajes_c_array, num_trajs, &
          selected_atom_indices, num_selected_atoms) &
          bind(C, name="crd_convert_c")
    implicit none
    type(s_molecule_c), intent(inout) :: molecule
    character(kind=c_char), intent(in) :: ctrl_path(*)
    type(c_ptr), intent(out) :: s_trajes_c_array
    integer(c_int), intent(out) :: num_trajs
    type(c_ptr), intent(out) :: selected_atom_indices
    integer(c_int), intent(out) :: num_selected_atoms
    character(len=:), allocatable :: fort_ctrl_path
    type(s_molecule) :: f_molecule

    call c2f_string_allocate(ctrl_path, fort_ctrl_path)
    call c2f_s_molecule(molecule, f_molecule)

    call crd_convert_main(f_molecule, fort_ctrl_path, s_trajes_c_array, num_trajs, &
                         selected_atom_indices, num_selected_atoms)
  end subroutine crd_convert_c

  subroutine crd_convert_main(molecule, ctrl_filename, s_trajes_c_array, num_trajs, &
                              selected_atom_indices, num_selected_atoms)
    implicit none
    type(s_molecule), intent(inout) :: molecule
    character(*), intent(in) :: ctrl_filename
    type(c_ptr), intent(out) :: s_trajes_c_array
    integer(c_int), intent(out) :: num_trajs
    type(c_ptr), intent(out) :: selected_atom_indices
    integer(c_int), intent(out) :: num_selected_atoms
    type(s_ctrl_data)      :: ctrl_data
    type(s_trj_list)       :: trj_list
    type(s_trajectory)     :: trajectory
    type(s_fitting)        :: fitting
    type(s_option)         :: option
    type(s_output)         :: output


    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.


    ! [Step1] Read control file
    !
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Convert'
    write(MsgOut,'(A)') ' '

    call control(ctrl_filename, ctrl_data)


    ! [Step2] Set relevant variables and structures
    !
    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '

    call setup(ctrl_data,  &
               molecule,   &
               trj_list,   &
               trajectory, &
               fitting,    &
               option,     &
               output)


    ! [Step3] Convert trajectory files
    !
    write(MsgOut,'(A)') '[STEP3] Convert trajectory files'
    write(MsgOut,'(A)') ' '

    call convert(molecule,   &
                 trj_list,   &
                 trajectory, &
                 fitting,    &
                 option,     &
                 output,     &
                 s_trajes_c_array, &
                 num_trajs)

    ! Extract selected atom indices from option%trjout_atom
    write(MsgOut,'(A)') 'About to extract selected atom indices...'
    write(MsgOut,'(A,I8)') 'Extract_Selected_Atom_Indices> Number of selected atoms: ', &
         size(option%trjout_atom%idx)
    if (size(option%trjout_atom%idx) > 0) then
      write(MsgOut,'(A,10I8)') 'Extract_Selected_Atom_Indices> First 10 indices: ', &
           option%trjout_atom%idx(1:min(10, size(option%trjout_atom%idx)))
    end if
    write(MsgOut,'(A)') 'Calling extract_selected_atom_indices...'
    call extract_selected_atom_indices(option%trjout_atom, selected_atom_indices, num_selected_atoms)
    write(MsgOut,'(A)') 'extract_selected_atom_indices completed.'


    ! [Step4] Deallocate memory
    !
    write(MsgOut,'(A)') '[STEP4] Deallocate memory'
    write(MsgOut,'(A)') ' '

    call dealloc_option(option)
    call dealloc_fitting(fitting)
    call dealloc_trajectory(trajectory)
    call dealloc_trj_list(trj_list)
  end subroutine crd_convert_main

  subroutine setup(ctrl_data,  &
                   molecule,   &
                   trj_list,   &
                   trajectory, &
                   fitting,    &
                   option,     &
                   output)
    use cc_control_mod
    use cc_option_mod
    use cc_option_str_mod
    use fitting_mod
    use fitting_str_mod
    use input_mod
    use output_mod
    use output_str_mod
    use trajectory_mod
    use trajectory_str_mod
    use select_mod
    use molecules_mod
    use molecules_str_mod
    use fileio_grocrd_mod
    use fileio_grotop_mod
    use fileio_ambcrd_mod
    use fileio_prmtop_mod
    use fileio_psf_mod
    use fileio_pdb_mod
    implicit none

    ! formal arguments
    type(s_ctrl_data),  intent(in)    :: ctrl_data
    type(s_molecule),   intent(inout) :: molecule
    type(s_trj_list),   intent(inout) :: trj_list
    type(s_trajectory), intent(inout) :: trajectory
    type(s_fitting),    intent(inout) :: fitting
    type(s_option),     intent(inout) :: option
    type(s_output),     intent(inout) :: output

    ! local variables
    type(s_psf)              :: psf
    type(s_pdb)              :: ref, ref_out
    type(s_prmtop)           :: prmtop
    type(s_ambcrd)           :: ambcrd
    type(s_grotop)           :: grotop
    type(s_grocrd)           :: grocrd

    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_prmtop_all(prmtop)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_grotop_all(grotop)
    call dealloc_grocrd_all(grocrd)


    ! setup trajectory
    !
    call setup_trajectory(ctrl_data%trj_info, &
                          molecule, trj_list, trajectory)

    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)

    ! setup fitting
    !
    call setup_fitting(ctrl_data%fit_info, ctrl_data%sel_info, &
                       molecule, fitting)

    ! setup option
    !
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, option)

    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)


    ! export reference molecules
    !
    if (output%pdbfile /= '') then

      call export_molecules(molecule, option%trjout_atom, ref_out)
      call output_pdb(output%pdbfile, ref_out)
      call dealloc_pdb_all(ref_out)

    end if

    return

  end subroutine setup

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    extract_selected_atom_indices
  !> @brief        extract selected atom indices from s_selatoms to C array
  !! @authors      Generated
  !! @param[in]    selatoms              : selected atoms structure
  !! @param[out]   selected_atom_indices : C pointer to integer array
  !! @param[out]   num_selected_atoms    : number of selected atoms
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine extract_selected_atom_indices(selatoms, selected_atom_indices, num_selected_atoms)
    use, intrinsic :: iso_c_binding
    use select_atoms_str_mod
    use conv_f_c_util
    use messages_mod
    implicit none

    ! formal arguments
    type(s_selatoms), intent(in) :: selatoms
    type(c_ptr), intent(out) :: selected_atom_indices
    integer(c_int), intent(out) :: num_selected_atoms

    ! local variables
    integer :: i, nsel
    integer(c_int), pointer :: c_array(:)

    nsel = size(selatoms%idx)
    num_selected_atoms = nsel
    
    write(MsgOut,'(A,I8)') 'Extract_Selected_Atom_Indices> Called with nsel = ', nsel

    if (nsel > 0) then
      write(MsgOut,'(A)') 'Extract_Selected_Atom_Indices> Allocating C array'
      ! Allocate C array using conv_f_c_util
      selected_atom_indices = allocate_c_int_array(int(nsel, c_int))
      call c_f_pointer(selected_atom_indices, c_array, [nsel])
      do i = 1, nsel
        c_array(i) = int(selatoms%idx(i), c_int)
      end do
      write(MsgOut,'(A)') 'Extract_Selected_Atom_Indices> C array allocated and filled'
    else
      write(MsgOut,'(A)') 'Extract_Selected_Atom_Indices> No atoms selected, returning null'
      selected_atom_indices = c_null_ptr
    end if

  end subroutine extract_selected_atom_indices

end module crd_convert_c_mod
