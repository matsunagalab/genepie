!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  ma_main
!! @brief   MSD analysis
!! @authors Donatas Surblys (DS), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module msd_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use msd_impl_mod

  use ma_control_mod
  use ma_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use string_mod
  use messages_mod
  use constants_mod
  implicit none

  public :: ma_analysis_c
  public :: deallocate_msd_results_c

  ! Module-level pointer for results (to be deallocated later)
  real(wp), pointer, save :: msd_ptr(:,:) => null()

contains
   subroutine ma_analysis_c(molecule, s_trajes_c, ana_period, &
                            ctrl_text, ctrl_len, &
                            result_msd, num_analysis_mols, num_delta) &
        bind(C, name="ma_analysis_c")
    use conv_f_c_util
    implicit none
    type(s_molecule_c), intent(in) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer, intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out) :: result_msd
    integer, intent(out) :: num_analysis_mols
    integer, intent(out) :: num_delta

    type(s_molecule) :: f_molecule

    ! Deallocate previous results if any
    if (associated(msd_ptr)) then
      deallocate(msd_ptr)
      nullify(msd_ptr)
    end if

    call c2f_s_molecule(molecule, f_molecule)
    call ma_analysis_main( &
        f_molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, &
        msd_ptr, num_analysis_mols, num_delta)
    result_msd = c_loc(msd_ptr)
  end subroutine ma_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_msd_results_c
  !> @brief        Deallocate MSD analysis results
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine deallocate_msd_results_c() bind(C, name="deallocate_msd_results_c")
    implicit none

    if (associated(msd_ptr)) then
      deallocate(msd_ptr)
      nullify(msd_ptr)
    end if

  end subroutine deallocate_msd_results_c

  subroutine ma_analysis_main( &
          molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, &
          result_msd, num_analysis_mols, num_delta)
    implicit none
    type(s_molecule), intent(inout) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer,                intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer,                intent(in) :: ctrl_len
    real(wp), pointer, intent(out) :: result_msd(:,:)
    integer, intent(out) :: num_analysis_mols
    integer, intent(out) :: num_delta


    ! local variables
    type(s_ctrl_data)      :: ctrl_data
    type(s_trajectory)     :: trajectory
    type(s_output)         :: output
    type(s_option)         :: option

    ! [Step1] Read control file
    !
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '

    call control_from_string(ctrl_text, ctrl_len, ctrl_data)


    ! [Step2] Set relevant variables and structures
    !
    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '

    call setup(molecule, ctrl_data, output, option)


    ! [Step3] Analyze trajectory
    !
    write(MsgOut,'(A)') '[STEP3] Analysis of trajectory files'
    write(MsgOut,'(A)') ' '

    call analyze(molecule, s_trajes_c, ana_period, output, option, &
                 result_msd, num_analysis_mols, num_delta)


    ! [Step4] Deallocate memory
    !
    write(MsgOut,'(A)') '[STEP4] Deallocate memory'
    write(MsgOut,'(A)') ' '

    call dealloc_trajectory(trajectory)
    call dealloc_option(option)
    call dealloc_molecules_all(molecule)
end subroutine ma_analysis_main

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures
  !! @authors      TM
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] output     : output information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(molecule, ctrl_data, output, option)
    use ma_control_mod
    use ma_option_mod
    use ma_option_str_mod
    use trajectory_mod
    use output_mod
    use input_mod
    use trajectory_str_mod
    use output_str_mod
    use select_mod
    use molecules_mod
    use molecules_str_mod
    use fileio_psf_mod
    use fileio_pdb_mod
    use fileio_ambcrd_mod
    use fileio_prmtop_mod
    use select_molecules_mod
    implicit none


    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_output),          intent(inout) :: output
    type(s_option),          intent(inout) :: option

    ! local variables
    type(s_selmols), dimension(:), allocatable      :: selmols
    type(s_one_molecule), dimension(:), allocatable :: allmols


    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)


    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)

    ! setup molecule selection
    !
    call setup_molselect(ctrl_data%molsel_info, ctrl_data%sel_info, &
      molecule, selmols, allmols=allmols)

    ! setup option
    !
    call setup_option(ctrl_data%opt_info, selmols, allmols, molecule, option)


    return

  end subroutine setup

end module msd_c_mod
