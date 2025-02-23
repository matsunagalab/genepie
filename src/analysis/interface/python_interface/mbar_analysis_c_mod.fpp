!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  ma_main
!! @brief   MBAR analysis tool
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

 module mbar_analysis_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use mbar_analysis_analyze_c_mod

  use mbar_control_mod
  use mbar_option_str_mod
  use output_str_mod
  use input_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

 contains
  ! subroutine mbar_analysis_c(molecule, s_trajes_c, ana_period, ctrl_path) &
  subroutine mbar_analysis_c(ctrl_path, result_fene) &
        bind(C, name="mbar_analysis_c")
    use conv_f_c_util
    implicit none
    ! type(s_molecule_c), intent(in) :: molecule
    ! type(s_trajectories_c), intent(in) :: s_trajes_c
    ! integer, intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_path(*)
    type(c_ptr), intent(out) :: result_fene

    ! type(s_molecule) :: f_molecule
    character(len=:), allocatable :: fort_ctrl_path
    real(wp), pointer :: fene_f(:,:) => null()

    call c2f_string_allocate(ctrl_path, fort_ctrl_path)
    ! call c2f_s_molecule(molecule, f_molecule)
    ! call wa_analysis_main( &
    !     f_molecule, s_trajes_c, ana_period, fort_ctrl_path)
    call mbar_analysis_main( &
        fort_ctrl_path, fene_f)
    result_fene = c_loc(fene_f)
  end subroutine mbar_analysis_c

  ! subroutine mbar_analysis_main( &
  !         molecule, s_trajes_c, ana_period, ctrl_filename)
  subroutine mbar_analysis_main( &
          ctrl_filename, result_fene)
    implicit none
    ! type(s_molecule), intent(inout) :: molecule
    ! type(s_trajectories_c), intent(in) :: s_trajes_c
    ! integer,                intent(in) :: ana_period
    character(*), intent(in) :: ctrl_filename
    real(wp), pointer, intent(out) :: result_fene(:,:)

    ! local variables
    type(s_ctrl_data)      :: ctrl_data
    type(s_molecule)       :: molecule
    type(s_option)         :: option
    type(s_input)          :: input
    type(s_output)         :: output


    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.


    ! [Step1] Read control file
    !
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '

    call control(ctrl_filename, ctrl_data)


    ! [Step2] Set relevant variables and structures 
    !
    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '

    call setup(ctrl_data, molecule, option, input, output)


    ! [Step3] Analyze trajectory
    !
    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '

    ! call analyze(molecule, s_trajes_c, ana_period, input, output, option)
    call analyze(molecule, input, output, option, result_fene)


    ! [Step4] Deallocate memory
    !
    write(MsgOut,'(A)') '[STEP4] Deallocate memory'
    write(MsgOut,'(A)') ' '

    call dealloc_option(option)
    call dealloc_molecules_all(molecule)
end subroutine mbar_analysis_main

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in MBAR_ANALYSIS
  !! @authors      NT
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] molecule   : molecule information
  !! @param[inout] option     : option information
  !! @param[inout] input      : input information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data, molecule, option, input, output)
    use mbar_control_mod
    use mbar_option_mod
    use mbar_option_str_mod
    use output_mod
    use input_mod
    use output_str_mod
    use input_str_mod
    use select_mod
    use molecules_mod
    use molecules_str_mod
    use fileio_grocrd_mod
    use fileio_grotop_mod
    use fileio_ambcrd_mod
    use fileio_prmtop_mod
    use fileio_psf_mod
    use fileio_pdb_mod
    use constants_mod
    implicit none


    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_option),          intent(inout) :: option
    type(s_input),           intent(inout) :: input
    type(s_output),          intent(inout) :: output


    ! setup input
    !
    call setup_input(ctrl_data%inp_info, input)


    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)


    ! setup option
    !
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, option)


    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)


    return

  end subroutine setup

end module mbar_analysis_c_mod
