!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  ta_main
!! @brief   analysis trajectory files
!! @authors Takaharu Mori (TM), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module trj_analysis_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use trj_analysis_analyze_c_mod

  use ta_control_mod
  use ta_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

contains
  subroutine trj_analysis_c(molecule, s_trajes_c, ana_period, ctrl_path, &
                            result_distance, num_distance, &
                            result_angle, num_angle, &
                            result_torsion, num_torsion, &
                            result_cdis, num_cdis, &
                            result_cang, num_cang, &
                            result_ctor, num_ctor) &
        bind(C, name="trj_analysis_c")
    use conv_f_c_util
    implicit none
    type(s_molecule_c), intent(in) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer, intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_path(*)
    type(c_ptr), intent(out) :: result_distance
    integer(c_int), intent(out) :: num_distance
    type(c_ptr), intent(out) :: result_angle
    integer(c_int), intent(out) :: num_angle
    type(c_ptr), intent(out) :: result_torsion
    integer(c_int), intent(out) :: num_torsion
    type(c_ptr), intent(out) :: result_cdis
    integer(c_int), intent(out) :: num_cdis
    type(c_ptr), intent(out) :: result_cang
    integer(c_int), intent(out) :: num_cang
    type(c_ptr), intent(out) :: result_ctor
    integer(c_int), intent(out) :: num_ctor

    type(s_molecule) :: f_molecule
    character(len=:), allocatable :: fort_ctrl_path
    real(wp), pointer :: distance_f(:,:) => null()
    integer :: num_distance_f
    real(wp), pointer :: angle_f(:,:) => null()
    integer :: num_angle_f
    real(wp), pointer :: torsion_f(:,:) => null()
    integer :: num_torsion_f
    real(wp), pointer :: cdis_f(:,:) => null()
    integer :: num_cdis_f
    real(wp), pointer :: cang_f(:,:) => null()
    integer :: num_cang_f
    real(wp), pointer :: ctor_f(:,:) => null()
    integer :: num_ctor_f

    call c2f_string_allocate(ctrl_path, fort_ctrl_path)
    call c2f_s_molecule(molecule, f_molecule)
    call trj_analysis_main( &
        f_molecule, s_trajes_c, ana_period, fort_ctrl_path, &
        distance_f, num_distance_f, &
        angle_f, num_angle_f, &
        torsion_f, num_torsion_f, &
        cdis_f, num_cdis_f, &
        cang_f, num_cang_f, &
        ctor_f, num_ctor_f)
    call dealloc_molecules_all(f_molecule)

    if (associated(distance_f)) then
      result_distance = c_loc(distance_f)
      num_distance = num_distance_f
    else
      result_distance = c_null_ptr
      num_distance = 0
    end if
    if (associated(angle_f)) then
      result_angle = c_loc(angle_f)
      num_angle = num_angle_f
    else
      result_angle = c_null_ptr
      num_angle = 0
    end if
    if (associated(torsion_f)) then
      result_torsion = c_loc(torsion_f)
      num_torsion = num_torsion_f
    else
      result_torsion = c_null_ptr
      num_torsion = 0
    end if
    if (associated(cdis_f)) then
      result_cdis = c_loc(cdis_f)
      num_cdis = num_cdis_f
    else
      result_cdis = c_null_ptr
      num_cdis = 0
    end if
    if (associated(cang_f)) then
      result_cang = c_loc(cang_f)
      num_cang = num_cang_f
    else
      result_cang = c_null_ptr
      num_cang = 0
    end if
    if (associated(ctor_f)) then
      result_ctor = c_loc(ctor_f)
      num_ctor = num_ctor_f
    else
      result_ctor = c_null_ptr
      num_ctor = 0
    end if
  end subroutine trj_analysis_c

  subroutine trj_analysis_main( &
          molecule, s_trajes_c, ana_period, ctrl_filename, &
          result_distance, num_distance, &
          result_angle, num_angle, &
          result_torsion, num_torsion, &
          result_cdis, num_cdis, &
          result_cang, num_cang, &
          result_ctor, num_ctor)
    implicit none
    type(s_molecule), intent(inout) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer,                intent(in) :: ana_period
    character(*), intent(in) :: ctrl_filename
    real(wp), pointer, intent(out) :: result_distance(:,:)
    integer, intent(out) :: num_distance
    real(wp), pointer, intent(out) :: result_angle(:,:)
    integer, intent(out) :: num_angle
    real(wp), pointer, intent(out) :: result_torsion(:,:)
    integer, intent(out) :: num_torsion
    real(wp), pointer, intent(out) :: result_cdis(:,:)
    integer, intent(out) :: num_cdis
    real(wp), pointer, intent(out) :: result_cang(:,:)
    integer, intent(out) :: num_cang
    real(wp), pointer, intent(out) :: result_ctor(:,:)
    integer, intent(out) :: num_ctor

    ! local variables
    type(s_ctrl_data)      :: ctrl_data
    type(s_trajectory)     :: trajectory
    type(s_output)         :: output
    type(s_option)         :: option


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

    call setup(molecule, ctrl_data, output, option)


    ! [Step3] Analyze trajectory
    !
    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '

    call analyze(molecule, s_trajes_c, ana_period, option, &
                 result_distance, num_distance, &
                 result_angle, num_angle, &
                 result_torsion, num_torsion, &
                 result_cdis, num_cdis, &
                 result_cang, num_cang, &
                 result_ctor, num_ctor)


    ! [Step4] Deallocate memory
    !
    write(MsgOut,'(A)') '[STEP4] Deallocate memory'
    write(MsgOut,'(A)') ' '

    call dealloc_trajectory(trajectory)
  end subroutine trj_analysis_main

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in TRJ_ANALYSIS
  !! @authors      NT, TM
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] output     : output information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(molecule, ctrl_data, output, option)
    use ta_control_mod
    use ta_option_mod
    use ta_option_str_mod
    use trajectory_mod
    use output_mod
    use input_mod
    use trajectory_str_mod
    use output_str_mod
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
    type(s_molecule),        intent(inout) :: molecule
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_option),          intent(inout) :: option

    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)


    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)


    ! setup option
    !
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, output, option)

    return

  end subroutine setup

end module trj_analysis_c_mod
