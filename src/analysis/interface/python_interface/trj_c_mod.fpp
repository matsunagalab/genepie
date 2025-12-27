!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> C language interface for trj_analysis
!! @brief   C language interface of trj_analysis
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module trj_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use trj_impl_mod

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

  public :: trj_analysis_c
  public :: deallocate_trj_results_c

  ! Module-level pointers for results (to be deallocated later)
  real(wp), pointer, save :: distance_ptr(:,:) => null()
  real(wp), pointer, save :: angle_ptr(:,:) => null()
  real(wp), pointer, save :: torsion_ptr(:,:) => null()
  real(wp), pointer, save :: cdis_ptr(:,:) => null()
  real(wp), pointer, save :: cang_ptr(:,:) => null()
  real(wp), pointer, save :: ctor_ptr(:,:) => null()

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_c
  !> @brief        C language interface for trj_analysis
  !! @param[in] molecule : structure of molecule information
  !! @param[in] s_trajes_c : structure of trajectories information
  !! @param[in] ana_period : interval for sampling frames
  !! @param[in] ctrl_path : path to the control file
  !! @param[out] result_distance : array storing calculated distances
  !! @param[out] num_distance : number of distance results stored in result_distance
  !! @param[out] result_angle : array storing calculated angles
  !! @param[out] num_angle : number of angle results stored in result_angle
  !! @param[out] result_torsion : array storing calculated torsion angles
  !! @param[out] num_torsion : number of torsion angle results stored in result_torsion
  !! @param[out] result_cdis : array storing calculated contact distances
  !! @param[out] num_cdis : number of contact distance results stored in result_cdis
  !! @param[out] result_cang : array storing calculated contact angles
  !! @param[out] num_cang : number of contact angle results stored in result_cang
  !! @param[out] result_ctor : array storing calculated contact torsion angles
  !! @param[out] num_ctor : number of contact torsion angle results stored in result_ctor
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine trj_analysis_c(molecule, s_trajes_c, ana_period, &
                            ctrl_text, ctrl_len, &
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
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
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

    ! Deallocate previous results if any
    call deallocate_trj_results_internal()

    call c2f_s_molecule(molecule, f_molecule)
    call trj_analysis_main( &
        f_molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, &
        distance_ptr, angle_ptr, torsion_ptr, cdis_ptr, cang_ptr, ctor_ptr)
    call dealloc_molecules_all(f_molecule)

    if (associated(distance_ptr)) then
      result_distance = c_loc(distance_ptr)
      num_distance = size(distance_ptr, dim=1)
    else
      result_distance = c_null_ptr
      num_distance = 0
    end if
    if (associated(angle_ptr)) then
      result_angle = c_loc(angle_ptr)
      num_angle = size(angle_ptr, dim=1)
    else
      result_angle = c_null_ptr
      num_angle = 0
    end if
    if (associated(torsion_ptr)) then
      result_torsion = c_loc(torsion_ptr)
      num_torsion = size(torsion_ptr, dim=1)
    else
      result_torsion = c_null_ptr
      num_torsion = 0
    end if
    if (associated(cdis_ptr)) then
      result_cdis = c_loc(cdis_ptr)
      num_cdis = size(cdis_ptr, dim=1)
    else
      result_cdis = c_null_ptr
      num_cdis = 0
    end if
    if (associated(cang_ptr)) then
      result_cang = c_loc(cang_ptr)
      num_cang = size(cang_ptr, dim=1)
    else
      result_cang = c_null_ptr
      num_cang = 0
    end if
    if (associated(ctor_ptr)) then
      result_ctor = c_loc(ctor_ptr)
      num_ctor = size(ctor_ptr, dim=1)
    else
      result_ctor = c_null_ptr
      num_ctor = 0
    end if
  end subroutine trj_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_trj_results_internal
  !> @brief        Internal helper to deallocate TRJ analysis results
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine deallocate_trj_results_internal()
    implicit none

    if (associated(distance_ptr)) then
      deallocate(distance_ptr)
      nullify(distance_ptr)
    end if
    if (associated(angle_ptr)) then
      deallocate(angle_ptr)
      nullify(angle_ptr)
    end if
    if (associated(torsion_ptr)) then
      deallocate(torsion_ptr)
      nullify(torsion_ptr)
    end if
    if (associated(cdis_ptr)) then
      deallocate(cdis_ptr)
      nullify(cdis_ptr)
    end if
    if (associated(cang_ptr)) then
      deallocate(cang_ptr)
      nullify(cang_ptr)
    end if
    if (associated(ctor_ptr)) then
      deallocate(ctor_ptr)
      nullify(ctor_ptr)
    end if

  end subroutine deallocate_trj_results_internal

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_trj_results_c
  !> @brief        Deallocate TRJ analysis results (C interface)
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine deallocate_trj_results_c() bind(C, name="deallocate_trj_results_c")
    implicit none

    call deallocate_trj_results_internal()

  end subroutine deallocate_trj_results_c

  subroutine trj_analysis_main( &
          molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, &
          result_distance, result_angle, result_torsion, &
          result_cdis, result_cang, result_ctor)
    implicit none
    type(s_molecule), intent(inout) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer,                intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer,                intent(in) :: ctrl_len
    real(wp), pointer, intent(out) :: result_distance(:,:)
    real(wp), pointer, intent(out) :: result_angle(:,:)
    real(wp), pointer, intent(out) :: result_torsion(:,:)
    real(wp), pointer, intent(out) :: result_cdis(:,:)
    real(wp), pointer, intent(out) :: result_cang(:,:)
    real(wp), pointer, intent(out) :: result_ctor(:,:)

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

    call control_from_string(ctrl_text, ctrl_len, ctrl_data)


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
                 result_distance, result_angle, result_torsion, &
                 result_cdis, result_cang, result_ctor)


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

end module trj_c_mod
