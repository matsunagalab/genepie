!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  ra_main
!! @brief   RMSD analysis
!! @authors Takaharu Mori (TM), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rmsd_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use rmsd_impl_mod

  use ra_control_mod
  use ra_option_str_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use error_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: ra_analysis_c
  public :: deallocate_rmsd_results_c

  ! Module-level pointer for results (to be deallocated later)
  real(wp), pointer, save :: ra_ptr(:) => null()

contains
  subroutine ra_analysis_c(molecule, s_trajes_c, ana_period, &
                           ctrl_text, ctrl_len, &
                           result_ra, status, msg, msglen) &
        bind(C, name="ra_analysis_c")
    use conv_f_c_util
    implicit none
    type(s_molecule_c), intent(in) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer, intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out) :: result_ra
    integer(c_int),          intent(out) :: status
    character(kind=c_char),  intent(out) :: msg(*)
    integer(c_int),          value       :: msglen

    type(s_molecule) :: f_molecule

    type(s_error) :: err

    ! Deallocate previous results if any
    if (associated(ra_ptr)) then
      deallocate(ra_ptr)
      nullify(ra_ptr)
    end if

    call error_init(err)
    call c2f_s_molecule(molecule, f_molecule)
    call ra_analysis_main( &
        f_molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, ra_ptr, err)

    if (error_has(err)) then
      call error_to_c(err, status, msg, msglen)
      return
    end if

    status = 0
    if (msglen > 0) msg(1) = c_null_char

    result_ra = c_loc(ra_ptr)
  end subroutine ra_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_rmsd_results_c
  !> @brief        Deallocate RMSD analysis results
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine deallocate_rmsd_results_c() bind(C, name="deallocate_rmsd_results_c")
    implicit none

    if (associated(ra_ptr)) then
      deallocate(ra_ptr)
      nullify(ra_ptr)
    end if

  end subroutine deallocate_rmsd_results_c

  subroutine ra_analysis_main( &
          molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, ra, err)
    use, intrinsic :: iso_c_binding
    implicit none
    type(s_molecule), intent(inout) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer,                intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer,                intent(in) :: ctrl_len
    real(wp), pointer, intent(out) :: ra(:)
    type(s_error),                   intent(inout) :: err

    ! local variables
    type(s_ctrl_data)      :: ctrl_data
    type(s_trajectory)     :: trajectory
    type(s_fitting)        :: fitting
    type(s_output)         :: output
    type(s_option)         :: option


    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.


    ! [Step1] Read control parameters from string
    !
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '

    call control_from_string(ctrl_text, ctrl_len, ctrl_data)


    ! [Step2] Set relevant variables and structures 
    !
    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '

    call setup(molecule, ctrl_data, fitting, output, option)


    ! [Step3] Analyze trajectory
    !
    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '

    call analyze(molecule, s_trajes_c, ana_period, output, option, &
                 fitting, ra, err)
    if (error_has(err)) return


    ! [Step4] Deallocate memory
    !
    write(MsgOut,'(A)') '[STEP4] Deallocate memory'
    write(MsgOut,'(A)') ' '

    call dealloc_trajectory(trajectory)
    call dealloc_molecules_all(molecule)
end subroutine ra_analysis_main

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in TRJ_ANALYSIS
  !! @authors      TM
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] output     : output information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(molecule, ctrl_data, fitting, output, option)
    use ra_control_mod
    use ra_option_mod
    use ra_option_str_mod
    use fitting_mod
    use fitting_str_mod
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
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_fitting),         intent(inout) :: fitting
    type(s_output),          intent(inout) :: output
    type(s_option),          intent(inout) :: option


    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)


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


    return

  end subroutine setup

end module rmsd_c_mod
