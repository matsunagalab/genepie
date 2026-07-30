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

module mbar_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use mbar_impl_mod

  use mbar_control_mod
  use mbar_option_str_mod
  use output_str_mod
  use input_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use error_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

 contains
  subroutine mbar_analysis_c(ctrl_text, ctrl_len, return_weights, result_fene, &
                             n_replica, n_blocks, result_weights, &
                             n_weight_replica, n_weight_step, &
                             status, msg, msglen) &
        bind(C, name="mbar_analysis_c")
    use conv_f_c_util
    implicit none
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    integer(c_int), value :: return_weights
    type(c_ptr), intent(out) :: result_fene
    integer(c_int), intent(out) :: n_replica
    integer(c_int), intent(out) :: n_blocks
    type(c_ptr), intent(out) :: result_weights
    integer(c_int), intent(out) :: n_weight_replica
    integer(c_int), intent(out) :: n_weight_step
    integer(c_int),          intent(out) :: status
    character(kind=c_char),  intent(out) :: msg(*)
    integer(c_int),          value       :: msglen


    real(wp), pointer :: fene_f(:,:) => null()
    real(wp), pointer :: weights_f(:,:) => null()

    type(s_error) :: err
    integer(c_int) :: grc

    result_fene = c_null_ptr
    result_weights = c_null_ptr
    n_replica = 0
    n_blocks = 0
    n_weight_replica = 0
    n_weight_step = 0

    call error_init(err)

    ! Run under the library-mode error guard (see error_mod / fileio_data_.c)
    ! so a fatal error_msg becomes a catchable error rather than exit(1).
    grc = fi_error_guard_run(c_funloc(run_body))
    if (grc /= 0) then
      call error_from_pending(err)
      call error_finish_to_c(err, status, msg, msglen)
      return
    end if

    call error_finish_to_c(err, status, msg, msglen)
    if (error_has(err)) return

    if (associated(fene_f)) result_fene = c_loc(fene_f)
    if (associated(weights_f)) result_weights = c_loc(weights_f)
  contains
    subroutine run_body() bind(C)
      call mbar_analysis_main( &
          ctrl_text, ctrl_len, return_weights /= 0, &
          fene_f, n_replica, n_blocks, weights_f, &
          n_weight_replica, n_weight_step, err)
    end subroutine run_body
  end subroutine mbar_analysis_c

  subroutine mbar_analysis_main( &
          ctrl_text, ctrl_len, return_weights, result_fene, &
          n_replica, n_blocks, result_weights, &
          n_weight_replica, n_weight_step, err)
    implicit none
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer,                intent(in) :: ctrl_len
    logical,                intent(in) :: return_weights
    real(wp), pointer, intent(out) :: result_fene(:,:)
    integer,           intent(out) :: n_replica
    integer,           intent(out) :: n_blocks
    real(wp), pointer, intent(out) :: result_weights(:,:)
    integer,           intent(out) :: n_weight_replica
    integer,           intent(out) :: n_weight_step
    type(s_error),     intent(inout) :: err


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

    call control_from_string(ctrl_text, ctrl_len, ctrl_data)


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
    call analyze(molecule, input, output, option, return_weights, &
                 result_fene, n_replica, n_blocks, result_weights, &
                 n_weight_replica, n_weight_step, err)
    if (error_has(err)) return


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

end module mbar_c_mod
