!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  kc_main
!! @brief   analysis trajectory files
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module kmeans_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use kmeans_impl_mod

  use kc_control_mod
  use kc_option_str_mod
  use fitting_str_mod
  use trajectory_str_mod
  use input_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use error_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

contains
  subroutine kc_analysis_c(molecule, s_trajes_c, ana_period, &
          ctrl_text, ctrl_len, &
          out_pdb_ptr, cluster_index, size_cluster_index, status, msg, msglen)  &
        bind(C, name="kc_analysis_c")
    use conv_f_c_util
    implicit none
    type(s_molecule_c), intent(in) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer(c_int), intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out) :: out_pdb_ptr
    type(c_ptr), intent(out) :: cluster_index
    integer(c_int), intent(out) :: size_cluster_index
    integer(c_int),          intent(out) :: status
    character(kind=c_char),  intent(out) :: msg(*)
    integer(c_int),          value       :: msglen

    type(s_molecule) :: f_molecule
    character(len=:), allocatable :: out_pdb_f
    character(kind=c_char), pointer :: out_pdb_c(:)
    integer, allocatable :: cluster_index_f(:)

    type(s_error) :: err

    call error_init(err)

    call c2f_s_molecule(molecule, f_molecule)
    call kc_analysis_main( &
        f_molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, &
        out_pdb_f, cluster_index_f, err)

    if (error_has(err)) then
      call error_to_c(err, status, msg, msglen)
      return
    end if

    status = 0
    if (msglen > 0) msg(1) = c_null_char

    if (allocated(out_pdb_f)) then
      call f2c_string(out_pdb_f, out_pdb_c)
      out_pdb_ptr = c_loc(out_pdb_c(1))
    else
      out_pdb_ptr = c_null_ptr
    end if
    if (allocated(cluster_index_f)) then
        cluster_index = f2c_int_array(cluster_index_f)
        size_cluster_index = size(cluster_index_f)
    else
        cluster_index = c_null_ptr
        size_cluster_index = 0
    end if
  end subroutine kc_analysis_c

  subroutine kc_analysis_main( &
          molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, &
          out_pdb, cluster_index, err)
    implicit none
    type(s_molecule), intent(inout) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer,                intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer,                intent(in) :: ctrl_len
    character(len=:), allocatable, intent(out) :: out_pdb
    integer, allocatable,     intent(out) :: cluster_index(:)
    type(s_error),                   intent(inout) :: err

    ! local variables
    type(s_ctrl_data)      :: ctrl_data
    type(s_input)          :: input
    type(s_trajectory)     :: trajectory
    type(s_fitting)        :: fitting
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

    call setup(ctrl_data,  &
               molecule,   &
               input,      &
               fitting,    &
               option,     &
               output)


    ! [Step3] Analyze trajectory
    !
    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '

    call analyze(molecule,   &
                 input,      &
                 s_trajes_c, &
                 ana_period, &
                 fitting,    &
                 option,     &
                 output,     &
                 out_pdb,    &
                 cluster_index, &
                 err)
    if (error_has(err)) return


    ! [Step4] Deallocate memory
    !
    write(MsgOut,'(A)') '[STEP4] Deallocate memory'
    write(MsgOut,'(A)') ' '

    call dealloc_trajectory(trajectory)
    call dealloc_molecules_all(molecule)
    call dealloc_option(option)
end subroutine kc_analysis_main

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures
  !! @authors      TM
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] input      : input information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data,  &
                   molecule,   &
                   input,      &
                   fitting,    &
                   option,     &
                   output)
    use kc_control_mod
    use kc_option_mod
    use kc_option_str_mod
    use fitting_mod
    use fitting_str_mod
    use input_mod
    use input_str_mod
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
    use fileio_grotop_mod
    use fileio_par_mod
    use fileio_psf_mod
    use fileio_pdb_mod
    implicit none

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_input),           intent(inout) :: input
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output


    ! setup input
    !
    call setup_input(ctrl_data%inp_info, input)

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


    return

  end subroutine setup

end module kmeans_c_mod

