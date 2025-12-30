!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  rg_main
!! @brief   RG analysis
!! @authors Motoshi Kamiya (MK), Takaharu Mori (TM), Yuji Sugita (YS)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rg_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use rg_impl_mod

  use rg_control_mod
  use rg_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: rg_analysis_c
  public :: rg_analysis_zerocopy_c
  public :: deallocate_rg_results_c

  ! Module-level pointer for results (to be deallocated later)
  real(wp), pointer, save :: rg_ptr(:) => null()

contains
  subroutine rg_analysis_c(molecule, s_trajes_c, ana_period, &
                           ctrl_text, ctrl_len, result_rg) &
        bind(C, name="rg_analysis_c")
    use conv_f_c_util
    implicit none
    type(s_molecule_c), intent(in) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer, intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out) :: result_rg

    type(s_molecule) :: f_molecule

    ! Deallocate previous results if any
    if (associated(rg_ptr)) then
      deallocate(rg_ptr)
      nullify(rg_ptr)
    end if

    call c2f_s_molecule(molecule, f_molecule)
    call rg_analysis_main( &
        f_molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, rg_ptr)
    result_rg = c_loc(rg_ptr)
  end subroutine rg_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rg_analysis_zerocopy_c
  !> @brief        RG analysis with true zero-copy (mass array passed directly)
  !! @authors      Claude Code
  !! @param[in]    mass_ptr        : pointer to mass array (from Python NumPy)
  !! @param[in]    n_atoms         : number of atoms
  !! @param[in]    s_trajes_c      : trajectories C structure
  !! @param[in]    ana_period      : analysis period
  !! @param[in]    analysis_idx    : analysis atom indices (1-indexed)
  !! @param[in]    n_analysis      : number of analysis atoms
  !! @param[in]    mass_weighted   : use mass weighting (0 or 1)
  !! @param[out]   result_rg       : pointer to result array
  !! @param[out]   status          : error status
  !! @param[out]   msg             : error message
  !! @param[in]    msglen          : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rg_analysis_zerocopy_c(mass_ptr, n_atoms, s_trajes_c, ana_period, &
                                    analysis_idx, n_analysis, mass_weighted, &
                                    result_rg, status, msg, msglen) &
        bind(C, name="rg_analysis_zerocopy_c")
    use error_mod, only: s_error, error_init, error_set, &
                         ERROR_INVALID_PARAM, error_to_c
    implicit none

    ! Arguments
    type(c_ptr), value :: mass_ptr
    integer(c_int), value :: n_atoms
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer(c_int), value :: ana_period
    type(c_ptr), value :: analysis_idx
    integer(c_int), value :: n_analysis
    integer(c_int), value :: mass_weighted
    type(c_ptr), intent(out) :: result_rg
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_error) :: err
    real(wp), pointer :: mass_f(:)
    integer, pointer :: idx_f(:)
    integer, allocatable :: idx_copy(:)
    logical :: use_mass

    ! Initialize
    call error_init(err)
    status = 0
    result_rg = c_null_ptr

    ! Deallocate previous results if any
    if (associated(rg_ptr)) then
      deallocate(rg_ptr)
      nullify(rg_ptr)
    end if

    ! Validate inputs
    if (n_analysis <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_zerocopy_c: n_analysis must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(mass_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_zerocopy_c: mass_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create zero-copy view of mass array from Python
    call C_F_POINTER(mass_ptr, mass_f, [n_atoms])

    ! Convert analysis indices from C pointer to Fortran array
    call C_F_POINTER(analysis_idx, idx_f, [n_analysis])
    allocate(idx_copy(n_analysis))
    idx_copy(:) = idx_f(:)

    ! Convert mass_weighted to logical
    use_mass = (mass_weighted /= 0)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Run analysis
    write(MsgOut,'(A)') '[STEP1] RG Analysis (zero-copy interface)'
    write(MsgOut,'(A)') ' '

    call analyze_zerocopy(mass_f, s_trajes_c, ana_period, &
                          idx_copy, n_analysis, use_mass, rg_ptr)

    result_rg = c_loc(rg_ptr)

    ! Cleanup local arrays (mass_f is a view, don't deallocate)
    deallocate(idx_copy)

  end subroutine rg_analysis_zerocopy_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_rg_results_c
  !> @brief        Deallocate RG analysis results
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine deallocate_rg_results_c() bind(C, name="deallocate_rg_results_c")
    implicit none

    if (associated(rg_ptr)) then
      deallocate(rg_ptr)
      nullify(rg_ptr)
    end if

  end subroutine deallocate_rg_results_c

  subroutine rg_analysis_main( &
          molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, rg)
    use, intrinsic :: iso_c_binding
    implicit none
    type(s_molecule), intent(inout) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer,                intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer,                intent(in) :: ctrl_len
    real(wp), pointer, intent(out) :: rg(:)

    ! local variables
    type(s_ctrl_data)      :: ctrl_data
    type(s_trajectory)     :: trajectory
    type(s_output)         :: output
    type(s_option)         :: option


    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! [STEP1] Read control parameters from string
    !
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '

    call control_from_string(ctrl_text, ctrl_len, ctrl_data)

    ! [STEP2] Set variables and structures
    !
    write(MsgOut,'(A)') '[STEP2] Set Variables and Structures'
    write(MsgOut,'(A)') ' '

    call setup(molecule, ctrl_data, output, option)

    ! [STEP3] Analyze trajectory
    !
    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '

    call analyze(molecule, s_trajes_c, ana_period, output, option, &
                 rg)


    ! [STEP4] Deallocate memory
    !
    write(MsgOut,'(A)') '[STEP4] Deallocate memory'
    write(MsgOut,'(A)') ' '

    call dealloc_trajectory(trajectory)
    call dealloc_molecules_all(molecule)
  end subroutine rg_analysis_main

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

  subroutine setup(molecule, ctrl_data, output, option)
    use rg_control_mod
    use rg_option_mod
    use rg_option_str_mod
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
                      molecule, option)

    return

  end subroutine setup

end module rg_c_mod
