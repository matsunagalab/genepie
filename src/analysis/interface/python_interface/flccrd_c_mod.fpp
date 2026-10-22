!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   flccrd_c_mod
!> @brief   Python (bind(C)) entry point of flccrd_analysis
!! @authors Claude Code
!
!  The science lives in fa_analyze_mod (shared with the CLI program): the
!  entry point parses the control text ([SELECTION], [FITTING], [OPTION]),
!  runs analyze_flccrd_unified on an in-memory or lazily read trajectory and
!  fills the caller's variance-covariance matrix and RMSF buffers.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module flccrd_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use source_args_mod
  use fa_analyze_mod
  use fa_control_mod
  use fa_option_mod
  use fa_option_str_mod
  use trj_source_mod
  use fitting_mod
  use fitting_str_mod
  use select_mod
  use molecules_mod
  use molecules_str_mod
  use error_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: flccrd_analysis_c

  !> State shared between flccrd_analysis_c and its guarded body.
  type :: t_flccrd_ctx
    ! inputs
    type(c_ptr) :: molecule_ptr = c_null_ptr     ! caller's s_molecule_c
    type(t_source_args) :: src
    type(c_ptr) :: ctrl_text_ptr = c_null_ptr
    integer     :: ctrl_len = 0
    type(c_ptr) :: ave_coord_ptr = c_null_ptr    ! (3, natom)
    type(c_ptr) :: vcv_ptr = c_null_ptr          ! (nvcv, nvcv)
    integer     :: nvcv = 0
    type(c_ptr) :: rmsf_ptr = c_null_ptr         ! (nsel)
    integer     :: nsel = 0
    ! outputs
    integer     :: nstru_out = 0
    ! error state and resources released by the wrapper
    type(s_error)       :: err
    type(s_molecule)    :: f_molecule
    type(s_trj_source)  :: source
  end type t_flccrd_ctx

  private :: t_flccrd_ctx

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    flccrd_analysis_c
  !> @brief        Variance-covariance matrix of an in-memory or lazy DCD
  !!               trajectory about an averaged structure
  !! @authors      Claude Code
  !! @param[in]    molecule          : molecule C structure
  !! @param[in]    s_trajes_c        : trajectories C structure (memory mode)
  !! @param[in]    is_lazy           : 1 = read dcd_filename lazily
  !! @param[in]    dcd_filename, filename_len, trj_type, dcd_natom_expected,
  !!               source_selection_ptr, n_source_selection : lazy DCD source
  !! @param[in]    ana_period        : analysis period
  !! @param[in]    ctrl_text, ctrl_len : control text
  !! @param[in]    ave_coord_ptr     : averaged coordinates (3, natom)
  !! @param[in]    vcv_ptr, nvcv     : pre-allocated (nvcv, nvcv) matrix
  !! @param[in]    rmsf_ptr, nsel    : pre-allocated (nsel = nvcv/3) RMSF
  !! @param[out]   nstru_out         : number of structures analysed
  !! @param[out]   status, msg       : error status and message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine flccrd_analysis_c(molecule, s_trajes_c, is_lazy, &
                               dcd_filename, filename_len, trj_type, &
                               dcd_natom_expected, source_selection_ptr, &
                               n_source_selection, ana_period, &
                               ctrl_text, ctrl_len, &
                               ave_coord_ptr, vcv_ptr, nvcv, &
                               rmsf_ptr, nsel, nstru_out, &
                               status, msg, msglen) &
        bind(C, name="flccrd_analysis_c")
    implicit none

    ! Arguments
    type(s_molecule_c), intent(in), target :: molecule
    type(s_trajectories_c), intent(in), target :: s_trajes_c
    integer(c_int), value :: is_lazy
    character(kind=c_char), intent(in) :: dcd_filename(*)
    integer(c_int), value :: filename_len
    integer(c_int), value :: trj_type
    integer(c_int), value :: dcd_natom_expected
    type(c_ptr), value :: source_selection_ptr
    integer(c_int), value :: n_source_selection
    integer(c_int), value :: ana_period
    character(kind=c_char), intent(in), target :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), value :: ave_coord_ptr
    type(c_ptr), value :: vcv_ptr
    integer(c_int), value :: nvcv
    type(c_ptr), value :: rmsf_ptr
    integer(c_int), value :: nsel
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(t_flccrd_ctx), target :: c

    c%molecule_ptr  = c_loc(molecule)
    call pack_source_args(c%src, c_loc(s_trajes_c), is_lazy, dcd_filename, &
                          filename_len, trj_type, dcd_natom_expected, &
                          source_selection_ptr, n_source_selection, ana_period)
    c%ctrl_text_ptr = c_loc(ctrl_text)
    c%ctrl_len      = ctrl_len
    c%ave_coord_ptr = ave_coord_ptr
    c%vcv_ptr       = vcv_ptr
    c%nvcv          = nvcv
    c%rmsf_ptr      = rmsf_ptr
    c%nsel          = nsel

    call run_guarded(flccrd_analysis_body, c, c%err, status, msg, msglen)

    nstru_out = c%nstru_out

    call finalize_source(c%source)
    call dealloc_molecules_all(c%f_molecule)

  end subroutine flccrd_analysis_c

  !> Guarded body of flccrd_analysis_c (see run_guarded in error_mod).
  subroutine flccrd_analysis_body(ctx) bind(C, name="flccrd_analysis_body")
    implicit none
    type(c_ptr), value :: ctx

    type(t_flccrd_ctx), pointer :: c
    type(s_molecule_c), pointer :: molecule
    character(kind=c_char), pointer :: ctrl_text(:)
    type(s_ctrl_data)   :: ctrl_data
    type(s_fitting)     :: fitting
    type(s_option)      :: option
    real(wp), pointer   :: ave_coord(:,:), vcv(:,:), rmsf(:)
    integer             :: natom, nvcv_expected, nstru

    call c_f_pointer(ctx, c)
    call c_f_pointer(c%molecule_ptr, molecule)
    call c_f_pointer(c%ctrl_text_ptr, ctrl_text, [c%ctrl_len])

    if (.not. c_associated(c%ave_coord_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "flccrd_analysis_c: ave_coord_ptr is null")
      return
    end if
    if (.not. c_associated(c%vcv_ptr) .or. .not. c_associated(c%rmsf_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "flccrd_analysis_c: result pointers are null")
      return
    end if
    if (c%nvcv <= 0 .or. c%nsel * 3 /= c%nvcv) then
      call error_set(c%err, ERROR_DIMENSION, &
                     "flccrd_analysis_c: nvcv must equal 3 * nsel")
      return
    end if

    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    call c2f_s_molecule(molecule, c%f_molecule)
    natom = c%f_molecule%num_atoms

    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '
    call control_from_string(ctrl_text, c%ctrl_len, ctrl_data)

    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '
    call setup_selection(ctrl_data%sel_info, c%f_molecule)
    call setup_fitting(ctrl_data%fit_info, ctrl_data%sel_info, &
                       c%f_molecule, fitting)
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      c%f_molecule, option)

    nvcv_expected = flccrd_vcv_size(fitting, option)
    if (nvcv_expected /= c%nvcv) then
      call error_set(c%err, ERROR_DIMENSION, &
                     "flccrd_analysis_c: result size does not match the selection")
      call dealloc_option(option)
      return
    end if

    call c_f_pointer(c%ave_coord_ptr, ave_coord, [3, natom])
    call c_f_pointer(c%vcv_ptr, vcv, [c%nvcv, c%nvcv])
    call c_f_pointer(c%rmsf_ptr, rmsf, [c%nsel])

    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '
    call init_source_from_args(c%src, natom, 'flccrd_analysis_c', &
                               c%source, c%err)
    if (error_has(c%err)) then
      call dealloc_option(option)
      return
    end if

    call analyze_flccrd_unified(c%f_molecule, c%source, fitting, option, &
                                ave_coord, vcv, rmsf, nstru)
    c%nstru_out = nstru

    call dealloc_option(option)

  end subroutine flccrd_analysis_body

end module flccrd_c_mod
