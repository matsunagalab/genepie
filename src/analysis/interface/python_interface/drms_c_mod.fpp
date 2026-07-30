!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  dr_main
!! @brief   analysis of Drms
!! @authors Chigusa Kobayashi (CK), Daisuke Matsuoka (DM), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module drms_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use dr_analyze_mod           ! Use unified analysis from CLI module
  use trj_source_mod
  use result_sink_mod

  use trajectory_str_mod
  use molecules_str_mod
  use string_mod
  use error_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: drms_analysis_c
  public :: drms_analysis_lazy_c

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    drms_analysis_c
  !> @brief        DRMS analysis (zerocopy, pre-allocated result array)
  !! @authors      Claude Code
  !! @param[in]    contact_list_ptr : pointer to contact atom pairs (2, n_contact)
  !! @param[in]    contact_dist_ptr : pointer to reference distances
  !! @param[in]    n_contact        : number of contacts
  !! @param[in]    s_trajes_c       : trajectories C structure
  !! @param[in]    ana_period       : analysis period
  !! @param[in]    pbc_correct      : apply PBC correction (0 or 1)
  !! @param[in]    result_ptr       : pointer to pre-allocated result array
  !! @param[in]    result_size      : size of pre-allocated result array
  !! @param[out]   nstru_out        : actual number of structures analyzed
  !! @param[out]   status           : error status
  !! @param[out]   msg              : error message
  !! @param[in]    msglen           : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine drms_analysis_c(contact_list_ptr, contact_dist_ptr, &
                             n_contact, s_trajes_c, ana_period, &
                             pbc_correct, result_ptr, result_size, &
                             nstru_out, status, msg, msglen) &
        bind(C, name="drms_analysis_c")
    implicit none

    ! Arguments
    type(c_ptr), value :: contact_list_ptr
    type(c_ptr), value :: contact_dist_ptr
    integer(c_int), value :: n_contact
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer(c_int), value :: ana_period
    integer(c_int), value :: pbc_correct
    type(c_ptr), value :: result_ptr
    integer(c_int), value :: result_size
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_error) :: err
    type(s_trj_source) :: source
    type(s_result_sink) :: sink
    integer, pointer :: contact_list_f(:,:)
    real(wp), pointer :: contact_dist_f(:)
    real(wp), pointer :: result_f(:)
    logical :: pbc_flag
    integer :: nstru

    ! Initialize
    call error_init(err)
    status = 0
    nstru_out = 0

    ! Validate inputs
    if (n_contact <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_c: n_contact must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(contact_list_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_c: contact_list_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(contact_dist_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_c: contact_dist_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(result_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_c: result_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (result_size <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_c: result_size must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create zero-copy views of arrays from Python
    call C_F_POINTER(contact_list_ptr, contact_list_f, [2, n_contact])
    call C_F_POINTER(contact_dist_ptr, contact_dist_f, [n_contact])
    call C_F_POINTER(result_ptr, result_f, [result_size])

    ! Convert pbc_correct to logical
    pbc_flag = (pbc_correct /= 0)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Initialize source (memory mode) and sink (array mode)
    call init_source_memory(source, s_trajes_c%coords, s_trajes_c%pbc_boxes, &
                            s_trajes_c%natom, s_trajes_c%nframe, ana_period)
    call init_sink_array(sink, result_f, result_size)

    ! Run unified DRMS analysis
    write(MsgOut,'(A)') '[STEP1] DRMS Analysis (unified)'
    write(MsgOut,'(A)') ' '

    call analyze_drms_unified(source, sink, contact_list_f, contact_dist_f, &
                              n_contact, pbc_flag, nstru)

    nstru_out = nstru

    ! Cleanup
    call finalize_sink(sink)
    call finalize_source(source)

  end subroutine drms_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    drms_analysis_lazy_c
  !> @brief        DRMS analysis with lazy DCD loading (memory efficient)
  !! @authors      Claude Code
  !! @param[in]    dcd_filename     : DCD file path (C string)
  !! @param[in]    filename_len     : length of filename
  !! @param[in]    trj_type         : trajectory type (1=COOR, 2=COOR+BOX)
  !! @param[in]    contact_list_ptr : pointer to contact atom pairs (2, n_contact)
  !! @param[in]    contact_dist_ptr : pointer to reference distances
  !! @param[in]    n_contact        : number of contacts
  !! @param[in]    n_atoms          : number of atoms (for DCD validation)
  !! @param[in]    ana_period       : analysis period
  !! @param[in]    pbc_correct      : apply PBC correction (0 or 1)
  !! @param[in]    result_ptr       : pointer to pre-allocated result array
  !! @param[in]    result_size      : size of result array
  !! @param[out]   nstru_out        : actual number of structures analyzed
  !! @param[out]   dcd_nframe_out   : total frames in DCD
  !! @param[out]   dcd_natom_out    : atoms per frame in DCD
  !! @param[out]   status           : error status
  !! @param[out]   msg              : error message
  !! @param[in]    msglen           : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine drms_analysis_lazy_c(dcd_filename, filename_len, trj_type, &
                                  dcd_natom_expected, source_selection_ptr, &
                                  n_source_selection, &
                                  contact_list_ptr, contact_dist_ptr, &
                                  n_contact, n_atoms, ana_period, &
                                  pbc_correct, result_ptr, result_size, &
                                  nstru_out, dcd_nframe_out, dcd_natom_out, &
                                  status, msg, msglen) &
        bind(C, name="drms_analysis_lazy_c")
    implicit none

    ! Arguments
    character(kind=c_char), intent(in) :: dcd_filename(*)
    integer(c_int), value :: filename_len
    integer(c_int), value :: trj_type
    integer(c_int), value :: dcd_natom_expected
    type(c_ptr), value :: source_selection_ptr
    integer(c_int), value :: n_source_selection
    type(c_ptr), value :: contact_list_ptr
    type(c_ptr), value :: contact_dist_ptr
    integer(c_int), value :: n_contact
    integer(c_int), value :: n_atoms
    integer(c_int), value :: ana_period
    integer(c_int), value :: pbc_correct
    type(c_ptr), value :: result_ptr
    integer(c_int), value :: result_size
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: dcd_nframe_out
    integer(c_int), intent(out) :: dcd_natom_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_error) :: err
    type(s_trj_source) :: source
    type(s_result_sink) :: sink
    character(MaxFilename) :: filename_f
    integer, pointer :: contact_list_f(:,:)
    integer, pointer :: source_selection_f(:)
    real(wp), pointer :: contact_dist_f(:)
    real(wp), pointer :: result_f(:)
    logical :: pbc_flag
    integer :: nstru, init_status
    integer(c_int) :: grc

    ! Guard the whole body: init_source_lazy_dcd opens the DCD file, and a
    ! missing/unreadable file calls error_msg -> exit(1) in CLI mode, which
    ! would kill the host Python process.
    grc = fi_error_guard_run(c_funloc(run_body))
    if (grc /= 0) then
      call error_from_pending(err)
      call error_to_c(err, status, msg, msglen)
      call finalize_sink(sink)
      call finalize_source(source)
    end if
    return

  contains
    subroutine run_body() bind(C)

    ! Initialize
    call error_init(err)
    status = 0
    nstru_out = 0
    dcd_nframe_out = 0
    dcd_natom_out = 0

    ! Convert C string to Fortran string
    call c_filename_to_fortran(dcd_filename, filename_len, filename_f)

    ! Validate inputs
    if (n_contact <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_lazy_c: n_contact must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (n_source_selection /= n_atoms .or. &
        .not. c_associated(source_selection_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_lazy_c: invalid source selection")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(contact_list_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_lazy_c: contact_list_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(contact_dist_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_lazy_c: contact_dist_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(result_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_lazy_c: result_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (result_size <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_lazy_c: result_size must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create zero-copy views of arrays from Python
    call C_F_POINTER(contact_list_ptr, contact_list_f, [2, n_contact])
    call C_F_POINTER(contact_dist_ptr, contact_dist_f, [n_contact])
    call C_F_POINTER(result_ptr, result_f, [result_size])
    call C_F_POINTER(source_selection_ptr, source_selection_f, &
                     [n_source_selection])
    if (any(source_selection_f < 1) .or. &
        any(source_selection_f > dcd_natom_expected)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_lazy_c: source selection out of range")
      call error_to_c(err, status, msg, msglen)
      return
    end if
    if (any(contact_list_f < 1) .or. any(contact_list_f > n_atoms)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_lazy_c: contact index out of range")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Convert pbc_correct to logical
    pbc_flag = (pbc_correct /= 0)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Initialize lazy DCD source
    write(MsgOut,'(A)') '[STEP1] Initialize Lazy DCD Source for DRMS'
    write(MsgOut,'(A)') ' '

    call init_source_lazy_dcd(source, trim(filename_f), trj_type, ana_period, &
                              source_selection_f, n_source_selection, init_status)
    if (init_status /= 0) then
      call error_set(err, init_status, &
                     "drms_analysis_lazy_c: unable to initialize DCD source")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Return DCD info
    dcd_nframe_out = source%dcd_nframe
    dcd_natom_out = source%dcd_natom

    ! Check atom count
    if (source%dcd_natom /= dcd_natom_expected) then
      call error_set(err, ERROR_ATOM_COUNT, &
                     "drms_analysis_lazy_c: atom count mismatch")
      call error_to_c(err, status, msg, msglen)
      call finalize_source(source)
      return
    end if

    ! Initialize sink (array mode)
    call init_sink_array(sink, result_f, result_size)

    ! Run unified DRMS analysis (lazy loading via source abstraction)
    write(MsgOut,'(A)') '[STEP2] DRMS Analysis (lazy loading, unified)'
    write(MsgOut,'(A)') ' '

    call analyze_drms_unified(source, sink, contact_list_f, contact_dist_f, &
                              n_contact, pbc_flag, nstru)

    nstru_out = nstru

    ! Cleanup
    call finalize_sink(sink)
    call finalize_source(source)

    end subroutine run_body
  end subroutine drms_analysis_lazy_c

end module drms_c_mod
