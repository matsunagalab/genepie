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
  use drms_impl_mod

  use trajectory_str_mod
  use molecules_str_mod
  use error_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: drms_analysis_c

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
    integer, pointer :: contact_list_f(:,:)
    real(wp), pointer :: contact_dist_f(:)
    real(wp), pointer :: result_f(:)
    logical :: pbc_flag
    integer :: nstru_local

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

    ! Run analysis
    write(MsgOut,'(A)') '[STEP1] DRMS Analysis'
    write(MsgOut,'(A)') ' '

    call analyze(contact_list_f, contact_dist_f, n_contact, &
                 s_trajes_c, ana_period, pbc_flag, &
                 result_f, nstru_local)

    nstru_out = nstru_local

  end subroutine drms_analysis_c

end module drms_c_mod
