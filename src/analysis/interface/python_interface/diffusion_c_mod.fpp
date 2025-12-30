!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  ma_main
!! @brief   MSD analysis
!! @authors Donatas Surblys (DS), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module diffusion_c_mod

  use, intrinsic :: iso_c_binding
  use constants_mod
  use diffusion_impl_mod
  use error_mod
  use messages_mod

  implicit none

  public :: diffusion_analysis_c

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    diffusion_analysis_c
  !> @brief        Diffusion analysis (zerocopy, pre-allocated results)
  !! @authors      Donatas Surblys, Claude Code
  !! @param[in]    msd_ptr        : pointer to MSD data (ncols, ndata)
  !! @param[in]    ndata          : number of data points
  !! @param[in]    ncols          : number of columns
  !! @param[in]    time_step      : time step in ps
  !! @param[in]    distance_unit  : distance unit conversion factor
  !! @param[in]    ndofs          : number of degrees of freedom
  !! @param[in]    start_step     : start step for fitting
  !! @param[in]    stop_step      : stop step for fitting
  !! @param[in]    out_data_ptr   : pre-allocated output data (out_ncols, ndata)
  !! @param[in]    out_data_size  : total size of out_data array
  !! @param[in]    diff_coeff_ptr : pre-allocated diffusion coefficients (n_sets)
  !! @param[in]    n_sets         : number of MSD sets (ncols - 1)
  !! @param[out]   status         : error status
  !! @param[out]   msg            : error message
  !! @param[in]    msglen         : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine diffusion_analysis_c(msd_ptr, ndata, ncols, &
                                  time_step, distance_unit, ndofs, &
                                  start_step, stop_step, &
                                  out_data_ptr, out_data_size, &
                                  diff_coeff_ptr, n_sets, &
                                  status, msg, msglen) &
        bind(C, name="diffusion_analysis_c")
    implicit none

    ! Arguments
    type(c_ptr), value :: msd_ptr
    integer(c_int), value :: ndata
    integer(c_int), value :: ncols
    real(c_double), value :: time_step
    real(c_double), value :: distance_unit
    integer(c_int), value :: ndofs
    integer(c_int), value :: start_step
    integer(c_int), value :: stop_step
    type(c_ptr), value :: out_data_ptr
    integer(c_int), value :: out_data_size
    type(c_ptr), value :: diff_coeff_ptr
    integer(c_int), value :: n_sets
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_error) :: err
    real(wp), pointer :: msd_f(:,:)
    real(wp), pointer :: out_data_f(:,:)
    real(wp), pointer :: diff_coeff_f(:)
    integer :: out_ncols

    ! Initialize
    call error_init(err)
    status = 0

    ! Validate inputs
    if (.not. c_associated(msd_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_c: msd_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (ndata <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_c: ndata must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (ncols < 2) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_c: ncols must be >= 2")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(out_data_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_c: out_data_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(diff_coeff_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_c: diff_coeff_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create views of arrays from Python
    call C_F_POINTER(msd_ptr, msd_f, [ncols, ndata])

    ! Calculate output dimensions
    out_ncols = 2 * n_sets + 1  ! time + (msd + fit) * n_sets
    call C_F_POINTER(out_data_ptr, out_data_f, [out_ncols, ndata])
    call C_F_POINTER(diff_coeff_ptr, diff_coeff_f, [n_sets])

    ! Run analysis
    write(MsgOut,'(A)') '[STEP1] Diffusion Analysis'
    write(MsgOut,'(A)') ' '

    call analyze(msd_f, ndata, ncols, &
                 time_step, distance_unit, ndofs, &
                 start_step, stop_step, &
                 out_data_f, diff_coeff_f)

  end subroutine diffusion_analysis_c

end module diffusion_c_mod
