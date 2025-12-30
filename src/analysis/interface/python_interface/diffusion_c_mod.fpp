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
  use da_setup_mod
  use da_control_mod
  use da_option_str_mod
  use input_str_mod
  use output_str_mod
  use error_mod
  use string_mod
  use messages_mod

  implicit none

  public :: diffusion_analysis_c
  public :: diffusion_analysis_zerocopy_c
  public :: deallocate_diffusion_results_c

  ! Module-level pointers for zerocopy results (to be deallocated later)
  real(wp), pointer, save :: diffusion_out_ptr(:,:) => null()
  real(wp), pointer, save :: diffusion_coeff_ptr(:) => null()

contains

  subroutine diffusion_analysis_c(c_msd, msd_dim1, msd_dim2, &
          ctrl_text, ctrl_len, &
          out_data, status, msg, msglen) &
        bind(C, name="diffusion_analysis_c")
    use conv_f_c_util
    implicit none
    type(c_ptr), intent(in) :: c_msd
    integer(c_int), intent(in) :: msd_dim1
    integer(c_int), intent(in) :: msd_dim2
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out) :: out_data
    integer(c_int),          intent(out) :: status
    character(kind=c_char),  intent(out) :: msg(*)
    integer(c_int),          value       :: msglen

    real(wp), pointer :: f_msd(:,:)
    real(wp), pointer :: f_out_data(:,:)

    type(s_error) :: err

    call error_init(err)
    call C_F_POINTER(c_msd, f_msd, [msd_dim2, msd_dim1])
    allocate(f_out_data( &
        (size(f_msd, dim=1) - 1) * 2 + 1, &
        size(f_msd, dim=2)))
    call diffusion_analysis_main(f_msd, ctrl_text, ctrl_len, f_out_data, err)
    if (error_has(err)) then
      call error_to_c(err, status, msg, msglen)
      return
    end if
    status = 0
    if (msglen > 0) msg(1) = c_null_char

    out_data = f2c_double_array(f_out_data)
  end subroutine diffusion_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    diffusion_analysis_zerocopy_c
  !> @brief        Diffusion analysis with direct parameters (zerocopy version)
  !! @authors      Claude Code
  !! @param[in]    msd_ptr       : pointer to MSD data array (from Python NumPy)
  !! @param[in]    ndata         : number of data points
  !! @param[in]    ncols         : number of columns (time + MSD sets)
  !! @param[in]    time_step     : time per step in ps
  !! @param[in]    distance_unit : distance unit factor (1.0 for Angstrom)
  !! @param[in]    ndofs         : degrees of freedom (3 for 3D)
  !! @param[in]    start_step    : start step for fitting (1-indexed)
  !! @param[in]    stop_step     : stop step for fitting (1-indexed)
  !! @param[out]   out_data_ptr  : pointer to output data array
  !! @param[out]   diff_coeff_ptr: pointer to diffusion coefficients array
  !! @param[out]   status        : error status
  !! @param[out]   msg           : error message
  !! @param[in]    msglen        : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine diffusion_analysis_zerocopy_c(msd_ptr, ndata, ncols, &
                                           time_step, distance_unit, ndofs, &
                                           start_step, stop_step, &
                                           out_data_ptr, diff_coeff_ptr, &
                                           status, msg, msglen) &
        bind(C, name="diffusion_analysis_zerocopy_c")
    use conv_f_c_util
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
    type(c_ptr), intent(out) :: out_data_ptr
    type(c_ptr), intent(out) :: diff_coeff_ptr
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_error) :: err
    real(wp), pointer :: msd_f(:,:)
    integer :: n_sets, out_ncols

    ! Initialize
    call error_init(err)
    status = 0
    out_data_ptr = c_null_ptr
    diff_coeff_ptr = c_null_ptr

    ! Deallocate previous results if any
    if (associated(diffusion_out_ptr)) then
      deallocate(diffusion_out_ptr)
      nullify(diffusion_out_ptr)
    end if
    if (associated(diffusion_coeff_ptr)) then
      deallocate(diffusion_coeff_ptr)
      nullify(diffusion_coeff_ptr)
    end if

    ! Validate inputs
    if (.not. c_associated(msd_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_zerocopy_c: msd_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (ndata <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_zerocopy_c: ndata must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (ncols < 2) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_zerocopy_c: ncols must be >= 2")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create view of MSD data from Python (Fortran order: ncols x ndata)
    call C_F_POINTER(msd_ptr, msd_f, [ncols, ndata])

    ! Calculate output dimensions
    n_sets = ncols - 1
    out_ncols = 2 * n_sets + 1  ! time + (msd + fit) * n_sets

    ! Allocate output arrays (using module-level pointers)
    allocate(diffusion_out_ptr(out_ncols, ndata))
    allocate(diffusion_coeff_ptr(n_sets))

    ! Run analysis
    write(MsgOut,'(A)') '[STEP1] Diffusion Analysis (zerocopy interface)'
    write(MsgOut,'(A)') ' '

    call analyze_zerocopy(msd_f, ndata, ncols, &
                          real(time_step, wp), real(distance_unit, wp), ndofs, &
                          start_step, stop_step, &
                          diffusion_out_ptr, diffusion_coeff_ptr)

    ! Return pointers to allocated arrays
    out_data_ptr = c_loc(diffusion_out_ptr)
    diff_coeff_ptr = c_loc(diffusion_coeff_ptr)

    ! Note: msd_f is a view, don't deallocate
    ! Module-level pointers will be deallocated by Python calling
    ! deallocate_diffusion_results_c

  end subroutine diffusion_analysis_zerocopy_c

  subroutine diffusion_analysis_main(msd_data, ctrl_text, ctrl_len, out_data, err)
    use, intrinsic :: iso_c_binding
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer, intent(in) :: ctrl_len
    real(wp), intent(in) :: msd_data(:,:)
    real(wp), intent(out) :: out_data(:,:)
    type(s_error), intent(inout) :: err
    ! local variables
    type(s_ctrl_data)      :: ctrl_data
    type(s_input)          :: input
    type(s_output)         :: output
    type(s_option)         :: option


    ! [Step1] Read control file
    !
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '

    call control_from_string(ctrl_text, ctrl_len, ctrl_data)


    ! [Step2] Set relevant variables and structures
    !
    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '

    call setup(ctrl_data, input, output, option)


    ! [Step3] Analysis of mean square dispacement
    !
    write(MsgOut,'(A)') '[STEP3] Analysis of mean square dispacement'
    call analyze(msd_data, input, option, out_data, err)


    ! [Step4] Deallocate memory
    !
    write(MsgOut,'(A)') '[STEP5] Deallocate memory'
    write(MsgOut,'(A)') ' '

    call dealloc_option(option)

  end subroutine diffusion_analysis_main

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_diffusion_results_c
  !> @brief        Deallocate diffusion analysis results
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine deallocate_diffusion_results_c() bind(C, name="deallocate_diffusion_results_c")
    implicit none

    if (associated(diffusion_out_ptr)) then
      deallocate(diffusion_out_ptr)
      nullify(diffusion_out_ptr)
    end if

    if (associated(diffusion_coeff_ptr)) then
      deallocate(diffusion_coeff_ptr)
      nullify(diffusion_coeff_ptr)
    end if

  end subroutine deallocate_diffusion_results_c

end module diffusion_c_mod
