!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   da_analyze_mod
!> @brief   analyze mean square displacement
!! @authors Donatas Surblys (DS), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif


module diffusion_impl_mod

  use da_option_str_mod
  use input_str_mod
  use output_str_mod
  use fileio_mod
  use error_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  private :: fit_least_squares

  ! fitting result for a + bx
  type :: s_fitting_result
    real(wp), dimension(2) :: coeff   ! coefficients of the fitting (a, b)
    real(wp), dimension(2) :: stderr  ! standard error of each coefficient
    real(wp)               :: corr    ! correlation coefficient, showing
                                      ! the quality of the fit
  end type s_fitting_result

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        Diffusion analysis (zerocopy, pre-allocated results)
  !! @authors      Donatas Surblys (DS), Takaharu Mori (TM), Claude Code
  !! @param[in]    msd_data       : MSD data array (ncols, ndata)
  !! @param[in]    ndata          : number of data points
  !! @param[in]    ncols          : number of columns
  !! @param[in]    time_step      : time step in ps
  !! @param[in]    distance_unit  : distance unit conversion factor
  !! @param[in]    ndofs          : number of degrees of freedom
  !! @param[in]    start_step     : start step for fitting
  !! @param[in]    stop_step      : stop step for fitting
  !! @param[inout] out_data       : pre-allocated output data (out_ncols, ndata)
  !! @param[inout] diffusion_coeff: pre-allocated diffusion coefficients (n_sets)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(msd_data, ndata, ncols, &
                     time_step, distance_unit, ndofs, &
                     start_step, stop_step, &
                     out_data, diffusion_coeff)

    ! formal arguments
    real(wp), intent(in)    :: msd_data(:,:)
    integer,  intent(in)    :: ndata
    integer,  intent(in)    :: ncols
    real(wp), intent(in)    :: time_step
    real(wp), intent(in)    :: distance_unit
    integer,  intent(in)    :: ndofs
    integer,  intent(in)    :: start_step
    integer,  intent(in)    :: stop_step
    real(wp), intent(inout) :: out_data(:,:)      ! pre-allocated
    real(wp), intent(inout) :: diffusion_coeff(:) ! pre-allocated

    ! local variables
    integer                                :: n_sets
    integer                                :: iset, idata
    integer                                :: idx_start_fit, idx_stop_fit
    real(wp)                               :: d_coeff
    real(wp), allocatable, dimension(:, :) :: xydata
    character(*), parameter                :: format_float = "es25.16e3"
    type(s_fitting_result), dimension(:), allocatable :: fittings


    ! Number of MSD sets (columns minus time column)
    n_sets = ncols - 1

    allocate(fittings(n_sets))
    allocate(xydata(ncols, ndata))

    ! Copy and convert data
    xydata = msd_data

    ! Convert to ps and angstroms
    xydata(1, :)  = xydata(1, :)  * time_step
    xydata(2:, :) = xydata(2:, :) * distance_unit ** 2
    do iset = 2, ncols
      xydata(iset, :) = xydata(iset, :) / (2.0_wp * ndofs)
    end do

    ! Determine fitting range
    idx_start_fit = max(start_step, 1)
    idx_stop_fit = min(stop_step, ndata)

    ! Analyze
    write(MsgOut, '()')
    write(MsgOut, '("Analyze> Starting fit at",es9.2e2," ps")') &
      xydata(1, idx_start_fit)
    write(MsgOut, '("Analyze> Ending fit at  ",es9.2e2," ps")') &
      xydata(1, idx_stop_fit)
    write(MsgOut, '()')

    ! Perform least squares fitting for each MSD set
    do iset = 1, n_sets
      fittings(iset) = fit_least_squares(xydata(1, idx_start_fit:idx_stop_fit), &
                                         xydata(iset+1, idx_start_fit:idx_stop_fit), .true.)

      ! Convert diffusion coefficient to cm^2/s
      d_coeff = fittings(iset)%coeff(2) * 1e-4_wp
      diffusion_coeff(iset) = d_coeff

      write(MsgOut, '("Analyze> Diffusion coefficient",i4," =",&
                      &'//format_float//'," cm^2/s, (corr.",f9.6,")")') &
        iset, d_coeff, fittings(iset)%corr
    end do

    write(MsgOut, '()')

    ! Copy data to output
    ! Column 1: time
    out_data(1, :) = xydata(1, :)
    ! Columns 2, 4, 6, ...: MSD data
    ! Columns 3, 5, 7, ...: fit data
    do iset = 1, n_sets
      out_data(2*iset, :) = xydata(iset+1, :)
      do idata = 1, ndata
        out_data(2*iset+1, idata) = &
            fittings(iset)%coeff(1)+fittings(iset)%coeff(2)*xydata(1, idata)
      end do
    end do

    ! Cleanup
    deallocate(xydata)
    deallocate(fittings)

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      fit_least_squares
  !> @brief        fit line to data using least squares method
  !! @authors      DS
  !! @param[in]    xdata      : one dimensinal array of floats
  !! @param[in]    ydata      : one dimensinal array of floats
  !! @param[in]    normalize  : whether to normalize date before fitting
  !! @return       a data type containing fitting results
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function fit_least_squares(xdata, ydata, normalize) result(res)

    ! Formulas were taken from
    ! http://mathworld.wolfram.com/LeastSquaresFitting.html

    real(wp), dimension(:),           intent(in) :: xdata
    real(wp), dimension(size(xdata)), intent(in) :: ydata
    logical, optional                            :: normalize

    type(s_fitting_result) :: res

    logical  :: norm
    integer  :: i, n
    real(wp) :: x_mean, y_mean, x_min, x_max, y_min, y_max, x_scale, y_scale
    real(wp) :: ss_xx, ss_xy, ss_yy, a, b, r2, s, se_a, se_b

    real(wp), dimension(:), allocatable :: xwork, ywork


    norm = .false.
    if (present(normalize)) then
      norm = normalize
    end if

    n = size(xdata)

    allocate(xwork(n), ywork(n))

    if (norm) then

      x_min = minval(xdata)
      x_max = maxval(xdata)
      x_scale = x_max - x_min
      y_min = minval(ydata)
      y_max = maxval(ydata)
      y_scale = y_max - y_min

      xwork = (xdata - x_min) / x_scale
      ywork = (ydata - y_min) / y_scale

    else

      xwork = xdata
      ywork = ydata

    end if

    x_mean = sum(xwork) / n
    y_mean = sum(ywork) / n

    ss_xx = sum([ ((xwork(i)-x_mean)**2,i=1,n) ])
    ss_yy = sum([ ((ywork(i)-y_mean)**2,i=1,n) ])
    ss_xy = sum([ ((xwork(i)-x_mean)*(ywork(i)-y_mean),i=1,n) ])

    b = ss_xy / ss_xx
    a = y_mean - b  * x_mean

    r2 = ss_xy**2/(ss_xx*ss_yy)

    s = sqrt((ss_yy - b * ss_xy)/(n - 2))

    se_b = s / sqrt(ss_xx)

    if (norm) then

      b = b * y_scale / x_scale
      a = y_scale * a - b * x_min + y_min

      se_b = se_b * y_scale / x_scale
      se_a = s * y_scale &
        * sqrt(1_wp / n + (x_scale * x_mean + x_min)**2 / (x_scale)**2 / ss_xx)
    else
      se_a = s * sqrt(1_wp / n + x_mean**2/ss_xx)
    end if

    res%coeff = [a, b]
    res%stderr = [se_a, se_b]
    res%corr = sqrt(r2)

    return

  end function fit_least_squares

end module diffusion_impl_mod
