!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pmf_impl_mod
!> @brief   PMF analysis that returns the result in memory (Python interface)
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pmf_impl_mod

  use pm_option_str_mod
  use output_str_mod
  use input_str_mod
  use fileio_mod
  use error_mod
  use string_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! constants
  real(wp),        parameter   :: KB   = 0.00198719168260038_wp

  ! subroutines
  public  :: analyze
  private :: calc_pmf1d
  private :: calc_pmf2d
  private :: check_file_lines
  private :: check_file_column
  private :: get_replicate_name1
  private :: periodic

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        compute the PMF and return it in memory
  !! @authors      NT
  !! @param[in]    input   : input information
  !! @param[in]    output  : output information
  !! @param[in]    option  : option information
  !! @param[out]   pmf     : result array. 1D -> (3, nbin) with rows
  !!                         (center, standard PMF, Gaussian PMF); 2D ->
  !!                         (nbin_y, nbin_x) free energy matrix
  !! @param[out]   n_out1  : first  (slow) dimension exported to Python
  !! @param[out]   n_out2  : second (fast) dimension exported to Python
  !! @param[inout] err     : error information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(input, output, option, pmf, n_out1, n_out2, err)

    ! formal arguments
    type(s_input),           intent(in)    :: input
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(in)    :: option
    real(wp), pointer,       intent(out)   :: pmf(:,:)
    integer,                 intent(out)   :: n_out1
    integer,                 intent(out)   :: n_out2
    type(s_error),           intent(inout) :: err

    ! local variables
    real(wp)                 :: w, d1, d2, distance
    integer                  :: ireplica, i, j, file_w, file_d, file_dist
    integer                  :: nline, nline_w, nline_d, nline_dist, ncol_d
    real(wp)                 :: tim

    real(wp),    allocatable :: weight(:), data(:, :), pmf1d(:,:), pmf2d(:, :)


    if (option%check_only) return
    if (error_has(err)) return

    ! setup data
    !
    call check_file_lines(get_replicate_name1(input%weightfile, 1), nline_w)
    call check_file_lines(get_replicate_name1(input%cvfile, 1), nline_d)
    call check_file_column(get_replicate_name1(input%cvfile, 1), ncol_d)

    if (nline_w /= 0 .and. nline_w /= nline_d) then
      call error_set(err, ERROR_DATA_MISMATCH, &
        'Analyze> # of weight file lines is different from cv file lines')
      return
    end if

    if (ncol_d < (option%dimension+1)) then
      call error_set(err, ERROR_DATA_MISMATCH, &
        'Analyze> # of column of cv file must be >= dimension+1')
      return
    end if

    if (input%distfile /= '') then
      call check_file_lines(get_replicate_name1(input%distfile, 1), nline_dist)

      if (nline_w /= 0 .and. nline_w /= nline_dist) then
        call error_set(err, ERROR_DATA_MISMATCH, &
          'Analyze> # of weight file lines is different from distance file lines')
        return
      end if
    end if

    allocate(weight  (nline_d * option%nreplicas), &
             data(option%dimension, nline_d * option%nreplicas))

    nline = 0
    do ireplica = 1, option%nreplicas

      file_w = 0
      file_d = 0
      file_dist = 0

      ! open weight file
      if (input%weightfile /= '') &
        call open_file(file_w, get_replicate_name1(input%weightfile, ireplica), &
                       IOFileInput)

      ! open cv file
      call open_file(file_d, get_replicate_name1(input%cvfile, ireplica), &
                     IOFileInput)

      ! open distance file
      if (input%distfile /= '') &
        call open_file(file_dist, get_replicate_name1(input%distfile, ireplica), &
                       IOFileInput)

      do i = 1, nline_d

        ! read weight file
        if (file_w /= 0) &
          read(file_w,*) tim, w

        ! read cv file
        if (option%dimension == 1) then
          read(file_d,*) tim, d1
          if (file_dist /= 0) then
            read(file_dist,*) tim, distance
            if (option%cutoff == 0.0_wp .or. distance < option%cutoff) then
              nline = nline + 1
              weight(nline) = w
              data(1, nline) = d1
            end if
          else
            nline = nline + 1
            weight(nline) = w
            data(1, nline) = d1
          end if
        else if (option%dimension == 2) then
          read(file_d,*) tim, d1, d2
          if (file_dist /= 0) then
            read(file_dist,*) tim, distance
            if (option%cutoff == 0.0_wp .or. distance < option%cutoff) then
              nline = nline + 1
              weight(nline) = w
              data(1, nline) = d1
              data(2, nline) = d2
            end if
          else
            nline = nline + 1
            weight(nline) = w
            data(1, nline) = d1
            data(2, nline) = d2
          end if
        end if

      end do

      ! close cv file
      call close_file(file_d)

      ! close weight file
      if (file_w /= 0) &
        call close_file(file_w)

      ! close distance file
      if (file_dist /= 0) &
        call close_file(file_dist)

    end do

    if (input%weightfile == '') &
      weight(1:nline) = 1.0_wp / real(nline,wp)


    ! calc pmf and store into the result array
    !
    if (option%dimension == 1) then
      call calc_pmf1d(option, data(:, 1:nline), weight(1:nline), pmf1d)

      n_out1 = size(pmf1d(1,:))    ! number of bins (Python rows)
      n_out2 = 3                   ! center, standard PMF, Gaussian PMF
      allocate(pmf(n_out2, n_out1))
      do i = 1, n_out1
        pmf(1, i) = option%center(1, i)
        pmf(2, i) = pmf1d(1, i)
        pmf(3, i) = pmf1d(2, i)
      end do

      deallocate(pmf1d)

    else if (option%dimension == 2) then
      call calc_pmf2d(option, data(:, 1:nline), weight(1:nline), pmf2d)

      n_out1 = size(pmf2d(:, 1))   ! nbin_x (Python rows)
      n_out2 = size(pmf2d(1, :))   ! nbin_y (Python cols)
      allocate(pmf(n_out2, n_out1))
      do i = 1, n_out1
        do j = 1, n_out2
          pmf(j, i) = pmf2d(i, j)
        end do
      end do

      deallocate(pmf2d)

    end if

    deallocate(weight, data)

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_pmf1d(option, data, weight, pmf)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    real(wp),                intent(in)    :: data(:, :)
    real(wp),                intent(in)    :: weight(:)
    real(wp),                allocatable   :: pmf(:,:)
    real(wp)                               :: KBT

    ! local variables
    integer                  :: nbin, ndata
    integer                  :: ibin, idata
    integer                  :: ipmf

    real(wp)                 :: ddata = 0.0_wp

    real(wp),    allocatable :: dx1(:,:)


    KBT = KB * option%temperature

    nbin  = size(option%center(1, :))
    ndata = size(data(1, :))

    allocate(pmf(2, nbin))
    allocate(dx1(nbin, ndata))

    ! ipmf = 1 > Standard PMF, 2 > Gaussian type PMF
    do ipmf = 1, 2
      dx1(1:nbin, 1:ndata) = 0.0_wp

      do ibin = 1, nbin
        if (option%is_periodic(1)) then
          if(ipmf == 2) then
            do idata = 1, ndata
              dx1(ibin, idata) = (periodic(data(1, idata), option%center(1, ibin), &
                option%box_size(1)) / option%band_width(1))**2.0_wp
            end do
          else
            do idata = 1, ndata
              ddata = data(1, idata) - option%center(1, ibin)
              if((-option%delta_grid(1) < ddata) .and. (option%delta_grid(1) > ddata)) &
                dx1(ibin, idata) = dx1(ibin, idata) + 1.0_wp
            end do
          end if
        else
          if(ipmf == 2) then
            dx1(ibin, 1:ndata) = ((data(1, 1:ndata) - option%center(1, ibin)) / &
              option%band_width(1)) ** 2.0_wp
          else
            do idata = 1, ndata
              ddata = data(1, idata) - option%center(1, ibin)
              if((-option%delta_grid(1) < ddata) .and. (option%delta_grid(1) > ddata)) &
                dx1(ibin, idata) = dx1(ibin, idata) + 1.0_wp
            end do
          end if
        end if

        if(ipmf == 2) dx1(ibin, 1:ndata) = exp(-0.5_wp * dx1(ibin, 1:ndata)) / &
                                           (sqrt(2.0_wp*PI)*option%band_width(1))
        do idata = 1, ndata
          dx1(ibin, idata) = dx1(ibin, idata)*weight(idata)
        end do
      end do

      do ibin = 1, nbin
        pmf(ipmf, ibin) = 0.0_wp
        do idata = 1, ndata
          pmf(ipmf, ibin) = pmf(ipmf, ibin) + dx1(ibin, idata)
        end do
      end do

      do ibin = 1, nbin
        pmf(ipmf, ibin) = -KBT*log(pmf(ipmf, ibin))
      end do

      ! Lowest PMF is set to be 0.0
      pmf(ipmf, 1:nbin) = pmf(ipmf, 1:nbin) - minval(pmf(ipmf, 1:nbin))

    end do

    deallocate(dx1)

    return

  end subroutine calc_pmf1d

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_pmf2d(option, data, weight, pmf)

    ! formal arguments
    type(s_option),          intent(in)    :: option
    real(wp),                intent(in)    :: data(:, :)
    real(wp),                intent(in)    :: weight(:)
    real(wp),                allocatable   :: pmf(:, :)
    real(wp)                               :: KBT

    ! local variables
    integer                  :: nbin1, nbin2
    integer                  :: idata, ndata
    integer                  :: i, j

    real(wp),    allocatable :: dx1(:,:), dx2(:,:)


    KBT = KB * option%temperature

    nbin1  = option%num_grids(1)-1
    nbin2  = option%num_grids(2)-1
    ndata  = size(data(1, :))

    allocate(pmf(nbin1, nbin2))
    allocate(dx1(nbin1, ndata), dx2(nbin2, ndata))

    do i = 1, nbin1
      if (option%is_periodic(1)) then
        do idata = 1, ndata
          dx1(i, idata) = (periodic(data(1, idata), option%center(1, i), option%box_size(1)) &
            /option%band_width(1))**2.0_wp
        end do
      else
        dx1(i, 1:ndata) = ((data(1, 1:ndata) - option%center(1, i)) / &
          option%band_width(1)) ** 2.0_wp
      end if
      dx1(i, 1:ndata) = exp(-0.5_wp * dx1(i, 1:ndata)) / &
                           (sqrt(2.0_wp*PI)*option%band_width(1))
    end do

    do j = 1, nbin2
      if (option%is_periodic(2)) then
        do idata = 1, ndata
          dx2(j, idata) = (periodic(data(2, idata), option%center(2, j), option%box_size(2)) &
            /option%band_width(2))**2.0_wp
        end do
      else
        dx2(j, 1:ndata) = ((data(2, 1:ndata) - option%center(2, j)) / &
          option%band_width(2)) ** 2.0_wp
      end if
      dx2(j, 1:ndata) = exp(-0.5_wp * dx2(j, 1:ndata)) / &
                           (sqrt(2.0_wp*PI)*option%band_width(2))
    end do

    pmf(:, :) = 0.0_wp
    do i = 1, nbin1
      do j = 1, nbin2
        do idata = 1, ndata
          pmf(i, j) = pmf(i, j) + dx1(i, idata)*dx2(j, idata)*weight(idata)
        end do
      end do
    end do

    do i = 1, nbin1
      do j = 1, nbin2
        pmf(i, j) = -KBT*log(pmf(i, j))
      end do
    end do

    pmf(1:nbin1, 1:nbin2) = pmf(1:nbin1, 1:nbin2) - minval(pmf(1:nbin1, 1:nbin2))

    deallocate(dx1, dx2)

    return

  end subroutine calc_pmf2d

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_file_lines(filename, nline)

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(inout) :: nline

    ! local variables
    integer                  :: file
    character(100)           :: line


    if (filename == '') then
      nline = 0
      return
    end if

    call open_file(file, filename, IOFileInput)

    nline = 0
    do while(.true.)
      read(file,'(a)',end=10) line
      nline = nline + 1
    end do

10  call close_file(file)

    return

  end subroutine check_file_lines

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_file_column(filename, ncol)

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(inout) :: ncol

    ! local variables
    integer                  :: file
    character(1000)          :: line


    call open_file(file, filename, IOFileInput)

    read(file,'(a)') line
    ncol = split_num(line)

    call close_file(file)

    return

  end subroutine check_file_column

  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_replicate_name1(filename, no)

    ! return
    character(Maxfilename)   :: get_replicate_name1

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: no

    ! local variables
    integer                  :: bl, br


    if (filename == '') then
      get_replicate_name1 = ''
      return
    end if

    bl = index(filename, '{', back=.true.)
    br = index(filename, '}', back=.true.)

    if (bl == 0 .or. br == 0 .or. bl > br) &
      call error_msg('Get_Replicate_Name1> Syntax error.')

    write(get_replicate_name1, '(a,i0,a)') &
         filename(:bl-1),no,filename(br+1:len_trim(filename))

    return

  end function get_replicate_name1

  !======1=========2=========3=========4=========5=========6=========7=========8

  function periodic(x, center, box_size)

    ! return value
    real(wp)                 :: periodic

    ! formal arguments
    real(wp)                 :: x
    real(wp)                 :: center
    real(wp)                 :: box_size


    periodic = x - center
    periodic = periodic - nint(periodic/box_size)*box_size

    return

  end function periodic

end module pmf_impl_mod
