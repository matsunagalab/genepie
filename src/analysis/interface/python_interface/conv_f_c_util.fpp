!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> @brief   utilities for convert fortran and c lang
!! @authors imsbio
!
!--------1---------2---------3---------4---------5---------6---------7---------8
module conv_f_c_util
  use, intrinsic :: iso_c_binding
  implicit none
  private

  public :: c2f_string
  public :: f2c_bool_array
  public :: f2c_int_array
  public :: f2c_string_array
  public :: f2c_double_array

  public :: f2c_bool_array_nullcheck
  public :: f2c_int_array_nullcheck
  public :: f2c_string_array_nullcheck
  public :: f2c_double_array_nullcheck

contains
  ! Helper function to convert C string to Fortran string
  subroutine c2f_string(c_string, f_string)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char), intent(in) :: c_string(*)
    character(*), intent(out) :: f_string
    integer :: i

    f_string = ' '  ! Use a space instead of an empty string
    do i = 1, len(f_string)
      if (c_string(i) == c_null_char) exit
      f_string(i:i) = c_string(i)
    end do
  end subroutine c2f_string

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !> @brief        copy src fortran string array
  !!               to allocated c lang char array
  !! @authors      imsbio
  !! @param[in]    f_src : fortran string array
  !! @return       allocated c lang 2 dimensional char array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  type(c_ptr) function f2c_string_array(f_src) result(c_dst)
    implicit none
    character(*), intent(in) :: f_src(:)
    character(kind=c_char), pointer :: buf(:,:)
    integer :: i, j
    allocate(buf(len(f_src(lbound(f_src, 1))), &
                 lbound(f_src, 1):ubound(f_src, 1)))
    do i = lbound(f_src, 1), ubound(f_src, 1)
      do j = 1, len(f_src(i))
        buf(j,i) = f_src(i)(j:j)
      end do
    end do
    c_dst = c_loc(buf)
  end function f2c_string_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !> @brief        copy src fortran string array
  !!               to allocated c lang char array
  !! @authors      imsbio
  !! @param[in]    f_src : allocatable fortran string array
  !! @return       allocated c lang 2 dimensional char array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  type(c_ptr) function f2c_string_array_nullcheck(f_src) result(c_dst)
    implicit none
    character(*), intent(in), allocatable :: f_src(:)
    if (allocated(f_src)) then
        c_dst = f2c_string_array(f_src)
    else
        c_dst = c_null_ptr
    end if
  end function f2c_string_array_nullcheck

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !> @brief        copy src fortran integer array
  !!               to allocated c lang integer array
  !! @authors      imsbio
  !! @param[in]    f_src : fortran integer array
  !! @return       allocated c lang integer array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  type(c_ptr) function f2c_bool_array(f_src) result(c_dst)
    implicit none
    logical, intent(in) :: f_src(:)
    logical(c_bool), pointer :: buf(:)
    allocate(buf(size(f_src)))
    buf(1:size(f_src)) = f_src(lbound(f_src, 1):ubound(f_src, 1))
    c_dst = c_loc(buf)
  end function f2c_bool_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !> @brief        copy src fortran integer array
  !!               to allocated c lang integer array
  !! @authors      imsbio
  !! @param[in]    f_src : allocatable fortran integer array
  !! @return       allocated c lang integer array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  type(c_ptr) function f2c_bool_array_nullcheck(f_src) result(c_dst)
    implicit none
    logical, intent(in), allocatable :: f_src(:)
    if (allocated(f_src)) then
        c_dst = f2c_bool_array(f_src)
    else
        c_dst = c_null_ptr
    end if
  end function f2c_bool_array_nullcheck

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !> @brief        copy src fortran integer array
  !!               to allocated c lang integer array
  !! @authors      imsbio
  !! @param[in]    f_src : fortran integer array
  !! @return       allocated c lang integer array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  type(c_ptr) function f2c_int_array(f_src) result(c_dst)
    use constants_mod
    implicit none
    integer, intent(in) :: f_src(..)
    integer(c_int), pointer :: buf(:)
    allocate(buf(size(f_src)))
    select rank(f_src)
    rank(1)
        buf = reshape(f_src, [size(f_src)])
    rank(2)
        buf = reshape(f_src, [size(f_src)])
    rank(3)
        buf = reshape(f_src, [size(f_src)])
    rank(4)
        buf = reshape(f_src, [size(f_src)])
    rank default
        stop
    end select
    c_dst = c_loc(buf)
  end function f2c_int_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !> @brief        copy src fortran integer array
  !!               to allocated c lang integer array
  !! @authors      imsbio
  !! @param[in]    f_src : allocatable fortran integer array
  !! @return       allocated c lang integer array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  type(c_ptr) function f2c_int_array_nullcheck(f_src) result(c_dst)
    implicit none
    integer, intent(in), allocatable :: f_src(..)
    if (allocated(f_src)) then
        c_dst = f2c_int_array(f_src)
    else
        c_dst = c_null_ptr
    end if
  end function f2c_int_array_nullcheck

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !> @brief        copy src fortran double array
  !!               to allocated c lang double array
  !! @authors      imsbio
  !! @param[in]    f_src : fortran double array
  !! @return       allocated c lang double array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  type(c_ptr) function f2c_double_array(f_src) result(c_dst)
    use constants_mod
    implicit none
    real(wp), intent(in) :: f_src(..)
    real(c_double), pointer :: buf(:)
    allocate(buf(size(f_src)))
    select rank(f_src)
    rank(1)
        buf = reshape(f_src, [size(f_src)])
    rank(2)
        buf = reshape(f_src, [size(f_src)])
    rank(3)
        buf = reshape(f_src, [size(f_src)])
    rank(4)
        buf = reshape(f_src, [size(f_src)])
    rank default
        stop
    end select
    c_dst = c_loc(buf)
  end function f2c_double_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !> @brief        copy src fortran double array
  !!               to allocated c lang double array
  !! @authors      imsbio
  !! @param[in]    f_src : allocatable fortran double array
  !! @return       allocated c lang double array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  type(c_ptr) function f2c_double_array_nullcheck(f_src) result(c_dst)
    use constants_mod
    implicit none
    real(wp), intent(in), allocatable :: f_src(..)
    if (allocated(f_src)) then
        c_dst = f2c_double_array(f_src)
    else
        c_dst = c_null_ptr
    end if
  end function f2c_double_array_nullcheck
end module
