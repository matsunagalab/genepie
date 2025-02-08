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
  public :: c2f_string_allocate
  public :: f2c_bool_array
  public :: f2c_int_array
  public :: f2c_string_array
  public :: f2c_double_array

  public :: f2c_bool_array_nullcheck
  public :: f2c_int_array_nullcheck
  public :: f2c_string_array_nullcheck
  public :: f2c_double_array_nullcheck

  public :: allocate_c_bool_array
  public :: allocate_c_int_array
  public :: allocate_c_double_array
  public :: allocate_c_str_array

  public :: c2f_bool_array
  public :: c2f_int_array
  public :: c2f_int_array_static
  public :: c2f_double_array
  public :: c2f_string_array

  public :: deallocate_double
  public :: deallocate_double2

contains
  ! Helper function to convert C string to Fortran string
  subroutine c2f_string(c_string, f_string)
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

  subroutine c2f_string_allocate(c_string, f_string, limit_len)
    implicit none
    character(kind=c_char), intent(in) :: c_string(*)
    character(:), intent(inout), allocatable :: f_string
    integer, intent(in), optional :: limit_len
    integer :: i
    integer :: len_str
    integer :: limit

    if (.not. present(limit_len)) then
      limit = limit_len
    else
      limit = 10000000
    end if

    len_str = max(len_c_str(c_string, limit), 1) ! Use a space instead of an empty string
    if (allocated(f_string)) then
      if (len(f_string) < len_str) then
        deallocate(f_string)
        allocate(character(len_str) :: f_string)
      end if
    else
      allocate(character(len_str) :: f_string)
    end if
    call c2f_string(c_string, f_string)
  end subroutine c2f_string_allocate

  integer function len_c_str(c_string, limit_len) result(len_str)
    implicit none
    character(kind=c_char), intent(in) :: c_string(*)
    integer, intent(in), optional :: limit_len
    integer :: limit
    integer :: i

    if (.not. present(limit_len)) then
      limit = limit_len
    else
      limit = 10000000
    end if
    len_str = 0
    do i = 1, limit
      if (c_string(i) == c_null_char) then
        len_str = max(i - 1, 1)
        exit
      end if
    end do
  end function len_c_str

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

  !> foobar
  type(c_ptr) function allocate_c_bool_array(size) result(c_dst)
    implicit none
    integer, intent(in) :: size
    logical(c_bool), pointer :: buf(:)

    allocate(buf(size))
    c_dst = c_loc(buf)
  end function allocate_c_bool_array

  type(c_ptr) function allocate_c_int_array(size) result(c_dst)
    implicit none
    integer, intent(in) :: size
    integer(c_int), pointer :: buf(:)

    allocate(buf(size))
    c_dst = c_loc(buf)
  end function allocate_c_int_array

  type(c_ptr) function allocate_c_double_array(size) result(c_dst)
    implicit none
    integer, intent(in) :: size
    real(c_double), pointer :: buf(:)

    allocate(buf(size))
    c_dst = c_loc(buf)
  end function allocate_c_double_array

  type(c_ptr) function allocate_c_str_array(size, len_str) result(c_dst)
    implicit none
    integer, intent(in) :: size
    integer, intent(in) :: len_str
    character(kind=c_char), pointer :: buf(:,:)

    allocate(buf(len_str, size))
    c_dst = c_loc(buf)
  end function allocate_c_str_array

  subroutine c2f_bool_array(f_dst, c_src, v_shape)
    implicit none
    logical, intent(out), allocatable :: f_dst(..)
    type(c_ptr), intent(in) :: c_src
    integer, intent(in) :: v_shape(:)
    logical(c_bool), pointer :: buf(:)

    select rank(f_dst)
    rank(1)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        allocate(f_dst(v_shape(1)))
        f_dst = reshape(buf, [v_shape(1)])
    rank(2)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        allocate(f_dst(v_shape(1), v_shape(2)))
        f_dst = reshape(buf, [v_shape(1), v_shape(2)])
    rank(3)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        allocate(f_dst(v_shape(1), v_shape(2), v_shape(3)))
        f_dst = reshape(buf, [v_shape(1), v_shape(2), v_shape(3)])
    rank(4)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        allocate(f_dst(v_shape(1), v_shape(2), v_shape(3), v_shape(4)))
        f_dst = reshape(buf, [v_shape(1), v_shape(2), v_shape(3), v_shape(4)])
    rank default
        stop
    end select
  end subroutine c2f_bool_array

  subroutine c2f_int_array(f_dst, c_src, v_shape)
    implicit none
    integer, intent(out), allocatable :: f_dst(..)
    type(c_ptr), intent(in) :: c_src
    integer, intent(in) :: v_shape(:)
    integer(c_int), pointer :: buf(:)

    select rank(f_dst)
    rank(1)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        allocate(f_dst(v_shape(1)))
        f_dst = reshape(buf, [v_shape(1)])
    rank(2)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        allocate(f_dst(v_shape(1), v_shape(2)))
        f_dst = reshape(buf, [v_shape(1), v_shape(2)])
    rank(3)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        allocate(f_dst(v_shape(1), v_shape(2), v_shape(3)))
        f_dst = reshape(buf, [v_shape(1), v_shape(2), v_shape(3)])
    rank(4)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        allocate(f_dst(v_shape(1), v_shape(2), v_shape(3), v_shape(4)))
        f_dst = reshape(buf, [v_shape(1), v_shape(2), v_shape(3), v_shape(4)])
    rank default
        stop
    end select
  end subroutine c2f_int_array

  subroutine c2f_int_array_static(f_dst, c_src, v_shape)
    implicit none
    integer, intent(out) :: f_dst(..)
    type(c_ptr), intent(in) :: c_src
    integer, intent(in) :: v_shape(:)
    integer(c_int), pointer :: buf(:)

    select rank(f_dst)
    rank(1)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        f_dst = reshape(buf, [v_shape(1)])
    rank(2)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        f_dst = reshape(buf, [v_shape(1), v_shape(2)])
    rank(3)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        f_dst = reshape(buf, [v_shape(1), v_shape(2), v_shape(3)])
    rank(4)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        f_dst = reshape(buf, [v_shape(1), v_shape(2), v_shape(3), v_shape(4)])
    rank default
        stop
    end select
  end subroutine c2f_int_array_static

  subroutine c2f_double_array(f_dst, c_src, v_shape)
    use constants_mod
    implicit none
    real(wp), intent(out), allocatable :: f_dst(..)
    type(c_ptr), intent(in) :: c_src
    integer, intent(in) :: v_shape(:)
    real(c_double), pointer :: buf(:)

    select rank(f_dst)
    rank(1)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        allocate(f_dst(v_shape(1)))
        f_dst = reshape(buf, [v_shape(1)])
    rank(2)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        allocate(f_dst(v_shape(1), v_shape(2)))
        f_dst = reshape(buf, [v_shape(1), v_shape(2)])
    rank(3)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        allocate(f_dst(v_shape(1), v_shape(2), v_shape(3)))
        f_dst = reshape(buf, [v_shape(1), v_shape(2), v_shape(3)])
    rank(4)
        call C_F_pointer(c_src, buf, [product(v_shape)])
        allocate(f_dst(v_shape(1), v_shape(2), v_shape(3), v_shape(4)))
        f_dst = reshape(buf, [v_shape(1), v_shape(2), v_shape(3), v_shape(4)])
    rank default
        stop
    end select
  end subroutine c2f_double_array

  subroutine c2f_string_array(f_dst, c_src, size)
    use constants_mod
    implicit none
    character(*), intent(out), allocatable :: f_dst(:)
    type(c_ptr), intent(in) :: c_src
    integer, intent(in) :: size
    character(kind=c_char), pointer :: buf(:,:)
    integer :: i, j

    allocate(f_dst(size))
    call C_F_pointer(c_src, buf, [len(f_dst(0)), size])
    do i = 1, size
      do j = 1, len(f_dst(0))
        f_dst(i)(j:j) = buf(j,i)
      end do
    end do
  end subroutine c2f_string_array

  subroutine deallocate_double(c_src, size) &
          bind(C, name="deallocate_double")
    implicit none
    type(c_ptr), intent(in) :: c_src
    integer(c_int), intent(in) :: size
    real(c_double), pointer :: buf(:)

    call C_F_pointer(c_src, buf, [size])
    deallocate(buf)
  end subroutine

  subroutine deallocate_double2(c_src, size1, size2) &
          bind(C, name="deallocate_double2")
    implicit none
    type(c_ptr), intent(in) :: c_src
    integer(c_int), intent(in) :: size1
    integer(c_int), intent(in) :: size2
    real(c_double), pointer :: buf(:,:)

    call C_F_pointer(c_src, buf, [size2, size1])
    deallocate(buf)
  end subroutine
end module
