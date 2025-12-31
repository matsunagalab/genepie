!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   result_sink_mod
!> @brief   Abstract result sink for unified analysis output
!! @authors Claude Code
!
!  (c) Copyright 2024 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module result_sink_mod

  use, intrinsic :: iso_c_binding
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! Sink type enumeration
  integer, parameter, public :: SINK_FILE  = 1
  integer, parameter, public :: SINK_ARRAY = 2

  ! Abstract result sink type
  type, public :: s_result_sink
    integer :: sink_type = SINK_FILE

    ! File-based sink (CLI mode)
    integer            :: file_unit = -1
    logical            :: file_open = .false.

    ! Array-based sink (Python mode - zerocopy)
    real(wp), pointer  :: results(:) => null()
    integer            :: array_size = 0

    ! Common
    integer            :: current_index = 0
  end type s_result_sink

  ! Public subroutines
  public :: init_sink_file
  public :: init_sink_array
  public :: write_result
  public :: write_result_with_index
  public :: get_result_count
  public :: finalize_sink

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_sink_file
  !> @brief        Initialize file-based result sink (CLI mode)
  !! @authors      Claude Code
  !! @param[inout] sink     : result sink
  !! @param[in]    filename : output file name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_sink_file(sink, filename)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink
    character(*),        intent(in)    :: filename

    sink%sink_type = SINK_FILE
    sink%current_index = 0

    if (len_trim(filename) > 0) then
      call open_file(sink%file_unit, filename, IOFileOutputNew)
      sink%file_open = .true.
    else
      sink%file_open = .false.
    end if

    return

  end subroutine init_sink_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_sink_array
  !> @brief        Initialize array-based result sink (Python zerocopy mode)
  !! @authors      Claude Code
  !! @param[inout] sink       : result sink
  !! @param[in]    results    : pre-allocated results array pointer
  !! @param[in]    array_size : size of results array
  !! @note         The results array must be pre-allocated by caller (Python).
  !!               Fortran writes directly to this array (zerocopy).
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_sink_array(sink, results, array_size)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink
    real(wp), target,    intent(in)    :: results(:)
    integer,             intent(in)    :: array_size

    sink%sink_type = SINK_ARRAY
    sink%results => results
    sink%array_size = array_size
    sink%current_index = 0

    return

  end subroutine init_sink_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_result
  !> @brief        Write a single result value (auto-increments index)
  !! @authors      Claude Code
  !! @param[inout] sink  : result sink
  !! @param[in]    value : result value to write
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_result(sink, value)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink
    real(wp),            intent(in)    :: value

    sink%current_index = sink%current_index + 1
    call write_result_with_index(sink, sink%current_index, value)

    return

  end subroutine write_result

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_result_with_index
  !> @brief        Write a result value with explicit index
  !! @authors      Claude Code
  !! @param[inout] sink  : result sink
  !! @param[in]    idx   : result index (1-indexed)
  !! @param[in]    value : result value to write
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_result_with_index(sink, idx, value)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink
    integer,             intent(in)    :: idx
    real(wp),            intent(in)    :: value

    select case(sink%sink_type)

    case(SINK_FILE)
      if (sink%file_open) then
        write(sink%file_unit, '(i10,1x,f10.5)') idx, value
      end if

    case(SINK_ARRAY)
      if (idx >= 1 .and. idx <= sink%array_size) then
        sink%results(idx) = value
      else
        call error_msg('write_result_with_index> Index out of bounds')
      end if

    end select

    return

  end subroutine write_result_with_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_result_count
  !> @brief        Get number of results written
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_result_count(sink) result(count)

    ! formal arguments
    type(s_result_sink), intent(in) :: sink
    integer :: count

    count = sink%current_index

    return

  end function get_result_count

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    finalize_sink
  !> @brief        Finalize and cleanup result sink
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine finalize_sink(sink)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink

    select case(sink%sink_type)

    case(SINK_FILE)
      if (sink%file_open) then
        call close_file(sink%file_unit)
        sink%file_open = .false.
      end if
      sink%file_unit = -1

    case(SINK_ARRAY)
      sink%results => null()
      sink%array_size = 0

    end select

    sink%current_index = 0

    return

  end subroutine finalize_sink

end module result_sink_mod
