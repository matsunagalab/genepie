!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!  Module   messages_mod
!> @brief   utilities for taking messages to Std or Err output
!! @authors Yuji Sugita (YS)
!! @note    MsgOut and ErrOut are defined here.
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module messages_mod

  use, intrinsic :: iso_c_binding
  use mpi_parallel_mod

  implicit none
  private

  ! parameters
  integer, public, parameter :: MsgOut  = 6
  integer, public, parameter :: ErrOut  = 6

  ! Pending error handed to the Python interface when the setjmp guard is armed
  ! (see fileio_data_.c). error_msg records the message/code here and longjmps
  ! back to the bind(C) wrapper instead of aborting the process.
  character(len=:), allocatable, public, save :: fi_pending_msg
  integer,                       public, save :: fi_pending_code = 0

  ! C helpers (compiled into lib.a) implementing the library-mode error guard.
  interface
    function fi_error_is_armed() bind(C, name="fi_error_is_armed") result(armed)
      import :: c_int
      integer(c_int) :: armed
    end function fi_error_is_armed

    subroutine fi_error_signal() bind(C, name="fi_error_signal")
    end subroutine fi_error_signal
  end interface

  ! subroutines and functions
  public :: error_msg
  public :: error_msg_alloc
  public :: error_msg_dealloc
  public :: error_msg_fileio

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    error_msg
  !> @brief        write error message and stop execution
  !! @param[in]    message : Error message (optional)
  !! @authors      YS
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine error_msg(message, code)

#if defined(INTEL)
    use ifcore, only:tracebackqq

#elif defined(KCOMP)
    use service_routines, only:errtra

#endif

    ! formal arguments
    character(*), optional,  intent(in)    :: message
    integer,      optional,  intent(in)    :: code

    ! local variables
    character(1000)          :: msg


    if (present(message)) then
      write(msg,'(A,A12,I5)') trim(message), '  rank_no = ', my_world_rank
    else
      write(msg,'(A25,I5)') '               rank_no = ', my_world_rank
    end if

    ! Library mode (Python interface): record the error and unwind to the
    ! bind(C) wrapper instead of aborting. In CLI/MD-engine mode the guard is
    ! never armed, so execution falls through to the original abort below.
    if (fi_error_is_armed() /= 0) then
      if (present(message)) then
        fi_pending_msg = trim(message)
      else
        fi_pending_msg = 'GENESIS aborted'
      end if
      if (present(code)) then
        fi_pending_code = code
      else
        fi_pending_code = 699   ! ERROR_GENERIC (see error_mod)
      end if
      write(ErrOut,*) trim(msg)
      call fi_error_signal()    ! longjmp; does not return while armed
    end if

#if defined(INTEL)
    call TracebackQQ(trim(msg),-1)
    call exit(1)

#elif defined(KCOMP)
    write(ErrOut,*) trim(msg)
    call errtra
    call setrcd(1)
    call exit

#else
    write(ErrOut,*) trim(msg)
    call exit(1)

#endif

  end subroutine error_msg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    error_msg_alloc
  !> @brief        write error message in allocation and stop execution
  !! @param[in]    message : Error message
  !! @authors      YS
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine error_msg_alloc

    call error_msg('Memory allocation error', code=101)   ! ERROR_ALLOC

  end subroutine error_msg_alloc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    error_msg_dealloc
  !> @brief        write error message in deallocation and stop execution
  !! @param[in]    message : Error message
  !! @authors      YS
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine error_msg_dealloc

    call error_msg('Memory deallocation error', code=102)   ! ERROR_DEALLOC

  end subroutine error_msg_dealloc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    error_msg_fileio
  !> @brief        write error message in fileio and stop execution
  !! @param[in]    message : Error message
  !! @authors      YS
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine error_msg_fileio

    call error_msg('File I/O error', code=203)   ! ERROR_FILE_READ

  end subroutine error_msg_fileio

end module messages_mod
