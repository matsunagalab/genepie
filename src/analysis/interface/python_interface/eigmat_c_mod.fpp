!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   eigmat_c_mod
!> @brief   Python (bind(C)) entry point of eigmat_analysis
!! @authors Claude Code
!
!  The science lives in ea_analyze_mod (shared with the CLI program): the
!  entry point hands the caller's variance-covariance matrix to
!  analyze_eigmat_unified (LAPACK DSYEV) and fills the eigenvalue and
!  eigenvector buffers.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module eigmat_c_mod
  use, intrinsic :: iso_c_binding
  use ea_analyze_mod
  use error_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: eigmat_analysis_c

  !> State shared between eigmat_analysis_c and its guarded body.
  type :: t_eigmat_ctx
    type(c_ptr) :: vcv_ptr = c_null_ptr      ! (n, n) input
    integer     :: n = 0
    type(c_ptr) :: eigval_ptr = c_null_ptr   ! (n) output
    type(c_ptr) :: eigvec_ptr = c_null_ptr   ! (n, n) output
    type(s_error) :: err
  end type t_eigmat_ctx

  private :: t_eigmat_ctx

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    eigmat_analysis_c
  !> @brief        Eigen decomposition of a variance-covariance matrix
  !! @authors      Claude Code
  !! @param[in]    vcv_ptr    : symmetric (n, n) matrix
  !! @param[in]    n          : matrix dimension
  !! @param[in]    eigval_ptr : pre-allocated (n) eigenvalues, descending
  !! @param[in]    eigvec_ptr : pre-allocated (n, n) eigenvectors; element
  !!                            (j, k) is the j-th component of mode k
  !! @param[out]   status, msg : error status and message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine eigmat_analysis_c(vcv_ptr, n, eigval_ptr, eigvec_ptr, &
                               status, msg, msglen) &
        bind(C, name="eigmat_analysis_c")
    implicit none

    ! Arguments
    type(c_ptr), value :: vcv_ptr
    integer(c_int), value :: n
    type(c_ptr), value :: eigval_ptr
    type(c_ptr), value :: eigvec_ptr
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(t_eigmat_ctx), target :: c

    c%vcv_ptr    = vcv_ptr
    c%n          = n
    c%eigval_ptr = eigval_ptr
    c%eigvec_ptr = eigvec_ptr

    call run_guarded(eigmat_analysis_body, c, c%err, status, msg, msglen)

  end subroutine eigmat_analysis_c

  !> Guarded body of eigmat_analysis_c (see run_guarded in error_mod).
  subroutine eigmat_analysis_body(ctx) bind(C, name="eigmat_analysis_body")
    implicit none
    type(c_ptr), value :: ctx

    type(t_eigmat_ctx), pointer :: c
    real(wp), pointer :: vcv(:,:), eigval(:), eigvec(:,:)

    call c_f_pointer(ctx, c)

    if (c%n <= 0) then
      call error_set(c%err, ERROR_DIMENSION, &
                     "eigmat_analysis_c: matrix dimension must be positive")
      return
    end if
    if (.not. c_associated(c%vcv_ptr) .or. &
        .not. c_associated(c%eigval_ptr) .or. &
        .not. c_associated(c%eigvec_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "eigmat_analysis_c: array pointers are null")
      return
    end if

    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    call c_f_pointer(c%vcv_ptr, vcv, [c%n, c%n])
    call c_f_pointer(c%eigval_ptr, eigval, [c%n])
    call c_f_pointer(c%eigvec_ptr, eigvec, [c%n, c%n])

    write(MsgOut,'(A)') '[STEP1] Eigen decomposition (unified)'
    write(MsgOut,'(A)') ' '
    call analyze_eigmat_unified(vcv, c%n, eigval, eigvec)

  end subroutine eigmat_analysis_body

end module eigmat_c_mod
