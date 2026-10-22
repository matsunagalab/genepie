!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ea_analyze_mod
!> @brief   run analyzing pca file
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ea_analyze_mod

  use output_str_mod
  use input_str_mod
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  public  :: analyze_eigmat_unified

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_eigmat_unified
  !> @brief        eigen decomposition of a variance-covariance matrix
  !!               (shared by the CLI and the Python interface)
  !! @authors      NT, Claude Code
  !! @param[in]    vcv_matrix : symmetric variance-covariance matrix (n, n)
  !! @param[in]    n          : matrix dimension (3 * number of atoms)
  !! @param[out]   eigval     : eigenvalues in descending order (n)
  !! @param[out]   eigvec     : eigenvectors (n, n); eigvec(j, k) is the j-th
  !!                            component of the k-th principal mode
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_eigmat_unified(vcv_matrix, n, eigval, eigvec)

    ! formal arguments
    real(wp),                intent(in)    :: vcv_matrix(:,:)
    integer,                 intent(in)    :: n
    real(wp),                intent(out)   :: eigval(:)
    real(wp),                intent(out)   :: eigvec(:,:)

    ! local variables
    integer                  :: i, j, k, alloc_stat
    integer                  :: INFO, LDA, LWORK
    character(1)             :: JOBZ, UPLO

    real(wp), allocatable    :: a(:,:), w(:), work(:)


    if (n < 1) &
      call error_msg('Analyze_Eigmat_Unified> matrix dimension must be positive')
    if (size(vcv_matrix,1) < n .or. size(vcv_matrix,2) < n) &
      call error_msg('Analyze_Eigmat_Unified> vcv_matrix is too small')
    if (size(eigval) < n .or. size(eigvec,1) < n .or. size(eigvec,2) < n) &
      call error_msg('Analyze_Eigmat_Unified> eigval/eigvec are too small')

    LWORK = max(1, 3*n-1)

    allocate(a(n,n), w(n), work(LWORK), stat=alloc_stat)
    if (alloc_stat /= 0) &
       call error_msg_alloc

    do i = 1, n
      do j = 1, n
        a(i, j) = vcv_matrix(i, j)
      end do
    end do

    ! DSYEV is LAPACK subroutine
    !
    JOBZ  = 'V'
    UPLO  = 'U'
    LDA   = n
    INFO  = 0

#ifdef LAPACK
    call dsyev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
#else
    call error_msg('Analyze> ERROR: This subroutine needs LAPACK.')
#endif

    if (INFO /= 0) &
      call error_msg('Analyze_Eigmat_Unified> DSYEV failed to converge')

    ! DSYEV returns ascending eigenvalues: reverse to descending order
    !
    k = 1
    do i = n, 1, -1
      eigval(k) = w(i)
      do j = 1, n
        eigvec(j, k) = a(j, i)
      end do
      k = k + 1
    end do

    deallocate(work, w, a)

    return

  end subroutine analyze_eigmat_unified

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing pca file (CLI driver)
  !! @authors      NT, Claude Code
  !! @param[inout] input      : input information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(input, output)

    ! formal arguments
    type(s_input),           intent(inout) :: input
    type(s_output),          intent(inout) :: output

    ! local variables
    real(wp)                 :: rmsf, msf,  acu, acu2, r, s
    integer                  :: pca_in, val_out, vec_out, cnt_out
    integer                  :: i, j, k, natom, nvcv, alloc_stat

    real(wp), allocatable    :: vcv_matrix(:,:)
    real(wp), allocatable    :: eigval(:), eigvec(:,:)


    ! read PCA file
    !
    call open_file(pca_in, input%pcafile, IOFileInput)

    rmsf = 0.0_wp
    read(pca_in,*) natom

    nvcv = natom * 3

    allocate(vcv_matrix(nvcv, nvcv), stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    do i = 1, nvcv
      do j = 1, nvcv
        read(pca_in,*) vcv_matrix(i,j)
      end do
      rmsf = rmsf + vcv_matrix(i,i)
    end do

    write(MsgOut,'(a,i12)')    'INFO> natom = ', natom
    write(MsgOut,'(a,2f10.5)') '      rmsf  = ', rmsf, sqrt(rmsf/real(natom,wp))
    write(MsgOut,'(a)') ''


    ! calculate eigen-value
    !
    call open_file(val_out, output%valfile, IOFileOutputNew)
    call open_file(vec_out, output%vecfile, IOFileOutputNew)
    call open_file(cnt_out, output%cntfile, IOFileOutputNew)

    allocate(eigval(nvcv),       &
             eigvec(nvcv, nvcv), stat=alloc_stat)
    if (alloc_stat /= 0) &
       call error_msg_alloc

    call analyze_eigmat_unified(vcv_matrix, nvcv, eigval, eigvec)

    do k = 1, nvcv
      write(val_out, *) eigval(k)
    end do

    do k = 1, nvcv
      do j = 1, nvcv
        ! eigvec(j, k) j-th component of k-th PC mode
        write(vec_out, *) eigvec(j, k)
      end do
    end do

    ! calculate the contribution of each mode 
    !
    msf = 0.0_wp
    do i = 1, nvcv-7
      msf = msf + eigval(i)
    end do

    j    = 0
    acu  = 0.0_wp
    acu2 = 0.0_wp

    do i = 1, nvcv-7
      j   = j + 1
      r   = eigval(i)
      s   = eigval(i)/msf
      acu = acu + s
      if (r > 0) then
        write(cnt_out, '(i5,1x,3(f10.3,1x))') &
               j, sqrt(r/natom), s*100.0_wp, acu*100.0_wp
      else
        write(cnt_out, '(i5,1x,3(f10.3,1x))') &
               j, sqrt(-r/natom), s*100.0_wp, acu*100.0_wp
      end if
    end do


    write(MsgOut,'(a)')        'Analyze> RMSF'
    write(MsgOut,'(a,2f10.5)') '         rmsf = ', msf, sqrt(msf/real(natom,wp))
    write(MsgOut,'(a)')        ''

    call close_file(val_out)
    call close_file(vec_out)
    call close_file(cnt_out)

    ! Output summary
    !
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [valfile] ' // trim(output%valfile)
    write(MsgOut,'(A)') '    Column 1: Eigen values of the variance-covariance matrix'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [valfile] ' // trim(output%vecfile)
    write(MsgOut,'(A)') '    Column 1: Eigen vectors of the variance-covariance matrix'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [cntfile] ' // trim(output%cntfile)
    write(MsgOut,'(A)') '    Column 1: Mode index'
    write(MsgOut,'(A)') '    Column 2: sqrt(eigval/n)'
    write(MsgOut,'(A)') '    Column 3: Contribution of each mode to the total fluctuation (%)'
    write(MsgOut,'(A)') '    Column 4: Accumulated contribution (%)'
    write(MsgOut,'(A)') ''

    ! deallocate memory
    !
    deallocate(eigvec, eigval)
    deallocate(vcv_matrix)

    return

  end subroutine analyze

end module ea_analyze_mod
