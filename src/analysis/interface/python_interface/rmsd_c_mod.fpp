!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  ra_main
!! @brief   RMSD analysis
!! @authors Takaharu Mori (TM), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rmsd_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use rmsd_impl_mod
  use trj_source_mod

  use fitting_mod
  use fitting_str_mod
  use trajectory_str_mod
  use molecules_str_mod
  use string_mod
  use error_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: rmsd_analysis_c
  public :: rmsd_analysis_fitting_c
  public :: rmsd_analysis_lazy_c

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rmsd_analysis_c
  !> @brief        RMSD analysis (zerocopy, pre-allocated result array)
  !! @authors      Claude Code
  !! @param[in]    mass_ptr       : pointer to mass array (from Python NumPy)
  !! @param[in]    ref_coord_ptr  : pointer to reference coordinates (3, n_atoms)
  !! @param[in]    n_atoms        : number of atoms
  !! @param[in]    s_trajes_c     : trajectories C structure
  !! @param[in]    ana_period     : analysis period
  !! @param[in]    analysis_idx   : analysis atom indices (1-indexed)
  !! @param[in]    n_analysis     : number of analysis atoms
  !! @param[in]    mass_weighted  : use mass weighting (0 or 1)
  !! @param[in]    result_ptr     : pointer to pre-allocated result array
  !! @param[in]    result_size    : size of pre-allocated result array
  !! @param[out]   nstru_out      : number of frames analyzed
  !! @param[out]   status         : error status
  !! @param[out]   msg            : error message
  !! @param[in]    msglen         : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rmsd_analysis_c(mass_ptr, ref_coord_ptr, n_atoms, &
                             s_trajes_c, ana_period, &
                             analysis_idx, n_analysis, mass_weighted, &
                             result_ptr, result_size, nstru_out, &
                             status, msg, msglen) &
        bind(C, name="rmsd_analysis_c")
    implicit none

    ! Arguments
    type(c_ptr), value :: mass_ptr
    type(c_ptr), value :: ref_coord_ptr
    integer(c_int), value :: n_atoms
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer(c_int), value :: ana_period
    type(c_ptr), value :: analysis_idx
    integer(c_int), value :: n_analysis
    integer(c_int), value :: mass_weighted
    type(c_ptr), value :: result_ptr
    integer(c_int), value :: result_size
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_error) :: err
    real(wp), pointer :: mass_f(:)
    real(wp), pointer :: ref_coord_f(:,:)
    real(wp), pointer :: result_f(:)
    integer, pointer :: idx_f(:)
    integer, allocatable :: idx_copy(:)
    logical :: use_mass
    integer :: nstru

    ! Initialize
    call error_init(err)
    status = 0
    nstru_out = 0

    ! Validate inputs
    if (n_analysis <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_c: n_analysis must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(mass_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_c: mass_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(ref_coord_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_c: ref_coord_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(result_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_c: result_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (result_size <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_c: result_size must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create zero-copy views of arrays from Python
    call C_F_POINTER(mass_ptr, mass_f, [n_atoms])
    call C_F_POINTER(ref_coord_ptr, ref_coord_f, [3, n_atoms])
    call C_F_POINTER(result_ptr, result_f, [result_size])

    ! Convert analysis indices from C pointer to Fortran array
    call C_F_POINTER(analysis_idx, idx_f, [n_analysis])
    allocate(idx_copy(n_analysis))
    idx_copy(:) = idx_f(:)

    ! Convert mass_weighted to logical
    use_mass = (mass_weighted /= 0)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Run analysis
    write(MsgOut,'(A)') '[STEP1] RMSD Analysis (no fitting)'
    write(MsgOut,'(A)') ' '

    call analyze(mass_f, ref_coord_f, s_trajes_c, ana_period, &
                 idx_copy, n_analysis, use_mass, result_f, nstru)

    nstru_out = nstru

    ! Cleanup local arrays (views don't need deallocation)
    deallocate(idx_copy)

  end subroutine rmsd_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rmsd_analysis_fitting_c
  !> @brief        RMSD analysis with fitting (zerocopy, pre-allocated)
  !! @authors      Claude Code
  !! @param[in]    mass_ptr       : pointer to mass array (from Python NumPy)
  !! @param[in]    ref_coord_ptr  : pointer to reference coordinates (3, n_atoms)
  !! @param[in]    n_atoms        : number of atoms
  !! @param[in]    s_trajes_c     : trajectories C structure
  !! @param[in]    ana_period     : analysis period
  !! @param[in]    fitting_idx_ptr: pointer to fitting atom indices (1-indexed)
  !! @param[in]    n_fitting      : number of fitting atoms
  !! @param[in]    analysis_idx_ptr: pointer to analysis atom indices (1-indexed)
  !! @param[in]    n_analysis     : number of analysis atoms
  !! @param[in]    fitting_method : fitting method (1-6)
  !! @param[in]    mass_weighted  : use mass weighting (0 or 1)
  !! @param[in]    result_ptr     : pointer to pre-allocated result array
  !! @param[in]    result_size    : size of pre-allocated result array
  !! @param[out]   nstru_out      : number of frames analyzed
  !! @param[out]   status         : error status
  !! @param[out]   msg            : error message
  !! @param[in]    msglen         : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rmsd_analysis_fitting_c(mass_ptr, ref_coord_ptr, n_atoms, &
                                     s_trajes_c, ana_period, &
                                     fitting_idx_ptr, n_fitting, &
                                     analysis_idx_ptr, n_analysis, &
                                     fitting_method, mass_weighted, &
                                     result_ptr, result_size, nstru_out, &
                                     status, msg, msglen) &
        bind(C, name="rmsd_analysis_fitting_c")
    implicit none

    ! Arguments
    type(c_ptr), value :: mass_ptr
    type(c_ptr), value :: ref_coord_ptr
    integer(c_int), value :: n_atoms
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer(c_int), value :: ana_period
    type(c_ptr), value :: fitting_idx_ptr
    integer(c_int), value :: n_fitting
    type(c_ptr), value :: analysis_idx_ptr
    integer(c_int), value :: n_analysis
    integer(c_int), value :: fitting_method
    integer(c_int), value :: mass_weighted
    type(c_ptr), value :: result_ptr
    integer(c_int), value :: result_size
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_error) :: err
    real(wp), pointer :: mass_f(:)
    real(wp), pointer :: ref_coord_f(:,:)
    real(wp), pointer :: result_f(:)
    integer, pointer :: fitting_idx_f(:)
    integer, pointer :: analysis_idx_f(:)
    integer, allocatable :: fitting_idx_copy(:)
    integer, allocatable :: analysis_idx_copy(:)
    logical :: use_mass
    integer :: nstru

    ! Initialize
    call error_init(err)
    status = 0
    nstru_out = 0

    ! Validate inputs
    if (n_fitting <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: n_fitting must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (n_analysis <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: n_analysis must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(mass_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: mass_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(ref_coord_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: ref_coord_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(fitting_idx_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: fitting_idx_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(analysis_idx_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: analysis_idx_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(result_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: result_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (result_size <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_fitting_c: result_size must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create zero-copy views of arrays from Python
    call C_F_POINTER(mass_ptr, mass_f, [n_atoms])
    call C_F_POINTER(ref_coord_ptr, ref_coord_f, [3, n_atoms])
    call C_F_POINTER(result_ptr, result_f, [result_size])

    ! Convert indices from C pointer to Fortran array
    call C_F_POINTER(fitting_idx_ptr, fitting_idx_f, [n_fitting])
    call C_F_POINTER(analysis_idx_ptr, analysis_idx_f, [n_analysis])

    allocate(fitting_idx_copy(n_fitting))
    allocate(analysis_idx_copy(n_analysis))
    fitting_idx_copy(:) = fitting_idx_f(:)
    analysis_idx_copy(:) = analysis_idx_f(:)

    ! Convert mass_weighted to logical
    use_mass = (mass_weighted /= 0)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Run analysis
    write(MsgOut,'(A)') '[STEP1] RMSD Analysis (with fitting)'
    write(MsgOut,'(A)') ' '

    call analyze_with_fitting(mass_f, ref_coord_f, n_atoms, &
                              s_trajes_c, ana_period, &
                              fitting_idx_copy, n_fitting, &
                              analysis_idx_copy, n_analysis, &
                              fitting_method, use_mass, &
                              result_f, nstru)

    nstru_out = nstru

    ! Cleanup local arrays (views don't need deallocation)
    deallocate(fitting_idx_copy)
    deallocate(analysis_idx_copy)

  end subroutine rmsd_analysis_fitting_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rmsd_analysis_lazy_c
  !> @brief        RMSD analysis with lazy DCD loading (memory efficient)
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rmsd_analysis_lazy_c(dcd_filename, filename_len, trj_type, &
                                  mass_ptr, ref_coord_ptr, n_atoms, &
                                  ana_period, &
                                  fitting_idx_ptr, n_fitting, &
                                  analysis_idx_ptr, n_analysis, &
                                  fitting_method, mass_weighted, &
                                  result_ptr, result_size, &
                                  nstru_out, dcd_nframe_out, dcd_natom_out, &
                                  status, msg, msglen) &
        bind(C, name="rmsd_analysis_lazy_c")
    implicit none

    ! Arguments
    character(kind=c_char), intent(in) :: dcd_filename(*)
    integer(c_int), value :: filename_len
    integer(c_int), value :: trj_type
    type(c_ptr), value :: mass_ptr
    type(c_ptr), value :: ref_coord_ptr
    integer(c_int), value :: n_atoms
    integer(c_int), value :: ana_period
    type(c_ptr), value :: fitting_idx_ptr
    integer(c_int), value :: n_fitting
    type(c_ptr), value :: analysis_idx_ptr
    integer(c_int), value :: n_analysis
    integer(c_int), value :: fitting_method
    integer(c_int), value :: mass_weighted
    type(c_ptr), value :: result_ptr
    integer(c_int), value :: result_size
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: dcd_nframe_out
    integer(c_int), intent(out) :: dcd_natom_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_error) :: err
    type(s_trj_source) :: source
    type(s_trajectory) :: trajectory
    type(s_fitting) :: fitting
    character(MaxFilename) :: filename_f
    real(wp), pointer :: mass_f(:)
    real(wp), pointer :: ref_coord_f(:,:)
    real(wp), pointer :: result_f(:)
    integer, pointer :: fitting_idx_f(:)
    integer, pointer :: analysis_idx_f(:)
    integer, allocatable :: fitting_idx_copy(:)
    integer, allocatable :: analysis_idx_copy(:)
    real(wp), allocatable :: mass_fitting(:)
    real(wp), allocatable :: coord_work(:,:)
    logical :: use_mass, use_fitting
    integer :: nstru, frame_status, idx, iatom
    ! Fitting-related variables
    real(wp) :: rot_matrix(3,3)
    real(wp) :: com_ref(3), com_mov(3)
    real(wp) :: fit_rmsd
    integer :: ierr
    real(wp) :: rmsd, tot_mass, weight
    integer :: i

    ! Initialize
    call error_init(err)
    status = 0
    nstru_out = 0
    dcd_nframe_out = 0
    dcd_natom_out = 0

    ! Convert C string to Fortran string
    filename_f = ''
    do i = 1, min(filename_len, MaxFilename)
      if (dcd_filename(i) == c_null_char) exit
      filename_f(i:i) = dcd_filename(i)
    end do

    ! Validate inputs
    if (n_analysis <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_lazy_c: n_analysis must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(mass_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_lazy_c: mass_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(ref_coord_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_lazy_c: ref_coord_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(result_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_lazy_c: result_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (result_size <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rmsd_analysis_lazy_c: result_size must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create zero-copy views of arrays from Python
    call C_F_POINTER(mass_ptr, mass_f, [n_atoms])
    call C_F_POINTER(ref_coord_ptr, ref_coord_f, [3, n_atoms])
    call C_F_POINTER(result_ptr, result_f, [result_size])

    ! Convert analysis indices
    call C_F_POINTER(analysis_idx_ptr, analysis_idx_f, [n_analysis])
    allocate(analysis_idx_copy(n_analysis))
    analysis_idx_copy(:) = analysis_idx_f(:)

    ! Check for fitting
    use_fitting = (n_fitting > 0 .and. fitting_method > 0 .and. &
                   c_associated(fitting_idx_ptr))

    if (use_fitting) then
      call C_F_POINTER(fitting_idx_ptr, fitting_idx_f, [n_fitting])
      allocate(fitting_idx_copy(n_fitting))
      fitting_idx_copy(:) = fitting_idx_f(:)
    end if

    ! Convert mass_weighted to logical
    use_mass = (mass_weighted /= 0)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Initialize lazy DCD source
    write(MsgOut,'(A)') '[STEP1] Initialize Lazy DCD Source'
    write(MsgOut,'(A)') ' '

    call init_source_lazy_dcd(source, trim(filename_f), trj_type, ana_period)

    ! Return DCD info
    dcd_nframe_out = source%dcd_nframe
    dcd_natom_out = source%dcd_natom

    ! Check atom count
    if (source%dcd_natom /= n_atoms) then
      call error_set(err, ERROR_ATOM_COUNT, &
                     "rmsd_analysis_lazy_c: atom count mismatch")
      call error_to_c(err, status, msg, msglen)
      call finalize_source(source)
      deallocate(analysis_idx_copy)
      if (allocated(fitting_idx_copy)) deallocate(fitting_idx_copy)
      return
    end if

    ! Allocate work arrays
    allocate(coord_work(3, n_atoms))
    allocate(mass_fitting(n_atoms))

    ! Prepare mass array for fitting
    if (use_fitting .and. use_mass) then
      mass_fitting(:) = mass_f(:)
    else
      mass_fitting(:) = 1.0_wp
    end if

    ! Main analysis loop
    write(MsgOut,'(A)') '[STEP2] RMSD Analysis (lazy loading)'
    write(MsgOut,'(A)') ' '

    nstru = 0

    do while (has_more_frames(source))

      ! Get next frame via lazy loading
      call get_next_frame(source, trajectory, frame_status)
      if (frame_status /= 0) exit

      nstru = nstru + 1
      write(MsgOut,*) '      number of structures = ', nstru

      ! Copy coordinates for fitting
      coord_work(:,:) = trajectory%coord(:,:)

      ! Apply fitting if requested
      if (use_fitting) then
        select case(fitting_method)
        case(FittingMethodNO)
          ! no fitting

        case(FittingMethodTR_ROT)
          call fit_trrot(n_fitting, fitting_idx_copy, ref_coord_f, mass_fitting, &
                         coord_work, rot_matrix, com_ref, com_mov, fit_rmsd, ierr)
          call transform(n_atoms, rot_matrix, com_ref, com_mov, coord_work)

        case(FittingMethodTR)
          call fit_trans(n_fitting, fitting_idx_copy, ref_coord_f, coord_work, mass_fitting, &
                         .true., rot_matrix, com_ref, com_mov, fit_rmsd, ierr)
          call transform(n_atoms, rot_matrix, com_ref, com_mov, coord_work)

        case(FittingMethodXYTR)
          call fit_trans(n_fitting, fitting_idx_copy, ref_coord_f, coord_work, mass_fitting, &
                         .false., rot_matrix, com_ref, com_mov, fit_rmsd, ierr)
          call transform(n_atoms, rot_matrix, com_ref, com_mov, coord_work)

        case default
          ! For other methods (TR_ZROT, XYTR_ZROT), fall back to TR+ROT
          call fit_trrot(n_fitting, fitting_idx_copy, ref_coord_f, mass_fitting, &
                         coord_work, rot_matrix, com_ref, com_mov, fit_rmsd, ierr)
          call transform(n_atoms, rot_matrix, com_ref, com_mov, coord_work)

        end select
      end if

      ! Compute RMSD for analysis atoms
      rmsd = 0.0_wp
      tot_mass = 0.0_wp

      do iatom = 1, n_analysis
        idx = analysis_idx_copy(iatom)

        if (use_mass) then
          weight = mass_f(idx)
        else
          weight = 1.0_wp
        end if

        rmsd = rmsd + weight * ( &
               (ref_coord_f(1,idx) - coord_work(1,idx))**2 + &
               (ref_coord_f(2,idx) - coord_work(2,idx))**2 + &
               (ref_coord_f(3,idx) - coord_work(3,idx))**2)

        tot_mass = tot_mass + weight
      end do

      if (tot_mass > EPS) then
        rmsd = sqrt(rmsd / tot_mass)
      else
        rmsd = 0.0_wp
      end if

      ! Write result to zerocopy array
      if (nstru <= result_size) then
        result_f(nstru) = rmsd
      end if

      write(MsgOut,'(a,f10.5)') '              RMSD of analysis atoms = ', rmsd
      write(MsgOut,*) ''

    end do

    nstru_out = nstru

    ! Cleanup
    call finalize_source(source)
    deallocate(coord_work)
    deallocate(mass_fitting)
    deallocate(analysis_idx_copy)
    if (allocated(fitting_idx_copy)) deallocate(fitting_idx_copy)
    if (allocated(trajectory%coord)) deallocate(trajectory%coord)

  end subroutine rmsd_analysis_lazy_c

end module rmsd_c_mod
