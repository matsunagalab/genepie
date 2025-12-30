!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> C language interface for trj_analysis (zerocopy)
!! @brief   C language interface of trj_analysis
!! @authors Claude Code
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module trj_c_mod
  use, intrinsic :: iso_c_binding
  use s_trajectories_c_mod
  use trj_impl_mod

  use trajectory_str_mod
  use molecules_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: trj_analysis_c
  public :: trj_analysis_com_c

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_c
  !> @brief        Trajectory analysis (zerocopy, pre-allocated results)
  !! @authors      Claude Code
  !! @param[in]    s_trajes_c   : trajectories C structure
  !! @param[in]    ana_period   : analysis period
  !! @param[in]    dist_list_ptr: pointer to distance atom pairs (2, n_dist)
  !! @param[in]    n_dist       : number of distance measurements
  !! @param[in]    angl_list_ptr: pointer to angle atom triplets (3, n_angl)
  !! @param[in]    n_angl       : number of angle measurements
  !! @param[in]    tors_list_ptr: pointer to torsion atom quadruplets (4, n_tors)
  !! @param[in]    n_tors       : number of torsion measurements
  !! @param[in]    dist_ptr     : pre-allocated distance results (n_dist, n_frames)
  !! @param[in]    dist_size    : total size of dist array (n_dist * n_frames)
  !! @param[in]    angl_ptr     : pre-allocated angle results (n_angl, n_frames)
  !! @param[in]    angl_size    : total size of angl array
  !! @param[in]    tors_ptr     : pre-allocated torsion results (n_tors, n_frames)
  !! @param[in]    tors_size    : total size of tors array
  !! @param[out]   nstru_out    : actual number of structures analyzed
  !! @param[out]   status       : error status
  !! @param[out]   msg          : error message
  !! @param[in]    msglen       : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trj_analysis_c(s_trajes_c, ana_period, &
                            dist_list_ptr, n_dist, &
                            angl_list_ptr, n_angl, &
                            tors_list_ptr, n_tors, &
                            dist_ptr, dist_size, &
                            angl_ptr, angl_size, &
                            tors_ptr, tors_size, &
                            nstru_out, status, msg, msglen) &
        bind(C, name="trj_analysis_c")
    implicit none

    ! Arguments
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer(c_int), value :: ana_period
    type(c_ptr), value :: dist_list_ptr
    integer(c_int), value :: n_dist
    type(c_ptr), value :: angl_list_ptr
    integer(c_int), value :: n_angl
    type(c_ptr), value :: tors_list_ptr
    integer(c_int), value :: n_tors
    type(c_ptr), value :: dist_ptr
    integer(c_int), value :: dist_size
    type(c_ptr), value :: angl_ptr
    integer(c_int), value :: angl_size
    type(c_ptr), value :: tors_ptr
    integer(c_int), value :: tors_size
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    integer, pointer :: dist_list_f(:,:)
    integer, pointer :: angl_list_f(:,:)
    integer, pointer :: tors_list_f(:,:)
    real(wp), pointer :: dist_f(:,:)
    real(wp), pointer :: angl_f(:,:)
    real(wp), pointer :: tors_f(:,:)
    integer, allocatable :: dist_list_copy(:,:)
    integer, allocatable :: angl_list_copy(:,:)
    integer, allocatable :: tors_list_copy(:,:)
    integer :: n_frames, nstru_local

    ! Initialize
    status = 0
    nstru_out = 0
    n_frames = s_trajes_c%nframe / ana_period

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    write(MsgOut,'(A)') '[STEP1] Trajectory Analysis'
    write(MsgOut,'(A)') ' '

    ! Convert C pointers to Fortran arrays (input lists)
    if (n_dist > 0 .and. c_associated(dist_list_ptr)) then
      call C_F_POINTER(dist_list_ptr, dist_list_f, [2, n_dist])
      allocate(dist_list_copy(2, n_dist))
      dist_list_copy = dist_list_f
    else
      allocate(dist_list_copy(2, 0))
    end if

    if (n_angl > 0 .and. c_associated(angl_list_ptr)) then
      call C_F_POINTER(angl_list_ptr, angl_list_f, [3, n_angl])
      allocate(angl_list_copy(3, n_angl))
      angl_list_copy = angl_list_f
    else
      allocate(angl_list_copy(3, 0))
    end if

    if (n_tors > 0 .and. c_associated(tors_list_ptr)) then
      call C_F_POINTER(tors_list_ptr, tors_list_f, [4, n_tors])
      allocate(tors_list_copy(4, n_tors))
      tors_list_copy = tors_list_f
    else
      allocate(tors_list_copy(4, 0))
    end if

    ! Create views of pre-allocated result arrays
    if (n_dist > 0 .and. c_associated(dist_ptr) .and. dist_size > 0) then
      call C_F_POINTER(dist_ptr, dist_f, [n_dist, n_frames])
    else
      nullify(dist_f)
    end if

    if (n_angl > 0 .and. c_associated(angl_ptr) .and. angl_size > 0) then
      call C_F_POINTER(angl_ptr, angl_f, [n_angl, n_frames])
    else
      nullify(angl_f)
    end if

    if (n_tors > 0 .and. c_associated(tors_ptr) .and. tors_size > 0) then
      call C_F_POINTER(tors_ptr, tors_f, [n_tors, n_frames])
    else
      nullify(tors_f)
    end if

    ! Run analysis with pre-allocated arrays
    call analyze(s_trajes_c, ana_period, &
                 dist_list_copy, n_dist, &
                 angl_list_copy, n_angl, &
                 tors_list_copy, n_tors, &
                 dist_f, angl_f, tors_f, nstru_local)

    nstru_out = nstru_local

    ! Cleanup local arrays
    deallocate(dist_list_copy)
    deallocate(angl_list_copy)
    deallocate(tors_list_copy)

  end subroutine trj_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_com_c
  !> @brief        Trajectory analysis with COM (zerocopy, pre-allocated)
  !! @authors      Claude Code
  !! @note         All result arrays are pre-allocated by Python
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trj_analysis_com_c(mass_ptr, n_atoms, &
                                s_trajes_c, ana_period, &
                                dist_list_ptr, n_dist, &
                                angl_list_ptr, n_angl, &
                                tors_list_ptr, n_tors, &
                                cdis_atoms_ptr, n_cdis_atoms, &
                                cdis_offsets_ptr, n_cdis_offsets, &
                                cdis_pairs_ptr, n_cdis, &
                                cang_atoms_ptr, n_cang_atoms, &
                                cang_offsets_ptr, n_cang_offsets, &
                                cang_triplets_ptr, n_cang, &
                                ctor_atoms_ptr, n_ctor_atoms, &
                                ctor_offsets_ptr, n_ctor_offsets, &
                                ctor_quads_ptr, n_ctor, &
                                dist_ptr, dist_size, &
                                angl_ptr, angl_size, &
                                tors_ptr, tors_size, &
                                cdis_result_ptr, cdis_size, &
                                cang_result_ptr, cang_size, &
                                ctor_result_ptr, ctor_size, &
                                nstru_out, status, msg, msglen) &
        bind(C, name="trj_analysis_com_c")
    implicit none

    ! Arguments - mass array
    type(c_ptr), value :: mass_ptr
    integer(c_int), value :: n_atoms

    ! Arguments - trajectories
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer(c_int), value :: ana_period

    ! Arguments - atom-based measurements
    type(c_ptr), value :: dist_list_ptr
    integer(c_int), value :: n_dist
    type(c_ptr), value :: angl_list_ptr
    integer(c_int), value :: n_angl
    type(c_ptr), value :: tors_list_ptr
    integer(c_int), value :: n_tors

    ! Arguments - COM distance
    type(c_ptr), value :: cdis_atoms_ptr
    integer(c_int), value :: n_cdis_atoms
    type(c_ptr), value :: cdis_offsets_ptr
    integer(c_int), value :: n_cdis_offsets
    type(c_ptr), value :: cdis_pairs_ptr
    integer(c_int), value :: n_cdis

    ! Arguments - COM angle
    type(c_ptr), value :: cang_atoms_ptr
    integer(c_int), value :: n_cang_atoms
    type(c_ptr), value :: cang_offsets_ptr
    integer(c_int), value :: n_cang_offsets
    type(c_ptr), value :: cang_triplets_ptr
    integer(c_int), value :: n_cang

    ! Arguments - COM torsion
    type(c_ptr), value :: ctor_atoms_ptr
    integer(c_int), value :: n_ctor_atoms
    type(c_ptr), value :: ctor_offsets_ptr
    integer(c_int), value :: n_ctor_offsets
    type(c_ptr), value :: ctor_quads_ptr
    integer(c_int), value :: n_ctor

    ! Pre-allocated output arrays
    type(c_ptr), value :: dist_ptr
    integer(c_int), value :: dist_size
    type(c_ptr), value :: angl_ptr
    integer(c_int), value :: angl_size
    type(c_ptr), value :: tors_ptr
    integer(c_int), value :: tors_size
    type(c_ptr), value :: cdis_result_ptr
    integer(c_int), value :: cdis_size
    type(c_ptr), value :: cang_result_ptr
    integer(c_int), value :: cang_size
    type(c_ptr), value :: ctor_result_ptr
    integer(c_int), value :: ctor_size

    ! Output
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables - Fortran pointers from C
    real(wp), pointer :: mass_f(:)
    integer, pointer :: dist_list_f(:,:)
    integer, pointer :: angl_list_f(:,:)
    integer, pointer :: tors_list_f(:,:)
    integer, pointer :: cdis_atoms_f(:)
    integer, pointer :: cdis_offsets_f(:)
    integer, pointer :: cdis_pairs_f(:)
    integer, pointer :: cang_atoms_f(:)
    integer, pointer :: cang_offsets_f(:)
    integer, pointer :: cang_triplets_f(:)
    integer, pointer :: ctor_atoms_f(:)
    integer, pointer :: ctor_offsets_f(:)
    integer, pointer :: ctor_quads_f(:)

    ! Pre-allocated result arrays (views)
    real(wp), pointer :: distance_f(:,:)
    real(wp), pointer :: angle_f(:,:)
    real(wp), pointer :: torsion_f(:,:)
    real(wp), pointer :: cdis_f(:,:)
    real(wp), pointer :: cang_f(:,:)
    real(wp), pointer :: ctor_f(:,:)

    ! Local copies for 1-indexed conversion
    integer, allocatable :: dist_list_copy(:,:)
    integer, allocatable :: angl_list_copy(:,:)
    integer, allocatable :: tors_list_copy(:,:)
    integer, allocatable :: cdis_atoms_copy(:)
    integer, allocatable :: cdis_offsets_copy(:)
    integer, allocatable :: cdis_pairs_copy(:)
    integer, allocatable :: cang_atoms_copy(:)
    integer, allocatable :: cang_offsets_copy(:)
    integer, allocatable :: cang_triplets_copy(:)
    integer, allocatable :: ctor_atoms_copy(:)
    integer, allocatable :: ctor_offsets_copy(:)
    integer, allocatable :: ctor_quads_copy(:)

    integer :: nstru_local
    integer :: num_out_frame
    integer :: n_cdis_groups, n_cang_groups, n_ctor_groups

    ! Initialize
    status = 0
    nstru_out = 0
    if (msglen > 0) msg(1) = c_null_char

    num_out_frame = s_trajes_c%nframe / ana_period

    ! Calculate number of groups
    n_cdis_groups = n_cdis_offsets - 1
    n_cang_groups = n_cang_offsets - 1
    n_ctor_groups = n_ctor_offsets - 1

    ! Get mass array
    call C_F_POINTER(mass_ptr, mass_f, [n_atoms])

    ! Create views of pre-allocated result arrays
    if (n_dist > 0 .and. c_associated(dist_ptr)) then
      call C_F_POINTER(dist_ptr, distance_f, [n_dist, num_out_frame])
    else
      allocate(distance_f(1,1))  ! Dummy
    end if

    if (n_angl > 0 .and. c_associated(angl_ptr)) then
      call C_F_POINTER(angl_ptr, angle_f, [n_angl, num_out_frame])
    else
      allocate(angle_f(1,1))  ! Dummy
    end if

    if (n_tors > 0 .and. c_associated(tors_ptr)) then
      call C_F_POINTER(tors_ptr, torsion_f, [n_tors, num_out_frame])
    else
      allocate(torsion_f(1,1))  ! Dummy
    end if

    if (n_cdis > 0 .and. c_associated(cdis_result_ptr)) then
      call C_F_POINTER(cdis_result_ptr, cdis_f, [n_cdis, num_out_frame])
    else
      allocate(cdis_f(1,1))  ! Dummy
    end if

    if (n_cang > 0 .and. c_associated(cang_result_ptr)) then
      call C_F_POINTER(cang_result_ptr, cang_f, [n_cang, num_out_frame])
    else
      allocate(cang_f(1,1))  ! Dummy
    end if

    if (n_ctor > 0 .and. c_associated(ctor_result_ptr)) then
      call C_F_POINTER(ctor_result_ptr, ctor_f, [n_ctor, num_out_frame])
    else
      allocate(ctor_f(1,1))  ! Dummy
    end if

    ! Convert atom-based lists (already 1-indexed from Python)
    if (n_dist > 0) then
      call C_F_POINTER(dist_list_ptr, dist_list_f, [2, n_dist])
      allocate(dist_list_copy(2, n_dist))
      dist_list_copy = dist_list_f
    else
      allocate(dist_list_copy(2, 1))
      dist_list_copy = 0
    end if

    if (n_angl > 0) then
      call C_F_POINTER(angl_list_ptr, angl_list_f, [3, n_angl])
      allocate(angl_list_copy(3, n_angl))
      angl_list_copy = angl_list_f
    else
      allocate(angl_list_copy(3, 1))
      angl_list_copy = 0
    end if

    if (n_tors > 0) then
      call C_F_POINTER(tors_list_ptr, tors_list_f, [4, n_tors])
      allocate(tors_list_copy(4, n_tors))
      tors_list_copy = tors_list_f
    else
      allocate(tors_list_copy(4, 1))
      tors_list_copy = 0
    end if

    ! Convert COM distance arrays (atoms already 1-indexed from Python)
    if (n_cdis_atoms > 0 .and. n_cdis > 0) then
      call C_F_POINTER(cdis_atoms_ptr, cdis_atoms_f, [n_cdis_atoms])
      call C_F_POINTER(cdis_offsets_ptr, cdis_offsets_f, [n_cdis_offsets])
      call C_F_POINTER(cdis_pairs_ptr, cdis_pairs_f, [2 * n_cdis])
      allocate(cdis_atoms_copy(n_cdis_atoms))
      allocate(cdis_offsets_copy(n_cdis_offsets))
      allocate(cdis_pairs_copy(2 * n_cdis))
      cdis_atoms_copy = cdis_atoms_f       ! Atoms already 1-indexed
      cdis_offsets_copy = cdis_offsets_f   ! Offsets stay 0-based
      cdis_pairs_copy = cdis_pairs_f       ! Group indices stay 0-based
    else
      allocate(cdis_atoms_copy(1))
      allocate(cdis_offsets_copy(2))
      allocate(cdis_pairs_copy(2))
      cdis_atoms_copy = 0
      cdis_offsets_copy = 0
      cdis_pairs_copy = 0
    end if

    ! Convert COM angle arrays (atoms already 1-indexed from Python)
    if (n_cang_atoms > 0 .and. n_cang > 0) then
      call C_F_POINTER(cang_atoms_ptr, cang_atoms_f, [n_cang_atoms])
      call C_F_POINTER(cang_offsets_ptr, cang_offsets_f, [n_cang_offsets])
      call C_F_POINTER(cang_triplets_ptr, cang_triplets_f, [3 * n_cang])
      allocate(cang_atoms_copy(n_cang_atoms))
      allocate(cang_offsets_copy(n_cang_offsets))
      allocate(cang_triplets_copy(3 * n_cang))
      cang_atoms_copy = cang_atoms_f       ! Atoms already 1-indexed
      cang_offsets_copy = cang_offsets_f
      cang_triplets_copy = cang_triplets_f
    else
      allocate(cang_atoms_copy(1))
      allocate(cang_offsets_copy(2))
      allocate(cang_triplets_copy(3))
      cang_atoms_copy = 0
      cang_offsets_copy = 0
      cang_triplets_copy = 0
    end if

    ! Convert COM torsion arrays (atoms already 1-indexed from Python)
    if (n_ctor_atoms > 0 .and. n_ctor > 0) then
      call C_F_POINTER(ctor_atoms_ptr, ctor_atoms_f, [n_ctor_atoms])
      call C_F_POINTER(ctor_offsets_ptr, ctor_offsets_f, [n_ctor_offsets])
      call C_F_POINTER(ctor_quads_ptr, ctor_quads_f, [4 * n_ctor])
      allocate(ctor_atoms_copy(n_ctor_atoms))
      allocate(ctor_offsets_copy(n_ctor_offsets))
      allocate(ctor_quads_copy(4 * n_ctor))
      ctor_atoms_copy = ctor_atoms_f       ! Atoms already 1-indexed
      ctor_offsets_copy = ctor_offsets_f
      ctor_quads_copy = ctor_quads_f
    else
      allocate(ctor_atoms_copy(1))
      allocate(ctor_offsets_copy(2))
      allocate(ctor_quads_copy(4))
      ctor_atoms_copy = 0
      ctor_offsets_copy = 0
      ctor_quads_copy = 0
    end if

    ! Run analysis
    write(MsgOut,'(A)') '[STEP1] Trajectory Analysis with COM'
    write(MsgOut,'(A)') ' '

    call analyze_com(mass_f, s_trajes_c, ana_period, &
                     dist_list_copy, n_dist, &
                     angl_list_copy, n_angl, &
                     tors_list_copy, n_tors, &
                     cdis_atoms_copy, cdis_offsets_copy, cdis_pairs_copy, &
                     n_cdis, n_cdis_groups, &
                     cang_atoms_copy, cang_offsets_copy, cang_triplets_copy, &
                     n_cang, n_cang_groups, &
                     ctor_atoms_copy, ctor_offsets_copy, ctor_quads_copy, &
                     n_ctor, n_ctor_groups, &
                     distance_f, angle_f, torsion_f, &
                     cdis_f, cang_f, ctor_f, nstru_local)

    nstru_out = nstru_local

    ! Clean up dummy allocations and copies
    if (.not. (n_dist > 0 .and. c_associated(dist_ptr))) deallocate(distance_f)
    if (.not. (n_angl > 0 .and. c_associated(angl_ptr))) deallocate(angle_f)
    if (.not. (n_tors > 0 .and. c_associated(tors_ptr))) deallocate(torsion_f)
    if (.not. (n_cdis > 0 .and. c_associated(cdis_result_ptr))) deallocate(cdis_f)
    if (.not. (n_cang > 0 .and. c_associated(cang_result_ptr))) deallocate(cang_f)
    if (.not. (n_ctor > 0 .and. c_associated(ctor_result_ptr))) deallocate(ctor_f)

    deallocate(dist_list_copy)
    deallocate(angl_list_copy)
    deallocate(tors_list_copy)
    deallocate(cdis_atoms_copy)
    deallocate(cdis_offsets_copy)
    deallocate(cdis_pairs_copy)
    deallocate(cang_atoms_copy)
    deallocate(cang_offsets_copy)
    deallocate(cang_triplets_copy)
    deallocate(ctor_atoms_copy)
    deallocate(ctor_offsets_copy)
    deallocate(ctor_quads_copy)

  end subroutine trj_analysis_com_c

end module trj_c_mod
