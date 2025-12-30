!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> C language interface for trj_analysis
!! @brief   C language interface of trj_analysis
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module trj_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use trj_impl_mod

  use ta_control_mod
  use ta_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use error_mod
  implicit none

  public :: trj_analysis_c
  public :: trj_analysis_zerocopy_c
  public :: trj_analysis_zerocopy_com_c
  public :: deallocate_trj_results_c

  ! Module-level pointers for results (to be deallocated later)
  real(wp), pointer, save :: distance_ptr(:,:) => null()
  real(wp), pointer, save :: angle_ptr(:,:) => null()
  real(wp), pointer, save :: torsion_ptr(:,:) => null()
  real(wp), pointer, save :: cdis_ptr(:,:) => null()
  real(wp), pointer, save :: cang_ptr(:,:) => null()
  real(wp), pointer, save :: ctor_ptr(:,:) => null()

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_c
  !> @brief        C language interface for trj_analysis
  !! @param[in] molecule : structure of molecule information
  !! @param[in] s_trajes_c : structure of trajectories information
  !! @param[in] ana_period : interval for sampling frames
  !! @param[in] ctrl_path : path to the control file
  !! @param[out] result_distance : array storing calculated distances
  !! @param[out] num_distance : number of distance results stored in result_distance
  !! @param[out] result_angle : array storing calculated angles
  !! @param[out] num_angle : number of angle results stored in result_angle
  !! @param[out] result_torsion : array storing calculated torsion angles
  !! @param[out] num_torsion : number of torsion angle results stored in result_torsion
  !! @param[out] result_cdis : array storing calculated contact distances
  !! @param[out] num_cdis : number of contact distance results stored in result_cdis
  !! @param[out] result_cang : array storing calculated contact angles
  !! @param[out] num_cang : number of contact angle results stored in result_cang
  !! @param[out] result_ctor : array storing calculated contact torsion angles
  !! @param[out] num_ctor : number of contact torsion angle results stored in result_ctor
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine trj_analysis_c(molecule, s_trajes_c, ana_period, &
                            ctrl_text, ctrl_len, &
                            result_distance, num_distance, &
                            result_angle, num_angle, &
                            result_torsion, num_torsion, &
                            result_cdis, num_cdis, &
                            result_cang, num_cang, &
                            result_ctor, num_ctor) &
        bind(C, name="trj_analysis_c")
    use conv_f_c_util
    implicit none
    type(s_molecule_c), intent(in) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer, intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out) :: result_distance
    integer(c_int), intent(out) :: num_distance
    type(c_ptr), intent(out) :: result_angle
    integer(c_int), intent(out) :: num_angle
    type(c_ptr), intent(out) :: result_torsion
    integer(c_int), intent(out) :: num_torsion
    type(c_ptr), intent(out) :: result_cdis
    integer(c_int), intent(out) :: num_cdis
    type(c_ptr), intent(out) :: result_cang
    integer(c_int), intent(out) :: num_cang
    type(c_ptr), intent(out) :: result_ctor
    integer(c_int), intent(out) :: num_ctor

    type(s_molecule) :: f_molecule

    ! Deallocate previous results if any
    call deallocate_trj_results_internal()

    call c2f_s_molecule(molecule, f_molecule)
    call trj_analysis_main( &
        f_molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, &
        distance_ptr, angle_ptr, torsion_ptr, cdis_ptr, cang_ptr, ctor_ptr)
    call dealloc_molecules_all(f_molecule)

    if (associated(distance_ptr)) then
      result_distance = c_loc(distance_ptr)
      num_distance = size(distance_ptr, dim=1)
    else
      result_distance = c_null_ptr
      num_distance = 0
    end if
    if (associated(angle_ptr)) then
      result_angle = c_loc(angle_ptr)
      num_angle = size(angle_ptr, dim=1)
    else
      result_angle = c_null_ptr
      num_angle = 0
    end if
    if (associated(torsion_ptr)) then
      result_torsion = c_loc(torsion_ptr)
      num_torsion = size(torsion_ptr, dim=1)
    else
      result_torsion = c_null_ptr
      num_torsion = 0
    end if
    if (associated(cdis_ptr)) then
      result_cdis = c_loc(cdis_ptr)
      num_cdis = size(cdis_ptr, dim=1)
    else
      result_cdis = c_null_ptr
      num_cdis = 0
    end if
    if (associated(cang_ptr)) then
      result_cang = c_loc(cang_ptr)
      num_cang = size(cang_ptr, dim=1)
    else
      result_cang = c_null_ptr
      num_cang = 0
    end if
    if (associated(ctor_ptr)) then
      result_ctor = c_loc(ctor_ptr)
      num_ctor = size(ctor_ptr, dim=1)
    else
      result_ctor = c_null_ptr
      num_ctor = 0
    end if
  end subroutine trj_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_zerocopy_c
  !> @brief        Trajectory analysis with zerocopy interface
  !! @authors      Claude Code
  !! @param[in]    s_trajes_c   : trajectories C structure
  !! @param[in]    ana_period   : analysis period
  !! @param[in]    dist_list_ptr: pointer to distance atom pairs (2, n_dist)
  !! @param[in]    n_dist       : number of distance measurements
  !! @param[in]    angl_list_ptr: pointer to angle atom triplets (3, n_angl)
  !! @param[in]    n_angl       : number of angle measurements
  !! @param[in]    tors_list_ptr: pointer to torsion atom quadruplets (4, n_tors)
  !! @param[in]    n_tors       : number of torsion measurements
  !! @param[out]   result_distance: pointer to distance results
  !! @param[out]   result_angle   : pointer to angle results
  !! @param[out]   result_torsion : pointer to torsion results
  !! @param[out]   n_frames       : number of output frames
  !! @param[out]   status         : error status
  !! @param[out]   msg            : error message
  !! @param[in]    msglen         : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trj_analysis_zerocopy_c(s_trajes_c, ana_period, &
                                     dist_list_ptr, n_dist, &
                                     angl_list_ptr, n_angl, &
                                     tors_list_ptr, n_tors, &
                                     result_distance, result_angle, result_torsion, &
                                     n_frames, status, msg, msglen) &
        bind(C, name="trj_analysis_zerocopy_c")
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
    type(c_ptr), intent(out) :: result_distance
    type(c_ptr), intent(out) :: result_angle
    type(c_ptr), intent(out) :: result_torsion
    integer(c_int), intent(out) :: n_frames
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_error) :: err
    integer, pointer :: dist_list_f(:,:)
    integer, pointer :: angl_list_f(:,:)
    integer, pointer :: tors_list_f(:,:)
    integer, allocatable :: dist_list_copy(:,:)
    integer, allocatable :: angl_list_copy(:,:)
    integer, allocatable :: tors_list_copy(:,:)

    ! Initialize
    call error_init(err)
    status = 0
    result_distance = c_null_ptr
    result_angle = c_null_ptr
    result_torsion = c_null_ptr
    n_frames = s_trajes_c%nframe / ana_period

    ! Deallocate previous results
    call deallocate_trj_results_internal()

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    write(MsgOut,'(A)') '[STEP1] Trajectory Analysis (zerocopy interface)'
    write(MsgOut,'(A)') ' '

    ! Convert C pointers to Fortran arrays (create copies for proper indexing)
    ! Note: Fortran array is (2, n_dist), C array is (n_dist, 2) row-major
    ! so C_F_POINTER with (2, n_dist) reads it correctly

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

    ! Run analysis
    call analyze_zerocopy(s_trajes_c, ana_period, &
                          dist_list_copy, n_dist, &
                          angl_list_copy, n_angl, &
                          tors_list_copy, n_tors, &
                          distance_ptr, angle_ptr, torsion_ptr)

    ! Return pointers to results
    if (associated(distance_ptr)) then
      result_distance = c_loc(distance_ptr)
    end if
    if (associated(angle_ptr)) then
      result_angle = c_loc(angle_ptr)
    end if
    if (associated(torsion_ptr)) then
      result_torsion = c_loc(torsion_ptr)
    end if

    ! Cleanup local arrays
    deallocate(dist_list_copy)
    deallocate(angl_list_copy)
    deallocate(tors_list_copy)

  end subroutine trj_analysis_zerocopy_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    trj_analysis_zerocopy_com_c
  !> @brief        Trajectory analysis with COM calculations (zerocopy C interface)
  !! @authors      Claude Code
  !! @param[in]    mass_ptr     : pointer to atomic masses
  !! @param[in]    n_atoms      : number of atoms
  !! @param[in]    s_trajes_c   : trajectories C structure
  !! @param[in]    ana_period   : analysis period
  !! @param[in]    dist_list_ptr: pointer to distance atom pairs (2, n_dist)
  !! @param[in]    n_dist       : number of distance measurements
  !! @param[in]    angl_list_ptr: pointer to angle atom triplets (3, n_angl)
  !! @param[in]    n_angl       : number of angle measurements
  !! @param[in]    tors_list_ptr: pointer to torsion atom quadruplets (4, n_tors)
  !! @param[in]    n_tors       : number of torsion measurements
  !! @param[in]    cdis_atoms_ptr: flat array of atom indices for COM distance
  !! @param[in]    n_cdis_atoms : total number of atoms in cdis_atoms
  !! @param[in]    cdis_offsets_ptr: offsets for each COM distance group
  !! @param[in]    n_cdis_offsets: number of offsets (n_cdis_groups + 1)
  !! @param[in]    cdis_pairs_ptr: group pair indices for COM distances
  !! @param[in]    n_cdis       : number of COM distance measurements
  !! @param[in]    cang_atoms_ptr: flat array of atom indices for COM angles
  !! @param[in]    n_cang_atoms : total number of atoms in cang_atoms
  !! @param[in]    cang_offsets_ptr: offsets for each COM angle group
  !! @param[in]    n_cang_offsets: number of offsets
  !! @param[in]    cang_triplets_ptr: group triplet indices for COM angles
  !! @param[in]    n_cang       : number of COM angle measurements
  !! @param[in]    ctor_atoms_ptr: flat array of atom indices for COM torsions
  !! @param[in]    n_ctor_atoms : total number of atoms in ctor_atoms
  !! @param[in]    ctor_offsets_ptr: offsets for each COM torsion group
  !! @param[in]    n_ctor_offsets: number of offsets
  !! @param[in]    ctor_quads_ptr: group quad indices for COM torsions
  !! @param[in]    n_ctor       : number of COM torsion measurements
  !! @param[out]   result_distance: pointer to distance results
  !! @param[out]   result_angle   : pointer to angle results
  !! @param[out]   result_torsion : pointer to torsion results
  !! @param[out]   result_cdis    : pointer to COM distance results
  !! @param[out]   result_cang    : pointer to COM angle results
  !! @param[out]   result_ctor    : pointer to COM torsion results
  !! @param[out]   n_frames       : number of output frames
  !! @param[out]   status         : error status
  !! @param[out]   msg            : error message
  !! @param[in]    msglen         : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine trj_analysis_zerocopy_com_c(mass_ptr, n_atoms, &
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
                                         result_distance, result_angle, result_torsion, &
                                         result_cdis, result_cang, result_ctor, &
                                         n_frames, status, msg, msglen) &
        bind(C, name="trj_analysis_zerocopy_com_c")
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

    ! Arguments - output
    type(c_ptr), intent(out) :: result_distance
    type(c_ptr), intent(out) :: result_angle
    type(c_ptr), intent(out) :: result_torsion
    type(c_ptr), intent(out) :: result_cdis
    type(c_ptr), intent(out) :: result_cang
    type(c_ptr), intent(out) :: result_ctor
    integer(c_int), intent(out) :: n_frames
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

    ! Local copies
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

    integer :: n_cdis_groups, n_cang_groups, n_ctor_groups

    ! Initialize
    status = 0
    result_distance = c_null_ptr
    result_angle = c_null_ptr
    result_torsion = c_null_ptr
    result_cdis = c_null_ptr
    result_cang = c_null_ptr
    result_ctor = c_null_ptr
    n_frames = s_trajes_c%nframe / ana_period

    ! Deallocate previous results
    call deallocate_trj_results_internal()

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    write(MsgOut,'(A)') '[STEP1] Trajectory Analysis with COM (zerocopy interface)'
    write(MsgOut,'(A)') ' '

    ! Get mass array (zerocopy)
    if (.not. c_associated(mass_ptr)) then
      status = 301
      call copy_error_msg("Mass pointer is null", msg, msglen)
      return
    end if
    call C_F_POINTER(mass_ptr, mass_f, [n_atoms])

    ! Convert atom-based measurement arrays
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

    ! Convert COM distance arrays
    n_cdis_groups = n_cdis_offsets - 1
    if (n_cdis > 0 .and. n_cdis_atoms > 0) then
      call C_F_POINTER(cdis_atoms_ptr, cdis_atoms_f, [n_cdis_atoms])
      allocate(cdis_atoms_copy(n_cdis_atoms))
      cdis_atoms_copy = cdis_atoms_f

      call C_F_POINTER(cdis_offsets_ptr, cdis_offsets_f, [n_cdis_offsets])
      allocate(cdis_offsets_copy(n_cdis_offsets))
      cdis_offsets_copy = cdis_offsets_f

      call C_F_POINTER(cdis_pairs_ptr, cdis_pairs_f, [2 * n_cdis])
      allocate(cdis_pairs_copy(2 * n_cdis))
      cdis_pairs_copy = cdis_pairs_f
    else
      allocate(cdis_atoms_copy(0))
      allocate(cdis_offsets_copy(1))
      cdis_offsets_copy(1) = 0
      allocate(cdis_pairs_copy(0))
      n_cdis_groups = 0
    end if

    ! Convert COM angle arrays
    n_cang_groups = n_cang_offsets - 1
    if (n_cang > 0 .and. n_cang_atoms > 0) then
      call C_F_POINTER(cang_atoms_ptr, cang_atoms_f, [n_cang_atoms])
      allocate(cang_atoms_copy(n_cang_atoms))
      cang_atoms_copy = cang_atoms_f

      call C_F_POINTER(cang_offsets_ptr, cang_offsets_f, [n_cang_offsets])
      allocate(cang_offsets_copy(n_cang_offsets))
      cang_offsets_copy = cang_offsets_f

      call C_F_POINTER(cang_triplets_ptr, cang_triplets_f, [3 * n_cang])
      allocate(cang_triplets_copy(3 * n_cang))
      cang_triplets_copy = cang_triplets_f
    else
      allocate(cang_atoms_copy(0))
      allocate(cang_offsets_copy(1))
      cang_offsets_copy(1) = 0
      allocate(cang_triplets_copy(0))
      n_cang_groups = 0
    end if

    ! Convert COM torsion arrays
    n_ctor_groups = n_ctor_offsets - 1
    if (n_ctor > 0 .and. n_ctor_atoms > 0) then
      call C_F_POINTER(ctor_atoms_ptr, ctor_atoms_f, [n_ctor_atoms])
      allocate(ctor_atoms_copy(n_ctor_atoms))
      ctor_atoms_copy = ctor_atoms_f

      call C_F_POINTER(ctor_offsets_ptr, ctor_offsets_f, [n_ctor_offsets])
      allocate(ctor_offsets_copy(n_ctor_offsets))
      ctor_offsets_copy = ctor_offsets_f

      call C_F_POINTER(ctor_quads_ptr, ctor_quads_f, [4 * n_ctor])
      allocate(ctor_quads_copy(4 * n_ctor))
      ctor_quads_copy = ctor_quads_f
    else
      allocate(ctor_atoms_copy(0))
      allocate(ctor_offsets_copy(1))
      ctor_offsets_copy(1) = 0
      allocate(ctor_quads_copy(0))
      n_ctor_groups = 0
    end if

    ! Run analysis
    call analyze_zerocopy_com(mass_f, s_trajes_c, ana_period, &
                              dist_list_copy, n_dist, &
                              angl_list_copy, n_angl, &
                              tors_list_copy, n_tors, &
                              cdis_atoms_copy, cdis_offsets_copy, cdis_pairs_copy, &
                              n_cdis, n_cdis_groups, &
                              cang_atoms_copy, cang_offsets_copy, cang_triplets_copy, &
                              n_cang, n_cang_groups, &
                              ctor_atoms_copy, ctor_offsets_copy, ctor_quads_copy, &
                              n_ctor, n_ctor_groups, &
                              distance_ptr, angle_ptr, torsion_ptr, &
                              cdis_ptr, cang_ptr, ctor_ptr)

    ! Return pointers to results
    if (associated(distance_ptr)) then
      result_distance = c_loc(distance_ptr)
    end if
    if (associated(angle_ptr)) then
      result_angle = c_loc(angle_ptr)
    end if
    if (associated(torsion_ptr)) then
      result_torsion = c_loc(torsion_ptr)
    end if
    if (associated(cdis_ptr)) then
      result_cdis = c_loc(cdis_ptr)
    end if
    if (associated(cang_ptr)) then
      result_cang = c_loc(cang_ptr)
    end if
    if (associated(ctor_ptr)) then
      result_ctor = c_loc(ctor_ptr)
    end if

    ! Cleanup local arrays
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

  end subroutine trj_analysis_zerocopy_com_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    copy_error_msg
  !> @brief        Copy error message to C buffer
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine copy_error_msg(fortran_msg, c_msg, msglen)
    implicit none
    character(len=*), intent(in) :: fortran_msg
    character(kind=c_char), intent(out) :: c_msg(*)
    integer(c_int), value :: msglen

    integer :: i, copy_len

    copy_len = min(len_trim(fortran_msg), msglen - 1)
    do i = 1, copy_len
      c_msg(i) = fortran_msg(i:i)
    end do
    c_msg(copy_len + 1) = c_null_char

  end subroutine copy_error_msg

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_trj_results_internal
  !> @brief        Internal helper to deallocate TRJ analysis results
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine deallocate_trj_results_internal()
    implicit none

    if (associated(distance_ptr)) then
      deallocate(distance_ptr)
      nullify(distance_ptr)
    end if
    if (associated(angle_ptr)) then
      deallocate(angle_ptr)
      nullify(angle_ptr)
    end if
    if (associated(torsion_ptr)) then
      deallocate(torsion_ptr)
      nullify(torsion_ptr)
    end if
    if (associated(cdis_ptr)) then
      deallocate(cdis_ptr)
      nullify(cdis_ptr)
    end if
    if (associated(cang_ptr)) then
      deallocate(cang_ptr)
      nullify(cang_ptr)
    end if
    if (associated(ctor_ptr)) then
      deallocate(ctor_ptr)
      nullify(ctor_ptr)
    end if

  end subroutine deallocate_trj_results_internal

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_trj_results_c
  !> @brief        Deallocate TRJ analysis results (C interface)
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine deallocate_trj_results_c() bind(C, name="deallocate_trj_results_c")
    implicit none

    call deallocate_trj_results_internal()

  end subroutine deallocate_trj_results_c

  subroutine trj_analysis_main( &
          molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, &
          result_distance, result_angle, result_torsion, &
          result_cdis, result_cang, result_ctor)
    implicit none
    type(s_molecule), intent(inout) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer,                intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer,                intent(in) :: ctrl_len
    real(wp), pointer, intent(out) :: result_distance(:,:)
    real(wp), pointer, intent(out) :: result_angle(:,:)
    real(wp), pointer, intent(out) :: result_torsion(:,:)
    real(wp), pointer, intent(out) :: result_cdis(:,:)
    real(wp), pointer, intent(out) :: result_cang(:,:)
    real(wp), pointer, intent(out) :: result_ctor(:,:)

    ! local variables
    type(s_ctrl_data)      :: ctrl_data
    type(s_trajectory)     :: trajectory
    type(s_output)         :: output
    type(s_option)         :: option


    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! [Step1] Read control file
    !
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '

    call control_from_string(ctrl_text, ctrl_len, ctrl_data)


    ! [Step2] Set relevant variables and structures
    !
    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '

    call setup(molecule, ctrl_data, output, option)


    ! [Step3] Analyze trajectory
    !
    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '

    call analyze(molecule, s_trajes_c, ana_period, option, &
                 result_distance, result_angle, result_torsion, &
                 result_cdis, result_cang, result_ctor)


    ! [Step4] Deallocate memory
    !
    write(MsgOut,'(A)') '[STEP4] Deallocate memory'
    write(MsgOut,'(A)') ' '

    call dealloc_trajectory(trajectory)
  end subroutine trj_analysis_main

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in TRJ_ANALYSIS
  !! @authors      NT, TM
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] output     : output information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(molecule, ctrl_data, output, option)
    use ta_control_mod
    use ta_option_mod
    use ta_option_str_mod
    use trajectory_mod
    use output_mod
    use input_mod
    use trajectory_str_mod
    use output_str_mod
    use select_mod
    use molecules_mod
    use molecules_str_mod
    use fileio_grocrd_mod
    use fileio_grotop_mod
    use fileio_ambcrd_mod
    use fileio_prmtop_mod
    use fileio_psf_mod
    use fileio_pdb_mod
    implicit none

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_option),          intent(inout) :: option

    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)


    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)


    ! setup option
    !
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, output, option)

    return

  end subroutine setup

end module trj_c_mod
