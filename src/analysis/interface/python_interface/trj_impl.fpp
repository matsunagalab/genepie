!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ta_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Norio Takase (NT), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module trj_impl_mod

  use ta_option_str_mod
  use fileio_trj_mod
  use measure_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  public  :: analyze_zerocopy
  public  :: analyze_zerocopy_full
  public  :: analyze_zerocopy_com
  public  :: analyze_zerocopy_full_com
  private :: analyze_dis
  private :: analyze_ang
  private :: analyze_tor
  private :: analyze_comdis
  private :: analyze_comang
  private :: analyze_comtor
  private :: out_result

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @param[in]    molecule
  !! @param[in]    trajs_c
  !! @param[in]    ana_period
  !! @param[inout] option     : option information
  !! @param[out] distance
  !! @param[out] angle
  !! @param[out] torsion
  !! @param[out] cdis : com distance
  !! @param[out] cang : com angle
  !! @param[out] ctor : com torsion
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine analyze(molecule, trajes_c, ana_period, option, &
                     distance, angle, torsion, &
                     cdis, cang, ctor)
    use s_trajectories_c_mod

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period
    type(s_option),          intent(inout) :: option
    real(wp), pointer,       intent(out)   :: distance(:,:)
    real(wp), pointer,       intent(out)   :: angle(:,:)
    real(wp), pointer,       intent(out)   :: torsion(:,:)
    real(wp), pointer,       intent(out)   :: cdis(:,:)
    real(wp), pointer,       intent(out)   :: cang(:,:)
    real(wp), pointer,       intent(out)   :: ctor(:,:)


    ! local variables
    type(s_trajectory) :: trajectory
    integer            :: nstru, istep
    integer            :: dis_unit, ang_unit, tor_unit
    integer            :: cdis_unit, cang_unit, ctor_unit
    integer            :: num_out_frame


    if (option%check_only) &
      return

    num_out_frame = trajes_c%nframe / ana_period
    if (option%out_dis) then
      allocate( distance(size(option%distance), num_out_frame) )
    end if
    if (option%out_ang) then
      allocate( angle(size(option%angle), num_out_frame) )
    end if
    if (option%out_tor) then
      allocate( torsion(size(option%torsion), num_out_frame) )
    end if
    if (option%out_cdis) then
      allocate( cdis(size(option%cdistance), num_out_frame) )
    end if
    if (option%out_cang) then
      allocate( cang(size(option%cangle), num_out_frame) )
    end if
    if (option%out_ctor) then
      allocate( ctor(size(option%ctorsion), num_out_frame) )
    end if

    ! analysis loop
    !
    nstru = 0

    do istep = 1, trajes_c%nframe

      ! read trajectory
      !   coordinates of one MD snapshot are saved in trajectory%coord)
      !
      call get_frame(trajes_c, istep, trajectory)

      if (mod(istep, ana_period) == 0) then

        nstru = nstru + 1
        write(MsgOut,*) '      number of structures = ', nstru

        if (option%out_dis) then
          call analyze_dis(trajectory, option)
          distance(:, nstru) = option%distance
        end if
        if (option%out_ang) then
          call analyze_ang(trajectory, option)
          angle(:, nstru) = option%angle
        end if
        if (option%out_tor) then
          call analyze_tor(trajectory, option)
          torsion(:, nstru) = option%torsion
        end if
        if (option%out_cdis) then
          call analyze_comdis(molecule, trajectory, option)
          cdis(:, nstru) = option%cdistance
        end if
        if (option%out_cang) then
          call analyze_comang(molecule, trajectory, option)
          cang(:, nstru) = option%cangle
        end if
        if (option%out_ctor) then
          call analyze_comang(molecule, trajectory, option)
          ctor(:, nstru) = option%ctorsion
        end if
      end if
    end do

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_zerocopy
  !> @brief        Simplified trajectory analysis (zerocopy version)
  !! @authors      Claude Code
  !! @param[in]    trajes_c     : trajectories C structure
  !! @param[in]    ana_period   : analysis period
  !! @param[in]    dist_list    : distance atom pairs (2, n_dist)
  !! @param[in]    n_dist       : number of distance measurements
  !! @param[in]    angl_list    : angle atom triplets (3, n_angl)
  !! @param[in]    n_angl       : number of angle measurements
  !! @param[in]    tors_list    : torsion atom quadruplets (4, n_tors)
  !! @param[in]    n_tors       : number of torsion measurements
  !! @param[out]   distance     : distance results (n_dist, nframe)
  !! @param[out]   angle        : angle results (n_angl, nframe)
  !! @param[out]   torsion      : torsion results (n_tors, nframe)
  !! @note         This version handles simple atom-based measurements only.
  !!               For COM-based measurements, use the legacy function.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_zerocopy(trajes_c, ana_period, &
                              dist_list, n_dist, &
                              angl_list, n_angl, &
                              tors_list, n_tors, &
                              distance, angle, torsion)
    use s_trajectories_c_mod

    ! formal arguments
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period
    integer,                 intent(in)    :: dist_list(:,:)  ! (2, n_dist)
    integer,                 intent(in)    :: n_dist
    integer,                 intent(in)    :: angl_list(:,:)  ! (3, n_angl)
    integer,                 intent(in)    :: n_angl
    integer,                 intent(in)    :: tors_list(:,:)  ! (4, n_tors)
    integer,                 intent(in)    :: n_tors
    real(wp), pointer,       intent(out)   :: distance(:,:)
    real(wp), pointer,       intent(out)   :: angle(:,:)
    real(wp), pointer,       intent(out)   :: torsion(:,:)

    ! local variables
    type(s_trajectory) :: trajectory
    integer            :: nstru, istep, i
    integer            :: idx1, idx2, idx3, idx4
    integer            :: num_out_frame

    num_out_frame = trajes_c%nframe / ana_period

    ! Allocate output arrays
    if (n_dist > 0) then
      allocate(distance(n_dist, num_out_frame))
    else
      nullify(distance)
    end if

    if (n_angl > 0) then
      allocate(angle(n_angl, num_out_frame))
    else
      nullify(angle)
    end if

    if (n_tors > 0) then
      allocate(torsion(n_tors, num_out_frame))
    else
      nullify(torsion)
    end if

    ! Analysis loop
    nstru = 0

    do istep = 1, trajes_c%nframe

      ! Get trajectory frame
      call get_frame(trajes_c, istep, trajectory)

      if (mod(istep, ana_period) == 0) then

        nstru = nstru + 1
        write(MsgOut,*) '      number of structures = ', nstru

        ! Compute distances
        if (n_dist > 0) then
          do i = 1, n_dist
            idx1 = dist_list(1, i)
            idx2 = dist_list(2, i)
            distance(i, nstru) = compute_dis(trajectory%coord(:, idx1), &
                                             trajectory%coord(:, idx2))
          end do
        end if

        ! Compute angles
        if (n_angl > 0) then
          do i = 1, n_angl
            idx1 = angl_list(1, i)
            idx2 = angl_list(2, i)
            idx3 = angl_list(3, i)
            angle(i, nstru) = compute_ang(trajectory%coord(:, idx1), &
                                          trajectory%coord(:, idx2), &
                                          trajectory%coord(:, idx3))
          end do
        end if

        ! Compute torsions
        if (n_tors > 0) then
          do i = 1, n_tors
            idx1 = tors_list(1, i)
            idx2 = tors_list(2, i)
            idx3 = tors_list(3, i)
            idx4 = tors_list(4, i)
            torsion(i, nstru) = compute_dih(trajectory%coord(:, idx1), &
                                            trajectory%coord(:, idx2), &
                                            trajectory%coord(:, idx3), &
                                            trajectory%coord(:, idx4))
          end do
        end if

      end if

    end do

    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze_zerocopy> Trajectory analysis completed'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze_zerocopy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_zerocopy_full
  !> @brief        Trajectory analysis with full zero-copy (pre-allocated results)
  !! @authors      Claude Code
  !! @param[in]    trajes_c     : trajectories C structure
  !! @param[in]    ana_period   : analysis period
  !! @param[in]    dist_list    : distance atom pairs (2, n_dist)
  !! @param[in]    n_dist       : number of distance measurements
  !! @param[in]    angl_list    : angle atom triplets (3, n_angl)
  !! @param[in]    n_angl       : number of angle measurements
  !! @param[in]    tors_list    : torsion atom quadruplets (4, n_tors)
  !! @param[in]    n_tors       : number of torsion measurements
  !! @param[inout] distance     : pre-allocated distance results (n_dist, nframe)
  !! @param[inout] angle        : pre-allocated angle results (n_angl, nframe)
  !! @param[inout] torsion      : pre-allocated torsion results (n_tors, nframe)
  !! @param[out]   nstru_out    : actual number of structures analyzed
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_zerocopy_full(trajes_c, ana_period, &
                                   dist_list, n_dist, &
                                   angl_list, n_angl, &
                                   tors_list, n_tors, &
                                   distance, angle, torsion, nstru_out)
    use s_trajectories_c_mod

    ! formal arguments
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period
    integer,                 intent(in)    :: dist_list(:,:)  ! (2, n_dist)
    integer,                 intent(in)    :: n_dist
    integer,                 intent(in)    :: angl_list(:,:)  ! (3, n_angl)
    integer,                 intent(in)    :: n_angl
    integer,                 intent(in)    :: tors_list(:,:)  ! (4, n_tors)
    integer,                 intent(in)    :: n_tors
    real(wp),                intent(inout) :: distance(:,:)  ! pre-allocated
    real(wp),                intent(inout) :: angle(:,:)     ! pre-allocated
    real(wp),                intent(inout) :: torsion(:,:)   ! pre-allocated
    integer,                 intent(out)   :: nstru_out

    ! local variables
    type(s_trajectory) :: trajectory
    integer            :: nstru, istep, i
    integer            :: idx1, idx2, idx3, idx4

    ! Analysis loop (NO allocation - arrays are pre-allocated)
    nstru = 0

    do istep = 1, trajes_c%nframe

      ! Get trajectory frame
      call get_frame(trajes_c, istep, trajectory)

      if (mod(istep, ana_period) == 0) then

        nstru = nstru + 1
        write(MsgOut,*) '      number of structures = ', nstru

        ! Compute distances
        if (n_dist > 0) then
          do i = 1, n_dist
            idx1 = dist_list(1, i)
            idx2 = dist_list(2, i)
            distance(i, nstru) = compute_dis(trajectory%coord(:, idx1), &
                                             trajectory%coord(:, idx2))
          end do
        end if

        ! Compute angles
        if (n_angl > 0) then
          do i = 1, n_angl
            idx1 = angl_list(1, i)
            idx2 = angl_list(2, i)
            idx3 = angl_list(3, i)
            angle(i, nstru) = compute_ang(trajectory%coord(:, idx1), &
                                          trajectory%coord(:, idx2), &
                                          trajectory%coord(:, idx3))
          end do
        end if

        ! Compute torsions
        if (n_tors > 0) then
          do i = 1, n_tors
            idx1 = tors_list(1, i)
            idx2 = tors_list(2, i)
            idx3 = tors_list(3, i)
            idx4 = tors_list(4, i)
            torsion(i, nstru) = compute_dih(trajectory%coord(:, idx1), &
                                            trajectory%coord(:, idx2), &
                                            trajectory%coord(:, idx3), &
                                            trajectory%coord(:, idx4))
          end do
        end if

      end if

    end do

    nstru_out = nstru

    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze_zerocopy_full> Trajectory analysis completed (full zero-copy)'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze_zerocopy_full

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_zerocopy_com
  !> @brief        Trajectory analysis with COM calculations (zerocopy version)
  !! @authors      Claude Code
  !! @param[in]    mass         : atomic masses (n_atoms)
  !! @param[in]    trajes_c     : trajectories C structure
  !! @param[in]    ana_period   : analysis period
  !! @param[in]    dist_list    : distance atom pairs (2, n_dist)
  !! @param[in]    n_dist       : number of distance measurements
  !! @param[in]    angl_list    : angle atom triplets (3, n_angl)
  !! @param[in]    n_angl       : number of angle measurements
  !! @param[in]    tors_list    : torsion atom quadruplets (4, n_tors)
  !! @param[in]    n_tors       : number of torsion measurements
  !! @param[in]    cdis_atoms   : flat array of atom indices for COM distance groups
  !! @param[in]    cdis_offsets : offsets into cdis_atoms for each group (n_cdis_groups+1)
  !! @param[in]    cdis_pairs   : group pair indices for COM distances (2*n_cdis)
  !! @param[in]    n_cdis       : number of COM distance measurements
  !! @param[in]    n_cdis_groups: number of COM distance groups
  !! @param[in]    cang_atoms   : flat array of atom indices for COM angle groups
  !! @param[in]    cang_offsets : offsets into cang_atoms for each group
  !! @param[in]    cang_triplets: group triplet indices for COM angles (3*n_cang)
  !! @param[in]    n_cang       : number of COM angle measurements
  !! @param[in]    n_cang_groups: number of COM angle groups
  !! @param[in]    ctor_atoms   : flat array of atom indices for COM torsion groups
  !! @param[in]    ctor_offsets : offsets into ctor_atoms for each group
  !! @param[in]    ctor_quads   : group quad indices for COM torsions (4*n_ctor)
  !! @param[in]    n_ctor       : number of COM torsion measurements
  !! @param[in]    n_ctor_groups: number of COM torsion groups
  !! @param[out]   distance     : distance results (n_dist, nframe)
  !! @param[out]   angle        : angle results (n_angl, nframe)
  !! @param[out]   torsion      : torsion results (n_tors, nframe)
  !! @param[out]   cdis         : COM distance results (n_cdis, nframe)
  !! @param[out]   cang         : COM angle results (n_cang, nframe)
  !! @param[out]   ctor         : COM torsion results (n_ctor, nframe)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_zerocopy_com(mass, trajes_c, ana_period, &
                                  dist_list, n_dist, &
                                  angl_list, n_angl, &
                                  tors_list, n_tors, &
                                  cdis_atoms, cdis_offsets, cdis_pairs, &
                                  n_cdis, n_cdis_groups, &
                                  cang_atoms, cang_offsets, cang_triplets, &
                                  n_cang, n_cang_groups, &
                                  ctor_atoms, ctor_offsets, ctor_quads, &
                                  n_ctor, n_ctor_groups, &
                                  distance, angle, torsion, &
                                  cdis, cang, ctor)
    use s_trajectories_c_mod

    ! formal arguments
    real(wp),                intent(in)    :: mass(:)
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period

    ! Atom-based measurements
    integer,                 intent(in)    :: dist_list(:,:)  ! (2, n_dist)
    integer,                 intent(in)    :: n_dist
    integer,                 intent(in)    :: angl_list(:,:)  ! (3, n_angl)
    integer,                 intent(in)    :: n_angl
    integer,                 intent(in)    :: tors_list(:,:)  ! (4, n_tors)
    integer,                 intent(in)    :: n_tors

    ! COM distance: flat atoms + offsets + pairs
    integer,                 intent(in)    :: cdis_atoms(:)
    integer,                 intent(in)    :: cdis_offsets(:)  ! (n_cdis_groups + 1)
    integer,                 intent(in)    :: cdis_pairs(:)    ! (2 * n_cdis)
    integer,                 intent(in)    :: n_cdis
    integer,                 intent(in)    :: n_cdis_groups

    ! COM angle: flat atoms + offsets + triplets
    integer,                 intent(in)    :: cang_atoms(:)
    integer,                 intent(in)    :: cang_offsets(:)  ! (n_cang_groups + 1)
    integer,                 intent(in)    :: cang_triplets(:) ! (3 * n_cang)
    integer,                 intent(in)    :: n_cang
    integer,                 intent(in)    :: n_cang_groups

    ! COM torsion: flat atoms + offsets + quads
    integer,                 intent(in)    :: ctor_atoms(:)
    integer,                 intent(in)    :: ctor_offsets(:)  ! (n_ctor_groups + 1)
    integer,                 intent(in)    :: ctor_quads(:)    ! (4 * n_ctor)
    integer,                 intent(in)    :: n_ctor
    integer,                 intent(in)    :: n_ctor_groups

    ! Output arrays
    real(wp), pointer,       intent(out)   :: distance(:,:)
    real(wp), pointer,       intent(out)   :: angle(:,:)
    real(wp), pointer,       intent(out)   :: torsion(:,:)
    real(wp), pointer,       intent(out)   :: cdis(:,:)
    real(wp), pointer,       intent(out)   :: cang(:,:)
    real(wp), pointer,       intent(out)   :: ctor(:,:)

    ! local variables
    type(s_trajectory) :: trajectory
    integer            :: nstru, istep, i
    integer            :: idx1, idx2, idx3, idx4
    integer            :: grp1, grp2, grp3, grp4
    integer            :: start1, end1, start2, end2, start3, end3, start4, end4
    integer            :: num_out_frame
    real(wp)           :: com1(3), com2(3), com3(3), com4(3)

    num_out_frame = trajes_c%nframe / ana_period

    ! Allocate output arrays for atom-based measurements
    if (n_dist > 0) then
      allocate(distance(n_dist, num_out_frame))
    else
      nullify(distance)
    end if

    if (n_angl > 0) then
      allocate(angle(n_angl, num_out_frame))
    else
      nullify(angle)
    end if

    if (n_tors > 0) then
      allocate(torsion(n_tors, num_out_frame))
    else
      nullify(torsion)
    end if

    ! Allocate output arrays for COM-based measurements
    if (n_cdis > 0) then
      allocate(cdis(n_cdis, num_out_frame))
    else
      nullify(cdis)
    end if

    if (n_cang > 0) then
      allocate(cang(n_cang, num_out_frame))
    else
      nullify(cang)
    end if

    if (n_ctor > 0) then
      allocate(ctor(n_ctor, num_out_frame))
    else
      nullify(ctor)
    end if

    ! Analysis loop
    nstru = 0

    do istep = 1, trajes_c%nframe

      ! Get trajectory frame
      call get_frame(trajes_c, istep, trajectory)

      if (mod(istep, ana_period) == 0) then

        nstru = nstru + 1
        write(MsgOut,*) '      number of structures = ', nstru

        ! ===== Atom-based measurements =====

        ! Compute distances
        if (n_dist > 0) then
          do i = 1, n_dist
            idx1 = dist_list(1, i)
            idx2 = dist_list(2, i)
            distance(i, nstru) = compute_dis(trajectory%coord(:, idx1), &
                                             trajectory%coord(:, idx2))
          end do
        end if

        ! Compute angles
        if (n_angl > 0) then
          do i = 1, n_angl
            idx1 = angl_list(1, i)
            idx2 = angl_list(2, i)
            idx3 = angl_list(3, i)
            angle(i, nstru) = compute_ang(trajectory%coord(:, idx1), &
                                          trajectory%coord(:, idx2), &
                                          trajectory%coord(:, idx3))
          end do
        end if

        ! Compute torsions
        if (n_tors > 0) then
          do i = 1, n_tors
            idx1 = tors_list(1, i)
            idx2 = tors_list(2, i)
            idx3 = tors_list(3, i)
            idx4 = tors_list(4, i)
            torsion(i, nstru) = compute_dih(trajectory%coord(:, idx1), &
                                            trajectory%coord(:, idx2), &
                                            trajectory%coord(:, idx3), &
                                            trajectory%coord(:, idx4))
          end do
        end if

        ! ===== COM-based measurements =====

        ! Compute COM distances
        if (n_cdis > 0) then
          do i = 1, n_cdis
            ! Get group indices (0-based from Python, convert to 1-based)
            grp1 = cdis_pairs(2*i - 1) + 1
            grp2 = cdis_pairs(2*i) + 1

            ! Get atom range for group 1 (offsets are 0-based)
            start1 = cdis_offsets(grp1) + 1
            end1 = cdis_offsets(grp1 + 1)

            ! Get atom range for group 2
            start2 = cdis_offsets(grp2) + 1
            end2 = cdis_offsets(grp2 + 1)

            ! Compute COMs
            com1 = compute_com(trajectory%coord, mass, cdis_atoms(start1:end1))
            com2 = compute_com(trajectory%coord, mass, cdis_atoms(start2:end2))

            ! Compute distance
            cdis(i, nstru) = compute_dis(com1, com2)
          end do
        end if

        ! Compute COM angles
        if (n_cang > 0) then
          do i = 1, n_cang
            ! Get group indices (0-based from Python, convert to 1-based)
            grp1 = cang_triplets(3*i - 2) + 1
            grp2 = cang_triplets(3*i - 1) + 1
            grp3 = cang_triplets(3*i) + 1

            ! Get atom ranges
            start1 = cang_offsets(grp1) + 1
            end1 = cang_offsets(grp1 + 1)
            start2 = cang_offsets(grp2) + 1
            end2 = cang_offsets(grp2 + 1)
            start3 = cang_offsets(grp3) + 1
            end3 = cang_offsets(grp3 + 1)

            ! Compute COMs
            com1 = compute_com(trajectory%coord, mass, cang_atoms(start1:end1))
            com2 = compute_com(trajectory%coord, mass, cang_atoms(start2:end2))
            com3 = compute_com(trajectory%coord, mass, cang_atoms(start3:end3))

            ! Compute angle
            cang(i, nstru) = compute_ang(com1, com2, com3)
          end do
        end if

        ! Compute COM torsions
        if (n_ctor > 0) then
          do i = 1, n_ctor
            ! Get group indices (0-based from Python, convert to 1-based)
            grp1 = ctor_quads(4*i - 3) + 1
            grp2 = ctor_quads(4*i - 2) + 1
            grp3 = ctor_quads(4*i - 1) + 1
            grp4 = ctor_quads(4*i) + 1

            ! Get atom ranges
            start1 = ctor_offsets(grp1) + 1
            end1 = ctor_offsets(grp1 + 1)
            start2 = ctor_offsets(grp2) + 1
            end2 = ctor_offsets(grp2 + 1)
            start3 = ctor_offsets(grp3) + 1
            end3 = ctor_offsets(grp3 + 1)
            start4 = ctor_offsets(grp4) + 1
            end4 = ctor_offsets(grp4 + 1)

            ! Compute COMs
            com1 = compute_com(trajectory%coord, mass, ctor_atoms(start1:end1))
            com2 = compute_com(trajectory%coord, mass, ctor_atoms(start2:end2))
            com3 = compute_com(trajectory%coord, mass, ctor_atoms(start3:end3))
            com4 = compute_com(trajectory%coord, mass, ctor_atoms(start4:end4))

            ! Compute torsion
            ctor(i, nstru) = compute_dih(com1, com2, com3, com4)
          end do
        end if

      end if

    end do

    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze_zerocopy_com> Trajectory analysis with COM completed'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze_zerocopy_com

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_zerocopy_full_com
  !> @brief        Trajectory analysis with COM (full zerocopy - pre-allocated)
  !! @authors      Claude Code
  !! @param[in]    mass         : atomic masses (n_atoms)
  !! @param[in]    trajes_c     : trajectories C structure
  !! @param[in]    ana_period   : analysis period
  !! @param[in]    dist_list    : distance atom pairs (2, n_dist)
  !! @param[in]    n_dist       : number of distance measurements
  !! @param[in]    angl_list    : angle atom triplets (3, n_angl)
  !! @param[in]    n_angl       : number of angle measurements
  !! @param[in]    tors_list    : torsion atom quadruplets (4, n_tors)
  !! @param[in]    n_tors       : number of torsion measurements
  !! @param[in]    cdis_atoms   : flat array of atom indices for COM distance
  !! @param[in]    cdis_offsets : offsets into cdis_atoms for each group
  !! @param[in]    cdis_pairs   : group pair indices for COM distances
  !! @param[in]    n_cdis       : number of COM distance measurements
  !! @param[in]    n_cdis_groups: number of COM distance groups
  !! @param[in]    cang_atoms   : flat array of atom indices for COM angles
  !! @param[in]    cang_offsets : offsets into cang_atoms for each group
  !! @param[in]    cang_triplets: group triplet indices for COM angles
  !! @param[in]    n_cang       : number of COM angle measurements
  !! @param[in]    n_cang_groups: number of COM angle groups
  !! @param[in]    ctor_atoms   : flat array of atom indices for COM torsions
  !! @param[in]    ctor_offsets : offsets into ctor_atoms for each group
  !! @param[in]    ctor_quads   : group quad indices for COM torsions
  !! @param[in]    n_ctor       : number of COM torsion measurements
  !! @param[in]    n_ctor_groups: number of COM torsion groups
  !! @param[inout] distance     : pre-allocated distance results (n_dist, nframe)
  !! @param[inout] angle        : pre-allocated angle results (n_angl, nframe)
  !! @param[inout] torsion      : pre-allocated torsion results (n_tors, nframe)
  !! @param[inout] cdis         : pre-allocated COM distance results (n_cdis, nframe)
  !! @param[inout] cang         : pre-allocated COM angle results (n_cang, nframe)
  !! @param[inout] ctor         : pre-allocated COM torsion results (n_ctor, nframe)
  !! @param[out]   nstru_out    : actual number of frames processed
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_zerocopy_full_com(mass, trajes_c, ana_period, &
                                       dist_list, n_dist, &
                                       angl_list, n_angl, &
                                       tors_list, n_tors, &
                                       cdis_atoms, cdis_offsets, cdis_pairs, &
                                       n_cdis, n_cdis_groups, &
                                       cang_atoms, cang_offsets, cang_triplets, &
                                       n_cang, n_cang_groups, &
                                       ctor_atoms, ctor_offsets, ctor_quads, &
                                       n_ctor, n_ctor_groups, &
                                       distance, angle, torsion, &
                                       cdis, cang, ctor, nstru_out)
    use s_trajectories_c_mod

    ! formal arguments
    real(wp),                intent(in)    :: mass(:)
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period

    ! Atom-based measurements
    integer,                 intent(in)    :: dist_list(:,:)  ! (2, n_dist)
    integer,                 intent(in)    :: n_dist
    integer,                 intent(in)    :: angl_list(:,:)  ! (3, n_angl)
    integer,                 intent(in)    :: n_angl
    integer,                 intent(in)    :: tors_list(:,:)  ! (4, n_tors)
    integer,                 intent(in)    :: n_tors

    ! COM distance: flat atoms + offsets + pairs
    integer,                 intent(in)    :: cdis_atoms(:)
    integer,                 intent(in)    :: cdis_offsets(:)  ! (n_cdis_groups + 1)
    integer,                 intent(in)    :: cdis_pairs(:)    ! (2 * n_cdis)
    integer,                 intent(in)    :: n_cdis
    integer,                 intent(in)    :: n_cdis_groups

    ! COM angle: flat atoms + offsets + triplets
    integer,                 intent(in)    :: cang_atoms(:)
    integer,                 intent(in)    :: cang_offsets(:)  ! (n_cang_groups + 1)
    integer,                 intent(in)    :: cang_triplets(:) ! (3 * n_cang)
    integer,                 intent(in)    :: n_cang
    integer,                 intent(in)    :: n_cang_groups

    ! COM torsion: flat atoms + offsets + quads
    integer,                 intent(in)    :: ctor_atoms(:)
    integer,                 intent(in)    :: ctor_offsets(:)  ! (n_ctor_groups + 1)
    integer,                 intent(in)    :: ctor_quads(:)    ! (4 * n_ctor)
    integer,                 intent(in)    :: n_ctor
    integer,                 intent(in)    :: n_ctor_groups

    ! Pre-allocated output arrays
    real(wp),                intent(inout) :: distance(:,:)
    real(wp),                intent(inout) :: angle(:,:)
    real(wp),                intent(inout) :: torsion(:,:)
    real(wp),                intent(inout) :: cdis(:,:)
    real(wp),                intent(inout) :: cang(:,:)
    real(wp),                intent(inout) :: ctor(:,:)
    integer,                 intent(out)   :: nstru_out

    ! local variables
    type(s_trajectory) :: trajectory
    integer            :: nstru, istep, i
    integer            :: idx1, idx2, idx3, idx4
    integer            :: grp1, grp2, grp3, grp4
    integer            :: start1, end1, start2, end2, start3, end3, start4, end4
    real(wp)           :: com1(3), com2(3), com3(3), com4(3)

    ! Analysis loop
    nstru = 0

    do istep = 1, trajes_c%nframe

      ! Get trajectory frame
      call get_frame(trajes_c, istep, trajectory)

      if (mod(istep, ana_period) == 0) then

        nstru = nstru + 1
        write(MsgOut,*) '      number of structures = ', nstru

        ! ===== Atom-based measurements =====

        ! Compute distances
        if (n_dist > 0) then
          do i = 1, n_dist
            idx1 = dist_list(1, i)
            idx2 = dist_list(2, i)
            distance(i, nstru) = compute_dis(trajectory%coord(:, idx1), &
                                             trajectory%coord(:, idx2))
          end do
        end if

        ! Compute angles
        if (n_angl > 0) then
          do i = 1, n_angl
            idx1 = angl_list(1, i)
            idx2 = angl_list(2, i)
            idx3 = angl_list(3, i)
            angle(i, nstru) = compute_ang(trajectory%coord(:, idx1), &
                                          trajectory%coord(:, idx2), &
                                          trajectory%coord(:, idx3))
          end do
        end if

        ! Compute torsions
        if (n_tors > 0) then
          do i = 1, n_tors
            idx1 = tors_list(1, i)
            idx2 = tors_list(2, i)
            idx3 = tors_list(3, i)
            idx4 = tors_list(4, i)
            torsion(i, nstru) = compute_dih(trajectory%coord(:, idx1), &
                                            trajectory%coord(:, idx2), &
                                            trajectory%coord(:, idx3), &
                                            trajectory%coord(:, idx4))
          end do
        end if

        ! ===== COM-based measurements =====

        ! Compute COM distances
        if (n_cdis > 0) then
          do i = 1, n_cdis
            ! Get group indices (0-based from Python, convert to 1-based)
            grp1 = cdis_pairs(2*i - 1) + 1
            grp2 = cdis_pairs(2*i) + 1

            ! Get atom range for group 1 (offsets are 0-based)
            start1 = cdis_offsets(grp1) + 1
            end1 = cdis_offsets(grp1 + 1)

            ! Get atom range for group 2
            start2 = cdis_offsets(grp2) + 1
            end2 = cdis_offsets(grp2 + 1)

            ! Compute COMs
            com1 = compute_com(trajectory%coord, mass, cdis_atoms(start1:end1))
            com2 = compute_com(trajectory%coord, mass, cdis_atoms(start2:end2))

            ! Compute distance
            cdis(i, nstru) = compute_dis(com1, com2)
          end do
        end if

        ! Compute COM angles
        if (n_cang > 0) then
          do i = 1, n_cang
            ! Get group indices (0-based from Python, convert to 1-based)
            grp1 = cang_triplets(3*i - 2) + 1
            grp2 = cang_triplets(3*i - 1) + 1
            grp3 = cang_triplets(3*i) + 1

            ! Get atom ranges
            start1 = cang_offsets(grp1) + 1
            end1 = cang_offsets(grp1 + 1)
            start2 = cang_offsets(grp2) + 1
            end2 = cang_offsets(grp2 + 1)
            start3 = cang_offsets(grp3) + 1
            end3 = cang_offsets(grp3 + 1)

            ! Compute COMs
            com1 = compute_com(trajectory%coord, mass, cang_atoms(start1:end1))
            com2 = compute_com(trajectory%coord, mass, cang_atoms(start2:end2))
            com3 = compute_com(trajectory%coord, mass, cang_atoms(start3:end3))

            ! Compute angle
            cang(i, nstru) = compute_ang(com1, com2, com3)
          end do
        end if

        ! Compute COM torsions
        if (n_ctor > 0) then
          do i = 1, n_ctor
            ! Get group indices (0-based from Python, convert to 1-based)
            grp1 = ctor_quads(4*i - 3) + 1
            grp2 = ctor_quads(4*i - 2) + 1
            grp3 = ctor_quads(4*i - 1) + 1
            grp4 = ctor_quads(4*i) + 1

            ! Get atom ranges
            start1 = ctor_offsets(grp1) + 1
            end1 = ctor_offsets(grp1 + 1)
            start2 = ctor_offsets(grp2) + 1
            end2 = ctor_offsets(grp2 + 1)
            start3 = ctor_offsets(grp3) + 1
            end3 = ctor_offsets(grp3 + 1)
            start4 = ctor_offsets(grp4) + 1
            end4 = ctor_offsets(grp4 + 1)

            ! Compute COMs
            com1 = compute_com(trajectory%coord, mass, ctor_atoms(start1:end1))
            com2 = compute_com(trajectory%coord, mass, ctor_atoms(start2:end2))
            com3 = compute_com(trajectory%coord, mass, ctor_atoms(start3:end3))
            com4 = compute_com(trajectory%coord, mass, ctor_atoms(start4:end4))

            ! Compute torsion
            ctor(i, nstru) = compute_dih(com1, com2, com3, com4)
          end do
        end if

      end if

    end do

    nstru_out = nstru

    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze_zerocopy_full_com> Trajectory analysis completed'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze_zerocopy_full_com

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_dis
  !> @brief        analyze distances
  !! @authors      NT, TM, SI
  !! @param[in]    trajectory : trajectory information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_dis(trajectory, option)

    ! formal argments
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: i, j, idx1, idx2


    do i = 1, size(option%distance)
      option%distance(i) = 0.0_wp
      
      do j = 1, option%dist_num(i)/2
        idx1 = option%dist_list(2*j-1,i)
        idx2 = option%dist_list(2*j,i)

        option%distance(i) = option%distance(i) + option%dist_weight(i,j) &
                           * compute_dis(trajectory%coord(:,idx1), &
                                          trajectory%coord(:,idx2))
      end do

    end do

    return

  end subroutine analyze_dis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_ang
  !> @brief        analize angles
  !! @authors      NT, TM
  !! @param[in]    trajectory : trajectory information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_ang(trajectory, option)

    ! formal argments
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: i, idx1, idx2, idx3


    do i = 1, size(option%angle)
      
      idx1 = option%angl_list(1, i)
      idx2 = option%angl_list(2, i)
      idx3 = option%angl_list(3, i)

      option%angle(i) = compute_ang(trajectory%coord(:,idx1), &
                                    trajectory%coord(:,idx2), &
                                    trajectory%coord(:,idx3))

    end do

    return

  end subroutine analyze_ang

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_tor
  !> @brief        analize torsions
  !! @authors      NT, TM
  !! @param[in]    trajectory : trajectory informatin
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_tor(trajectory, option)

    ! formal argments
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_option),          intent(inout) :: option

    ! local variables
    integer                  :: i, idx1, idx2, idx3, idx4


    do i = 1, size(option%torsion)
      
      idx1 = option%tors_list(1, i)
      idx2 = option%tors_list(2, i)
      idx3 = option%tors_list(3, i)
      idx4 = option%tors_list(4, i)

      option%torsion(i) = compute_dih(trajectory%coord(:,idx1), &
                                      trajectory%coord(:,idx2), &
                                      trajectory%coord(:,idx3), &
                                      trajectory%coord(:,idx4))

    end do

    return

  end subroutine analyze_tor

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_comdis
  !> @brief        analyze COM distances
  !! @authors      NT
  !! @param[in]    molecule   : molecule information
  !! @param[in]    trajectory : trajectory information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_comdis(molecule, trajectory, option)

    ! formal argments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_option),          intent(inout) :: option

    ! local variables
    real(wp)                 :: c1(3), c2(3)
    integer                  :: i


    do i = 1, size(option%cdist_group(1,:))

      ! atom group 1
      c1 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%cdist_group(1,i))%idx)

      ! atom group 2
      c2 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%cdist_group(2,i))%idx)

      ! compute distance
      option%cdistance(i) = compute_dis(c1, c2)

    end do

    return

  end subroutine analyze_comdis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_comang
  !> @brief        analyze COM angle
  !! @authors      NT
  !! @param[in]    molecule   : molecule information
  !! @param[in]    trajectory : trajectory information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_comang(molecule, trajectory, option)

    ! formal argments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_option), target,  intent(inout) :: option

    ! local variables
    real(wp)                 :: c1(3), c2(3), c3(3)
    integer                  :: i


    do i = 1, size(option%cangl_group(1,:))

      ! atom group 1
      c1 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%cangl_group(1,i))%idx)

      ! atom group 2
      c2 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%cangl_group(2,i))%idx)

      ! atom group 3
      c3 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%cangl_group(3,i))%idx)

      ! compute angle
      option%cangle(i) = compute_ang(c1, c2, c3)

    end do

    return

  end subroutine analyze_comang

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_comtor
  !> @brief        analyze COM torsion
  !! @authors      NT
  !! @param[in]    molecule   : molecule information
  !! @param[in]    trajectory : trajectory information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_comtor(molecule, trajectory, option)

    ! formal argments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectory),      intent(in)    :: trajectory
    type(s_option),          intent(inout) :: option

    ! local variables
    real(wp)                 :: c1(3), c2(3), c3(3), c4(3)
    integer                  :: i


    do i = 1, size(option%ctor_group(1,:))

      ! atom group 1
      c1 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%ctor_group(1,i))%idx)

      ! atom group 2
      c2 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%ctor_group(2,i))%idx)

      ! atom group 3
      c3 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%ctor_group(3,i))%idx)

      ! atom group 4
      c4 = compute_com(trajectory%coord, molecule%mass, &
                       option%selatoms(option%ctor_group(4,i))%idx)

      ! compute angle
      option%ctorsion(i) = compute_dih(c1, c2, c3, c4)

    end do

    return

  end subroutine analyze_comtor

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    out_result
  !> @brief        output analysis results
  !! @authors      NT
  !! @param[in]    sturct_no : sturcture number
  !! @param[in]    unit_no   : file unit number
  !! @param[in]    results   : analysis results
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine out_result(struct_no, unit_no, results)

    ! formal arguments
    integer,                 intent(in)    :: struct_no
    integer,                 intent(in)    :: unit_no
    real(wp),                intent(in)    :: results(:)

    ! local variables
    integer                  :: i


    write(unit_no, '(i10,1x,100(f8.3,1x))') &
         struct_no, (results(i), i=1, size(results))

    return

  end subroutine out_result

end module trj_impl_mod
