!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   trj_impl_mod
!> @brief   Trajectory analysis implementation (zerocopy)
!! @authors Norio Takase (NT), Takaharu Mori (TM), Claude Code
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module trj_impl_mod

  use fileio_trj_mod
  use measure_mod
  use trajectory_str_mod
  use molecules_str_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  public  :: analyze_com

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        Trajectory analysis (zerocopy, pre-allocated results)
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

  subroutine analyze(trajes_c, ana_period, &
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
    write(MsgOut,'(A)') 'Analyze> Trajectory analysis completed'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_com
  !> @brief        Trajectory analysis with COM (zerocopy, pre-allocated)
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

  subroutine analyze_com(mass, trajes_c, ana_period, &
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
    write(MsgOut,'(A)') 'Analyze_com> Trajectory analysis with COM completed'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze_com

end module trj_impl_mod
