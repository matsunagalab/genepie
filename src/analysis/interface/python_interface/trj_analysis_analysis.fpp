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

module trj_analysis_analyze_c_mod

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
  !! @authors      NT, TM
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trajes_c, ana_period, option, &
                     distance, num_distance, &
                     angle, num_angle, &
                     torsion, num_torsion, &
                     cdis, num_cdis, &
                     cang, num_cang, &
                     ctor, num_ctor)
    use s_trajectories_c_mod

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period
    type(s_option),          intent(inout) :: option
    real(wp), pointer,       intent(out)   :: distance(:,:)
    integer,                 intent(out)   :: num_distance
    real(wp), pointer,       intent(out)   :: angle(:,:)
    integer,                 intent(out)   :: num_angle
    real(wp), pointer,       intent(out)   :: torsion(:,:)
    integer,                 intent(out)   :: num_torsion
    real(wp), pointer,       intent(out)   :: cdis(:,:)
    integer,                 intent(out)   :: num_cdis
    real(wp), pointer,       intent(out)   :: cang(:,:)
    integer,                 intent(out)   :: num_cang
    real(wp), pointer,       intent(out)   :: ctor(:,:)
    integer,                 intent(out)   :: num_ctor


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
      num_distance = size(option%distance)
      allocate( distance(num_distance, num_out_frame) )
    end if
    if (option%out_ang) then
      num_angle = size(option%angle)
      allocate( angle(num_angle, num_out_frame) )
    end if
    if (option%out_tor) then
      num_torsion = size(option%torsion)
      allocate( torsion(num_torsion, num_out_frame) )
    end if
    if (option%out_cdis) then
      num_cdis = size(option%cdistance)
      allocate( cdis(num_cdis, num_out_frame) )
    end if
    if (option%out_cang) then
      num_cang = size(option%cangle)
      allocate( cang(num_cang, num_out_frame) )
    end if
    if (option%out_ctor) then
      num_ctor = size(option%ctorsion)
      allocate( ctor(num_ctor, num_out_frame) )
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

end module trj_analysis_analyze_c_mod
