!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ra_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rmsd_impl_mod

  use fileio_trj_mod
  use fitting_mod
  use fitting_str_mod
  use measure_mod
  use trajectory_str_mod
  use molecules_str_mod
  use error_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  public  :: analyze_with_fitting

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        RMSD analysis (zerocopy, pre-allocated result array)
  !! @authors      Claude Code
  !! @param[in]    mass           : mass array pointer (view of Python NumPy)
  !! @param[in]    ref_coord      : reference coordinates (3, n_atoms)
  !! @param[in]    trajes_c       : trajectory C structure
  !! @param[in]    ana_period     : analysis period
  !! @param[in]    analysis_idx   : analysis atom indices (1-indexed)
  !! @param[in]    n_analysis     : number of analysis atoms
  !! @param[in]    mass_weighted  : use mass weighting
  !! @param[inout] rmsd_results   : pre-allocated RMSD results array
  !! @param[out]   nstru_out      : actual number of structures analyzed
  !! @note         This version does NOT perform fitting. Result array must be
  !!               pre-allocated by caller.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(mass, ref_coord, trajes_c, ana_period, &
                     analysis_idx, n_analysis, mass_weighted, &
                     rmsd_results, nstru_out)
    use s_trajectories_c_mod

    ! formal arguments
    real(wp), pointer,       intent(in)    :: mass(:)
    real(wp), pointer,       intent(in)    :: ref_coord(:,:)
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period
    integer,                 intent(in)    :: analysis_idx(:)
    integer,                 intent(in)    :: n_analysis
    logical,                 intent(in)    :: mass_weighted
    real(wp),                intent(inout) :: rmsd_results(:)
    integer,                 intent(out)   :: nstru_out

    ! local variables
    type(s_trajectory) :: trajectory
    integer            :: nstru, istep
    integer            :: iatom, idx
    real(wp)           :: rmsd, tot_mass, weight

    ! analysis loop (NO allocation - results array is pre-allocated)
    nstru = 0

    do istep = 1, trajes_c%nframe

      ! read trajectory frame
      call get_frame(trajes_c, istep, trajectory)

      if (mod(istep, ana_period) == 0) then

        nstru = nstru + 1
        write(MsgOut,*) '      number of structures = ', nstru

        ! compute RMSD
        rmsd = 0.0_wp
        tot_mass = 0.0_wp

        do iatom = 1, n_analysis
          idx = analysis_idx(iatom)

          if (mass_weighted) then
            weight = mass(idx)
          else
            weight = 1.0_wp
          end if

          rmsd = rmsd + weight * ( &
            (ref_coord(1, idx) - trajectory%coord(1, idx))**2 + &
            (ref_coord(2, idx) - trajectory%coord(2, idx))**2 + &
            (ref_coord(3, idx) - trajectory%coord(3, idx))**2)

          tot_mass = tot_mass + weight
        end do

        if (tot_mass > EPS) rmsd = sqrt(rmsd / tot_mass)

        rmsd_results(nstru) = rmsd

        ! output results
        write(MsgOut,'(a,f10.5)') '              RMSD of analysis atoms = ', rmsd
        write(MsgOut,*) ''

      end if

    end do

    nstru_out = nstru

    ! Output summary
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> RMSD analysis completed (no fitting)'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_with_fitting
  !> @brief        RMSD analysis with fitting (zerocopy, pre-allocated)
  !! @authors      Claude Code
  !! @param[in]    mass           : mass array pointer (view of Python NumPy)
  !! @param[in]    ref_coord      : reference coordinates (3, n_atoms)
  !! @param[in]    n_atoms        : number of atoms
  !! @param[in]    trajes_c       : trajectory C structure
  !! @param[in]    ana_period     : analysis period
  !! @param[in]    fitting_idx    : fitting atom indices (1-indexed)
  !! @param[in]    n_fitting      : number of fitting atoms
  !! @param[in]    analysis_idx   : analysis atom indices (1-indexed)
  !! @param[in]    n_analysis     : number of analysis atoms
  !! @param[in]    fitting_method : fitting method (1-6)
  !! @param[in]    mass_weighted  : use mass weighting
  !! @param[inout] rmsd_results   : pre-allocated RMSD results array
  !! @param[out]   nstru_out      : number of structures analyzed
  !! @note         This version performs structural fitting before RMSD calculation.
  !!               Result array must be pre-allocated by caller.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_with_fitting(mass, ref_coord, n_atoms, &
                              trajes_c, ana_period, &
                              fitting_idx, n_fitting, &
                              analysis_idx, n_analysis, &
                              fitting_method, mass_weighted, &
                              rmsd_results, nstru_out)
    use s_trajectories_c_mod

    ! formal arguments
    real(wp), pointer,       intent(in)    :: mass(:)
    real(wp), pointer,       intent(in)    :: ref_coord(:,:)
    integer,                 intent(in)    :: n_atoms
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period
    integer,                 intent(in)    :: fitting_idx(:)
    integer,                 intent(in)    :: n_fitting
    integer,                 intent(in)    :: analysis_idx(:)
    integer,                 intent(in)    :: n_analysis
    integer,                 intent(in)    :: fitting_method
    logical,                 intent(in)    :: mass_weighted
    real(wp),                intent(inout) :: rmsd_results(:)
    integer,                 intent(out)   :: nstru_out

    ! local variables
    type(s_trajectory) :: trajectory
    integer            :: nstru, istep
    integer            :: iatom, idx, i
    real(wp)           :: rmsd, tot_mass, weight

    ! Fitting-related variables
    real(wp), allocatable :: coord_work(:,:)     ! work buffer (3, n_atoms)
    real(wp), allocatable :: mass_fitting(:)     ! mass for fitting atoms
    real(wp)              :: rot_matrix(3,3)
    real(wp)              :: com_ref(3), com_mov(3)
    real(wp)              :: fit_rmsd
    integer               :: ierr

    ! allocate work buffer (to avoid modifying original trajectory)
    allocate(coord_work(3, n_atoms))

    ! prepare mass array for fitting
    allocate(mass_fitting(n_atoms))
    if (mass_weighted) then
      mass_fitting(:) = mass(:)
    else
      mass_fitting(:) = 1.0_wp
    end if

    ! analysis loop (NO allocation of results - pre-allocated)
    nstru = 0

    do istep = 1, trajes_c%nframe

      ! read trajectory frame
      call get_frame(trajes_c, istep, trajectory)

      if (mod(istep, ana_period) == 0) then

        nstru = nstru + 1
        write(MsgOut,*) '      number of structures = ', nstru

        ! copy coordinates to work buffer (important: don't modify original)
        coord_work(:,:) = trajectory%coord(:,:)

        ! perform fitting
        select case(fitting_method)
        case(FittingMethodNO)
          ! no fitting

        case(FittingMethodTR_ROT)
          call fit_trrot(n_fitting, fitting_idx, ref_coord, mass_fitting, &
                         coord_work, rot_matrix, com_ref, com_mov, fit_rmsd, ierr)
          call transform(n_atoms, rot_matrix, com_ref, com_mov, coord_work)

        case(FittingMethodTR)
          call fit_trans(n_fitting, fitting_idx, ref_coord, coord_work, mass_fitting, &
                         .true., rot_matrix, com_ref, com_mov, fit_rmsd, ierr)
          call transform(n_atoms, rot_matrix, com_ref, com_mov, coord_work)

        case(FittingMethodXYTR)
          call fit_trans(n_fitting, fitting_idx, ref_coord, coord_work, mass_fitting, &
                         .false., rot_matrix, com_ref, com_mov, fit_rmsd, ierr)
          call transform(n_atoms, rot_matrix, com_ref, com_mov, coord_work)

        case default
          ! For other methods (TR_ZROT, XYTR_ZROT), fall back to TR+ROT
          call fit_trrot(n_fitting, fitting_idx, ref_coord, mass_fitting, &
                         coord_work, rot_matrix, com_ref, com_mov, fit_rmsd, ierr)
          call transform(n_atoms, rot_matrix, com_ref, com_mov, coord_work)

        end select

        ! compute RMSD for analysis atoms using fitted coordinates
        rmsd = 0.0_wp
        tot_mass = 0.0_wp

        do iatom = 1, n_analysis
          idx = analysis_idx(iatom)

          if (mass_weighted) then
            weight = mass(idx)
          else
            weight = 1.0_wp
          end if

          ! For XYTR_ZROT, only consider XY dimensions
          if (fitting_method == FittingMethodXYTR_ZROT) then
            rmsd = rmsd + weight * ( &
              (ref_coord(1, idx) - coord_work(1, idx))**2 + &
              (ref_coord(2, idx) - coord_work(2, idx))**2)
          else
            rmsd = rmsd + weight * ( &
              (ref_coord(1, idx) - coord_work(1, idx))**2 + &
              (ref_coord(2, idx) - coord_work(2, idx))**2 + &
              (ref_coord(3, idx) - coord_work(3, idx))**2)
          end if

          tot_mass = tot_mass + weight
        end do

        if (tot_mass > EPS) rmsd = sqrt(rmsd / tot_mass)

        rmsd_results(nstru) = rmsd

        ! output results
        write(MsgOut,'(a,f10.5)') '              RMSD of analysis atoms = ', rmsd
        write(MsgOut,*) ''

      end if

    end do

    nstru_out = nstru

    ! cleanup
    deallocate(coord_work)
    deallocate(mass_fitting)

    ! Output summary
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze_with_fitting> RMSD analysis completed (with fitting)'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze_with_fitting

end module rmsd_impl_mod
