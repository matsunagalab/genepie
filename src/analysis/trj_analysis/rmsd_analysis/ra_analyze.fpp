!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ra_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Takaharu Mori (TM), Claude Code
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ra_analyze_mod

  use ra_option_str_mod
  use trj_source_mod
  use result_sink_mod
  use fileio_trj_mod
  use fitting_mod
  use fitting_str_mod
  use measure_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use select_atoms_mod
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  public  :: analyze_rmsd_unified

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      TM, Claude Code
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, option, fitting)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule
    type(s_trj_list), target, intent(in)    :: trj_list
    type(s_output),           intent(in)    :: output
    type(s_option),           intent(inout) :: option
    type(s_fitting),          intent(inout) :: fitting

    ! local variables
    type(s_trj_source)       :: source
    type(s_result_sink)      :: sink
    integer                  :: nstru, n_atoms, n_analysis, n_fitting
    integer,  allocatable    :: fitting_idx(:)
    integer,  allocatable    :: analysis_idx(:)
    integer                  :: fitting_method


    if (option%check_only) &
      return

    if (fitting%mass_weight .and. .not. option%mass_weight) then
      write(MsgOut, *) 'Warning: mass-weighted fitting is enable while RMSD is not mass-weighted '
    else if (.not. fitting%mass_weight .and.  option%mass_weight) then
      write(MsgOut, *) 'Warning: mass-weighted RMSD is enable while fitting is not mass-weighted '
    endif

    ! Get atom counts
    n_atoms = size(molecule%mass)
    n_analysis = size(option%analysis_atom%idx)

    ! Get fitting atom indices (may be empty if no fitting)
    if (fitting%fitting_method /= FittingMethodNO .and. &
        allocated(fitting%fitting_atom%idx)) then
      n_fitting = size(fitting%fitting_atom%idx)
      allocate(fitting_idx(n_fitting))
      fitting_idx(:) = fitting%fitting_atom%idx(:)
      fitting_method = fitting%fitting_method
    else
      n_fitting = 0
      allocate(fitting_idx(1))  ! Dummy allocation for Fortran
      fitting_idx(1) = 1
      fitting_method = FittingMethodNO
    end if

    ! Get analysis atom indices
    allocate(analysis_idx(n_analysis))
    analysis_idx(:) = option%analysis_atom%idx(:)

    ! Initialize source and sink
    call init_source_file(source, trj_list, n_atoms)
    call init_sink_file(sink, output%rmsfile)

    ! Run unified RMSD analysis with primitive arguments
    call analyze_rmsd_unified(source, sink, &
                              molecule%atom_coord, molecule%mass, n_atoms, &
                              fitting_idx, n_fitting, &
                              analysis_idx, n_analysis, &
                              fitting_method, option%mass_weight, nstru)

    ! Cleanup
    call finalize_sink(sink)
    call finalize_source(source)
    deallocate(fitting_idx, analysis_idx)

    ! Output summary
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [rmsfile] ' // trim(output%rmsfile)
    write(MsgOut,'(A)') '    Column 1: Snapshot index'
    write(MsgOut,'(A)') '    Column 2: Root-mean-square deviation (RMSD)(angstrom)'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_rmsd_unified
  !> @brief        Unified RMSD analysis loop using source/sink abstractions
  !! @authors      Claude Code
  !! @param[inout] source           : trajectory source (file, memory, or lazy)
  !! @param[inout] sink             : result sink (file or array)
  !! @param[in]    ref_coord        : reference coordinates (3, n_atoms)
  !! @param[in]    mass             : mass array (n_atoms)
  !! @param[in]    n_atoms          : number of atoms
  !! @param[in]    fitting_idx      : fitting atom indices (1-indexed)
  !! @param[in]    n_fitting        : number of fitting atoms (0 = no fitting)
  !! @param[in]    analysis_idx     : analysis atom indices (1-indexed)
  !! @param[in]    n_analysis       : number of analysis atoms
  !! @param[in]    fitting_method   : fitting method (FittingMethodNO, etc.)
  !! @param[in]    mass_weighted    : use mass weighting for RMSD
  !! @param[out]   nstru_out        : number of frames analyzed
  !! @note         This function can be called from both CLI and Python interface.
  !!               The source abstraction handles FILE, MEMORY, and LAZY_DCD modes.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_rmsd_unified(source, sink, ref_coord, mass, n_atoms, &
                                  fitting_idx, n_fitting, &
                                  analysis_idx, n_analysis, &
                                  fitting_method, mass_weighted, nstru_out)

    ! formal arguments
    type(s_trj_source),  intent(inout) :: source
    type(s_result_sink), intent(inout) :: sink
    real(wp),            intent(in)    :: ref_coord(:,:)  ! (3, n_atoms)
    real(wp),            intent(in)    :: mass(:)         ! (n_atoms)
    integer,             intent(in)    :: n_atoms
    integer,             intent(in)    :: fitting_idx(:)  ! fitting atom indices
    integer,             intent(in)    :: n_fitting       ! 0 = no fitting
    integer,             intent(in)    :: analysis_idx(:) ! analysis atom indices
    integer,             intent(in)    :: n_analysis
    integer,             intent(in)    :: fitting_method  ! FittingMethodNO, etc.
    logical,             intent(in)    :: mass_weighted   ! RMSD mass weighting
    integer,             intent(out)   :: nstru_out

    ! local variables
    type(s_trajectory)    :: trajectory
    real(wp), allocatable :: mass_fitting(:)
    real(wp), allocatable :: coord_work(:,:)
    real(wp)              :: rmsd, tot_mass, weight
    integer               :: nstru, status
    integer               :: iatom, idx

    ! Fitting-related variables
    real(wp) :: rot_matrix(3,3)
    real(wp) :: com_ref(3), com_mov(3)
    real(wp) :: fit_rmsd
    integer  :: ierr

    ! Allocate work buffers
    allocate(mass_fitting(n_atoms))
    allocate(coord_work(3, n_atoms))

    ! Prepare mass array for fitting (always use mass for fitting if requested)
    if (mass_weighted) then
      mass_fitting(:) = mass(:)
    else
      mass_fitting(:) = 1.0_wp
    end if

    ! Main analysis loop
    nstru = 0

    do while (has_more_frames(source))

      ! Get next frame
      call get_next_frame(source, trajectory, status)
      if (status /= 0) exit

      nstru = nstru + 1

      write(MsgOut,*) '      number of structures = ', nstru

      ! Copy coordinates to work buffer (fitting modifies coordinates)
      coord_work(:,:) = trajectory%coord(:,:)

      ! Apply fitting if requested
      if (n_fitting > 0 .and. fitting_method /= FittingMethodNO) then
        select case(fitting_method)
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
      end if

      ! Compute RMSD for analysis atoms
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
                 (ref_coord(1,idx) - coord_work(1,idx))**2 + &
                 (ref_coord(2,idx) - coord_work(2,idx))**2)
        else
          rmsd = rmsd + weight * ( &
                 (ref_coord(1,idx) - coord_work(1,idx))**2 + &
                 (ref_coord(2,idx) - coord_work(2,idx))**2 + &
                 (ref_coord(3,idx) - coord_work(3,idx))**2)
        end if

        tot_mass = tot_mass + weight
      end do

      if (tot_mass > EPS) then
        rmsd = sqrt(rmsd / tot_mass)
      else
        rmsd = 0.0_wp
      end if

      ! Write result
      call write_result_with_index(sink, nstru, rmsd)

      ! Output progress
      write(MsgOut,'(a,f10.5)') '              RMSD of analysis atoms = ', rmsd
      write(MsgOut,*) ''

    end do

    nstru_out = nstru

    ! Cleanup
    deallocate(mass_fitting)
    deallocate(coord_work)
    if (allocated(trajectory%coord)) deallocate(trajectory%coord)

    return

  end subroutine analyze_rmsd_unified

end module ra_analyze_mod
