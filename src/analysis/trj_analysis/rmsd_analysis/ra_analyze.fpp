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

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      TM, Claude Code
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, option, trajectory, fitting)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule
    type(s_trj_list), target, intent(in)    :: trj_list
    type(s_output),           intent(in)    :: output
    type(s_option),           intent(inout) :: option
    type(s_trajectory),       intent(inout) :: trajectory
    type(s_fitting),          intent(inout) :: fitting

    ! local variables
    type(s_trj_source)       :: source
    type(s_result_sink)      :: sink
    integer                  :: nstru, n_atoms, n_analysis
    integer                  :: iatom, idx

    real(wp), allocatable    :: mass(:)
    integer,  allocatable    :: analysis_idx(:)


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

    ! Allocate and prepare arrays
    allocate(mass(n_atoms))
    allocate(analysis_idx(n_analysis))

    mass(:) = molecule%mass(:)
    analysis_idx(:) = option%analysis_atom%idx(:)

    ! Initialize source and sink
    call init_source_file(source, trj_list, n_atoms)
    call init_sink_file(sink, output%rmsfile)

    ! Run unified RMSD analysis
    call analyze_rmsd_unified(source, sink, molecule, fitting, option, &
                              mass, analysis_idx, n_analysis, nstru)

    ! Cleanup
    call finalize_sink(sink)
    call finalize_source(source)
    deallocate(mass, analysis_idx)

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
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_rmsd_unified(source, sink, molecule, fitting, option, &
                                  mass, analysis_idx, n_analysis, nstru_out)

    ! formal arguments
    type(s_trj_source),  intent(inout) :: source
    type(s_result_sink), intent(inout) :: sink
    type(s_molecule),    intent(inout) :: molecule
    type(s_fitting),     intent(inout) :: fitting
    type(s_option),      intent(in)    :: option
    real(wp),            intent(in)    :: mass(:)
    integer,             intent(in)    :: analysis_idx(:)
    integer,             intent(in)    :: n_analysis
    integer,             intent(out)   :: nstru_out

    ! local variables
    type(s_trajectory)    :: trajectory
    real(wp), allocatable :: mass_fitting(:)
    real(wp)              :: rmsd, tot_mass, weight
    integer               :: nstru, status, n_atoms
    integer               :: iatom, idx

    n_atoms = size(mass)

    ! Allocate work buffers
    allocate(mass_fitting(n_atoms))

    ! Prepare mass array for fitting
    if (fitting%mass_weight) then
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

      ! Apply fitting if requested
      if (fitting%fitting_method /= FittingMethodNO) then
        call run_fitting(fitting, molecule%atom_coord, trajectory%coord, &
                         trajectory%coord, mass_fitting)
      end if

      ! Compute RMSD for analysis atoms
      rmsd = 0.0_wp
      tot_mass = 0.0_wp

      do iatom = 1, n_analysis
        idx = analysis_idx(iatom)

        if (option%mass_weight) then
          weight = mass(idx)
        else
          weight = 1.0_wp
        end if

        rmsd = rmsd + weight * ( &
               (molecule%atom_coord(1,idx) - trajectory%coord(1,idx))**2 + &
               (molecule%atom_coord(2,idx) - trajectory%coord(2,idx))**2 + &
               (molecule%atom_coord(3,idx) - trajectory%coord(3,idx))**2)

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
    if (allocated(trajectory%coord)) deallocate(trajectory%coord)

    return

  end subroutine analyze_rmsd_unified

end module ra_analyze_mod
