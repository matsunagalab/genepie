!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rg_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Motoshi Kamiya (MK), Takaharu Mori (TM), Claude Code
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rg_analyze_mod

  use rg_option_str_mod
  use trj_source_mod
  use result_sink_mod
  use fileio_trj_mod
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
  public  :: analyze_rg_unified

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      MK, TM, Claude Code
  !! @param[in]    molecule   : molecule information
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, option)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule
    type(s_trj_list), target, intent(in)    :: trj_list
    type(s_output),           intent(in)    :: output
    type(s_option),           intent(inout) :: option
    ! local variables
    type(s_trj_source)        :: source
    type(s_result_sink)       :: sink
    integer                   :: nstru, n_atoms, n_analysis
    integer,  allocatable     :: analysis_idx(:)


    if (option%check_only) &
      return

    ! Get atom counts
    n_atoms = size(molecule%mass)
    n_analysis = size(option%analysis_atom%idx)

    ! Get analysis atom indices
    allocate(analysis_idx(n_analysis))
    analysis_idx(:) = option%analysis_atom%idx(:)

    ! Initialize source and sink
    call init_source_file(source, trj_list, n_atoms)
    call init_sink_file(sink, output%rgfile)

    ! Run unified RG analysis with primitive arguments
    call analyze_rg_unified(source, sink, &
                            molecule%mass, n_atoms, &
                            analysis_idx, n_analysis, &
                            option%mass_weighted, nstru)

    ! Cleanup
    call finalize_sink(sink)
    call finalize_source(source)
    deallocate(analysis_idx)

    ! Output summary
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [rgfile] ' // trim(output%rgfile)
    write(MsgOut,'(A)') '    Column 1: Snapshot index'
    write(MsgOut,'(A)') '    Column 2: Radius of gyration (angstrom)'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_rg_unified
  !> @brief        Unified RG analysis loop using source/sink abstractions
  !! @authors      Claude Code
  !! @param[inout] source          : trajectory source (file, memory, or lazy)
  !! @param[inout] sink            : result sink (file or array)
  !! @param[in]    mass            : mass array (n_atoms)
  !! @param[in]    n_atoms         : number of atoms
  !! @param[in]    analysis_idx    : analysis atom indices (1-indexed)
  !! @param[in]    n_analysis      : number of analysis atoms
  !! @param[in]    mass_weighted   : use mass weighting
  !! @param[out]   nstru_out       : number of frames analyzed
  !! @note         This function can be called from both CLI and Python interface.
  !!               The source abstraction handles FILE, MEMORY, and LAZY_DCD modes.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_rg_unified(source, sink, mass, n_atoms, &
                                analysis_idx, n_analysis, &
                                mass_weighted, nstru_out)

    ! formal arguments
    type(s_trj_source),  intent(inout) :: source
    type(s_result_sink), intent(inout) :: sink
    real(wp),            intent(in)    :: mass(:)         ! (n_atoms)
    integer,             intent(in)    :: n_atoms
    integer,             intent(in)    :: analysis_idx(:) ! analysis atom indices
    integer,             intent(in)    :: n_analysis
    logical,             intent(in)    :: mass_weighted
    integer,             intent(out)   :: nstru_out

    ! local variables
    type(s_trajectory)    :: trajectory
    real(wp)              :: com(3), weight, tot_weight, rg
    integer               :: nstru, status
    integer               :: i, iatom, idx

    ! Main analysis loop
    nstru = 0

    do while (has_more_frames(source))

      ! Get next frame
      call get_next_frame(source, trajectory, status)
      if (status /= 0) exit

      nstru = nstru + 1

      write(MsgOut,*) '      number of structures = ', nstru

      ! Compute center of mass
      com(1:3) = 0.0_wp
      tot_weight = 0.0_wp
      do iatom = 1, n_analysis
        idx = analysis_idx(iatom)
        weight = 1.0_wp
        if (mass_weighted) weight = mass(idx)
        com(:) = com(:) + weight * trajectory%coord(:,idx)
        tot_weight = tot_weight + weight
      end do
      com(:) = com(:) / tot_weight

      ! Compute radius of gyration
      rg = 0.0_wp
      do iatom = 1, n_analysis
        idx = analysis_idx(iatom)
        weight = 1.0_wp
        if (mass_weighted) weight = mass(idx)
        do i = 1, 3
          rg = rg + weight * ((trajectory%coord(i,idx) - com(i)) ** 2)
        end do
      end do
      rg = sqrt(rg / tot_weight)

      ! Write result
      call write_result_with_index(sink, nstru, rg)

      ! Output progress
      write(MsgOut,'(a,f10.5)') '              RG of analysis atoms = ', rg
      write(MsgOut,*) ''

    end do

    nstru_out = nstru

    ! Cleanup
    if (allocated(trajectory%coord)) deallocate(trajectory%coord)

    return

  end subroutine analyze_rg_unified

end module rg_analyze_mod
