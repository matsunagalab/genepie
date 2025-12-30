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

  use ra_option_str_mod
  use fileio_trj_mod
  use fitting_mod
  use fitting_str_mod
  use measure_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use select_atoms_mod
  use fileio_mod
  use error_mod
  use messages_mod
  use constants_mod
 
  implicit none
  private

  ! subroutines
  public  :: analyze
  public  :: analyze_zerocopy

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      TM
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trajes_c, ana_period, output, option, &
                     fitting, ra1, err)
    use s_trajectories_c_mod

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(inout) :: option
    type(s_fitting),         intent(inout) :: fitting
    real(wp), pointer,       intent(out)   :: ra1(:)
    type(s_error),                   intent(inout) :: err


    ! local variables
    type(s_trajectory) :: trajectory
    integer                  :: nstru, istep
    integer                  :: iatom, rms_out, idx, i
    real(wp)                 :: rmsd, tot_mass

    real(wp), allocatable    :: mass_fitting(:), mass_analysis(:)


    if (option%check_only) &
      return

    if (error_has(err)) return

!    if (fitting%mass_weight) &
!      call error_msg('Analyze> mass weighted is not allowed')

    if (fitting%mass_weight .and. .not. option%mass_weight) then
      write(MsgOut, *) 'Warning: mass-weighted fitting is enable while RMSD is not mass-weighted '
    else if (.not. fitting%mass_weight .and.  option%mass_weight) then
      write(MsgOut, *) 'Warning: mass-weighted RMSD is enable while fitting is not mass-weighted '
    endif

    allocate(mass_fitting(1:size(molecule%mass(:))), &
             mass_analysis(1:size(option%analysis_atom%idx)))

    if (fitting%mass_weight) then
      if (abs(molecule%mass(1)) < 1.0e-05_wp) then
         call error_set(err, ERROR_MASS_UNDEFINED, &
                        'Analyze> mass is not defined')
         return
      else
         mass_fitting(:) = molecule%mass(:)
      endif
    else
      mass_fitting(:)=1.0_wp
    endif

    if (option%mass_weight) then
      tot_mass = 0.0_wp
      do iatom = 1, size(option%analysis_atom%idx)
        idx = option%analysis_atom%idx(iatom)
        mass_analysis(iatom) = molecule%mass(idx)
        tot_mass = tot_mass + mass_analysis(iatom)
      end do
    else
      mass_analysis(:)=1.0_wp
      tot_mass = real(size(option%analysis_atom%idx),wp)
    endif
    

    ! open output file
    !
    !if (output%rmsfile /= '') &
    !  call open_file(rms_out, output%rmsfile, IOFileOutputNew)

    allocate( ra1(trajes_c%nframe / ana_period) )


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


        ! fitting
        !
        call run_fitting(fitting, &
                         molecule%atom_coord, &
                         trajectory%coord, &
                         trajectory%coord, &
                         mass_fitting)


        ! compute RMSD for selected atoms
        !
        rmsd = 0.0_wp
        do iatom = 1, size(option%analysis_atom%idx)

          idx = option%analysis_atom%idx(iatom)
          if (fitting%fitting_method == FittingMethodXYTR_ZROT) then
            rmsd = rmsd + mass_analysis(iatom) * (                                    &
              + (molecule%atom_coord(1,idx) - trajectory%coord(1,idx))**2  &
              + (molecule%atom_coord(2,idx) - trajectory%coord(2,idx))**2)
          else
            rmsd = rmsd + mass_analysis(iatom) * (                                    &
            + (molecule%atom_coord(1,idx) - trajectory%coord(1,idx))**2  &
            + (molecule%atom_coord(2,idx) - trajectory%coord(2,idx))**2  &
            + (molecule%atom_coord(3,idx) - trajectory%coord(3,idx))**2)
          endif

        end do
        rmsd = sqrt(rmsd/tot_mass)

        ra1(nstru) = rmsd

        ! output results
        !
        write(MsgOut,'(a,f10.5)') '              RMSD of analysis atoms = ',rmsd
        write(MsgOut,*) ''

      end if

    end do

    ! close output file
    !
    !if (output%rmsfile /= '') call close_file(rms_out)


    ! Output summary
    !
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
  !  Subroutine    analyze_zerocopy
  !> @brief        RMSD analysis with true zero-copy (arrays from Python)
  !! @authors      Claude Code
  !! @param[in]    mass           : mass array pointer (view of Python NumPy)
  !! @param[in]    ref_coord      : reference coordinates (3, n_atoms)
  !! @param[in]    trajes_c       : trajectory C structure
  !! @param[in]    ana_period     : analysis period
  !! @param[in]    analysis_idx   : analysis atom indices (1-indexed)
  !! @param[in]    n_analysis     : number of analysis atoms
  !! @param[in]    mass_weighted  : use mass weighting
  !! @param[out]   rmsd_results   : RMSD results
  !! @note         This version does NOT perform fitting. Use when coordinates
  !!               are already aligned or fitting is done in Python.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_zerocopy(mass, ref_coord, trajes_c, ana_period, &
                              analysis_idx, n_analysis, mass_weighted, &
                              rmsd_results)
    use s_trajectories_c_mod

    ! formal arguments
    real(wp), pointer,       intent(in)    :: mass(:)
    real(wp), pointer,       intent(in)    :: ref_coord(:,:)
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period
    integer,                 intent(in)    :: analysis_idx(:)
    integer,                 intent(in)    :: n_analysis
    logical,                 intent(in)    :: mass_weighted
    real(wp), pointer,       intent(out)   :: rmsd_results(:)

    ! local variables
    type(s_trajectory) :: trajectory
    integer            :: nstru, istep
    integer            :: iatom, idx
    real(wp)           :: rmsd, tot_mass, weight

    ! allocate results array
    allocate(rmsd_results(trajes_c%nframe / ana_period))

    ! analysis loop
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

    ! Output summary
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze_zerocopy> RMSD analysis completed (zero-copy, no fitting)'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze_zerocopy

end module rmsd_impl_mod
