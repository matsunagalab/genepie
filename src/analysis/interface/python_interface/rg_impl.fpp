!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rg_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Motoshi Kamiya (MK), Takaharu Mori (TM)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rg_impl_mod

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

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        RG analysis (zerocopy, pre-allocated result array)
  !! @authors      Claude Code
  !! @param[in]    mass           : mass array pointer (view of Python NumPy)
  !! @param[in]    trajes_c       : trajectory C structure
  !! @param[in]    ana_period     : analysis period
  !! @param[in]    analysis_idx   : analysis atom indices (1-indexed)
  !! @param[in]    n_analysis     : number of analysis atoms
  !! @param[in]    mass_weighted  : use mass weighting
  !! @param[inout] rg_results     : pre-allocated result array (zerocopy from Python)
  !! @param[out]   nstru_out      : number of structures analyzed
  !! @note         Result array is passed from Python (no allocation in Fortran)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(mass, trajes_c, ana_period, &
                     analysis_idx, n_analysis, mass_weighted, &
                     rg_results, nstru_out)
    use s_trajectories_c_mod

    ! formal arguments
    real(wp), pointer,       intent(in)    :: mass(:)
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period
    integer,                 intent(in)    :: analysis_idx(:)
    integer,                 intent(in)    :: n_analysis
    logical,                 intent(in)    :: mass_weighted
    real(wp),                intent(inout) :: rg_results(:)  ! pre-allocated
    integer,                 intent(out)   :: nstru_out

    ! local variables
    type(s_trajectory) :: trajectory
    integer            :: nstru, istep
    integer            :: i, iatom, idx
    real(wp)           :: com(3), weight, tot_weight, rg

    ! analysis loop
    nstru = 0

    do istep = 1, trajes_c%nframe

      ! read trajectory frame
      call get_frame(trajes_c, istep, trajectory)

      if (mod(istep, ana_period) == 0) then

        nstru = nstru + 1
        write(MsgOut,*) '      number of structures = ', nstru

        ! compute center of mass
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

        ! compute radius of gyration
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

        ! Write directly to pre-allocated result array (zerocopy)
        rg_results(nstru) = rg

        ! output results
        write(MsgOut,'(a,f10.5)') '              RG of analysis atoms = ',rg
        write(MsgOut,*) ''

      end if

    end do

    nstru_out = nstru

    ! Output summary
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> RG analysis completed'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze

end module rg_impl_mod
