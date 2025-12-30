!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   dr_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Chigusa Kobayashi (CK), Daisuke Matsuoka (DM), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module drms_impl_mod

  use fileio_trj_mod
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
  !> @brief        DRMS analysis (zerocopy, pre-allocated result array)
  !! @authors      Claude Code
  !! @param[in]    contact_list   : contact atom pairs (2, n_contact)
  !! @param[in]    contact_dist   : reference distances for contacts
  !! @param[in]    n_contact      : number of contacts
  !! @param[in]    trajes_c       : trajectory C structure
  !! @param[in]    ana_period     : analysis period
  !! @param[in]    pbc_correct    : apply PBC correction
  !! @param[inout] dr_results     : pre-allocated DRMS results array
  !! @param[out]   nstru_out      : actual number of structures analyzed
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(contact_list, contact_dist, n_contact, &
                     trajes_c, ana_period, pbc_correct, &
                     dr_results, nstru_out)
    use s_trajectories_c_mod

    ! formal arguments
    integer,                 intent(in)    :: contact_list(:,:)
    real(wp),                intent(in)    :: contact_dist(:)
    integer,                 intent(in)    :: n_contact
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period
    logical,                 intent(in)    :: pbc_correct
    real(wp),                intent(inout) :: dr_results(:)
    integer,                 intent(out)   :: nstru_out

    ! local variables
    type(s_trajectory) :: trajectory
    integer            :: nstru, istep
    integer            :: i, i1, i2
    real(wp)           :: d12(3), r12, tmp, dn, drms
    real(wp)           :: box(3)

    ! analysis loop (NO allocation - results array is pre-allocated)
    nstru = 0
    dn = real(n_contact, wp)

    do istep = 1, trajes_c%nframe

      ! read trajectory frame
      call get_frame(trajes_c, istep, trajectory)

      if (mod(istep, ana_period) == 0) then

        nstru = nstru + 1
        write(MsgOut,*) '      number of structures = ', nstru

        ! get box size for PBC
        box(1) = trajectory%pbc_box(1,1)
        box(2) = trajectory%pbc_box(2,2)
        box(3) = trajectory%pbc_box(3,3)

        ! compute DRMS
        drms = 0.0_wp

        !$omp parallel do default(none)                          &
        !$omp private(i, d12, r12, tmp, i1, i2)                  &
        !$omp shared(trajectory, contact_list, contact_dist, n_contact, box, pbc_correct) &
        !$omp reduction(+:drms)
        do i = 1, n_contact
          i1 = contact_list(1, i)
          i2 = contact_list(2, i)
          d12(1:3) = trajectory%coord(1:3, i1) - trajectory%coord(1:3, i2)
          if (pbc_correct) &
            d12(1:3) = d12(1:3) - anint(d12(1:3) / box(1:3)) * box(1:3)
          r12 = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
          tmp = r12 - contact_dist(i)
          drms = drms + tmp * tmp
        end do
        !$omp end parallel do

        if (dn > EPS) drms = sqrt(drms / dn)

        dr_results(nstru) = drms

        ! output results
        write(MsgOut,'(a,f10.5)') '              DRMS = ', drms
        write(MsgOut,*) ''

      end if

    end do

    nstru_out = nstru

    ! Output summary
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> DRMS analysis completed'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze

end module drms_impl_mod
