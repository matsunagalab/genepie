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

  use dr_option_str_mod
  use fileio_trj_mod
  use fitting_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_mod
  use measure_mod
  use messages_mod
  use constants_mod
  use atom_libs_mod
  use error_mod
 
  implicit none
  private

  ! structure
  type, private :: s_contact
    integer               :: n_pair
    real(wp), allocatable :: r0_ij(:)
    integer,  allocatable :: cnt_pair(:,:)
  end type s_contact


  ! subroutines
  public  :: analyze
  public  :: analyze_zerocopy
  private :: compute_drms

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      CK, DM, NT
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !! @note         JPCB (2017) 121, 3364 - 3375
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trajes_c, ana_period, output, option, &
                     dr1, err)
    use s_trajectories_c_mod

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(inout) :: option
    real(wp), pointer,       intent(out)   :: dr1(:)
    type(s_error),           intent(inout) :: err


    ! local variables
    type(s_trajectory) :: trajectory
    integer                  :: nstru, ifile, istep
    integer                  :: rms_out
    real(wp)                 :: drms, drms_cur


    if (option%check_only) &
      return

    ! open output file
    !
    !if (output%rmsfile /= '') &
    !  call open_file(rms_out, output%rmsfile, IOFileOutputNew)

    allocate( dr1(trajes_c%nframe / ana_period) )


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

        if (option%two_states) then
          call compute_drms_two_states(option, trajectory%coord,  &
                                       trajectory%pbc_box, drms, drms_cur, err)
          if (error_has(err)) return

          !if (output%rmsfile /= '') &
          !  write(rms_out, '(i10,1x,2f8.3)') nstru, drms, drms_cur
        else
          call compute_drms(option, trajectory%coord, trajectory%pbc_box, drms, err)
          if (error_has(err)) return

          !if (output%rmsfile /= '') &
          !  write(rms_out, '(i10,1x,f8.3)') nstru, drms
        endif

        dr1(nstru) = drms

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
    write(MsgOut,'(A)') '    Column 2: Distance RMSD (angstrom)'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine  compute_drms
  !> @brief      calculate drms
  !! @authors    CK
  !! @param[in]  enefunc potential energy functions [str]
  !! @param[in]  coord coordinates of target systems [dble]
  !! @date       2018/07/19 (CK)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_drms(option, coord, pbc_box, drms, err)

    ! formal arguments
    type(s_option), target, intent(in)    :: option
    real(wp),       intent(in)    :: coord(:,:)
    real(wp),       intent(in)    :: pbc_box(3,3)
    real(wp),       intent(out)   :: drms
    type(s_error),  intent(inout) :: err

    ! local variables
    integer                       :: i, i1, i2,  num_contact
    real(wp)                      :: d12(1:3), r12
    real(wp)                      :: t, dn, tmp
    real(wp)                      :: box(1:3)
    integer,              pointer :: contact_list(:,:)
    real(wp),             pointer :: contact_dist(:)
    logical                       :: pbc_flag

    contact_list   => option%contact_list
    contact_dist   => option%contact_dist

    pbc_flag = option%pbc_correct

    box(1) = pbc_box(1,1)
    box(2) = pbc_box(2,2)
    box(3) = pbc_box(3,3)
    if (pbc_flag) then
      if (box(1) <= 0.0_wp .or.  &
        box(2) <= 0.0_wp .or.    &
        box(3) <= 0.0_wp  ) then
        call error_set(err, ERROR_PBC_BOX, &
                       'Compute_Drms> Box is required for pbc_correct')
        return
      endif
    endif

    num_contact = option%num_contact
    dn          = real(num_contact,wp)
    drms = 0.0_wp

    !$omp parallel do default(none)                          &
    !$omp private(i,  d12, r12,  tmp, i1, i2)    &
    !$omp shared(coord, contact_list, contact_dist, num_contact, box, pbc_flag)   &
    !$omp reduction(+:drms)
    !

    do i = 1, num_contact

      i1 = contact_list(1,i)
      i2 = contact_list(2,i)
      d12(1:3) = coord(1:3,i1) - coord(1:3,i2)
      if (pbc_flag) &
        d12(1:3) = d12(1:3)-anint(d12(1:3)/box(1:3))*box(1:3)
      r12      = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
        
      tmp      = r12 - contact_dist(i)
      drms     = drms + tmp*tmp

    end do
    !$omp end parallel do

    if (dn > EPS) drms = sqrt(drms/dn)

    return
  end subroutine compute_drms


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine  compute_drms_two_states
  !> @brief      calculate drms
  !! @authors    CK
  !! @param[in]  enefunc potential energy functions [str]
  !! @param[in]  coord coordinates of target systems [dble]
  !! @date       2018/08/16 (CK)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_drms_two_states(option, coord, pbc_box, drms, drms_cur, err)

    ! formal arguments
    type(s_option), target, intent(in)    :: option
    real(wp),       intent(in)    :: coord(:,:)
    real(wp),       intent(in)    :: pbc_box(3,3)
    real(wp),       intent(out)   :: drms
    real(wp),       intent(out)   :: drms_cur
    type(s_error),  intent(inout) :: err

    ! local variables
    integer                       :: i, i1, i2,  num_contact
    real(wp)                      :: d12(1:3), r12
    real(wp)                      :: t, dn, tmp, tmp2
    real(wp)                      :: box(1:3)
    integer,              pointer :: contact_list(:,:)
    real(wp),             pointer :: contact_dist(:)
    real(wp),             pointer :: contact_cur_dist(:)
    logical                       :: pbc_flag

    contact_list     => option%contact_list
    contact_dist     => option%contact_dist
    contact_cur_dist => option%contact_cur_dist

    pbc_flag = option%pbc_correct

    box(1) = pbc_box(1,1)
    box(2) = pbc_box(2,2)
    box(3) = pbc_box(3,3)
    if (pbc_flag) then
      if (box(1) <= 0.0_wp .or.  &
        box(2) <= 0.0_wp .or.    &
        box(3) <= 0.0_wp  ) then
        call error_set(err, ERROR_PBC_BOX, &
                       'Compute_Drms_Two_States> Box is required for pbc_correct')
        return
      endif
    endif

    num_contact = option%num_contact
    dn          = real(num_contact,wp)
    drms        = 0.0_wp
    drms_cur    = 0.0_wp

    !$omp parallel do default(none)                          &
    !$omp private(i,  d12, r12,  tmp, tmp2, i1, i2)    &
    !$omp shared(coord, contact_list, contact_dist, contact_cur_dist, num_contact,  &
    !$omp        box, pbc_flag)   &
    !$omp reduction(+:drms) reduction(+:drms_cur)
    !

    do i = 1, num_contact

      i1 = contact_list(1,i)
      i2 = contact_list(2,i)
      d12(1:3) = coord(1:3,i1) - coord(1:3,i2)
      if (pbc_flag) &
        d12(1:3) = d12(1:3)-anint(d12(1:3)/box(1:3))*box(1:3)
      r12      = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
        
      tmp      = r12 - contact_dist(i)
      tmp2     = r12 - contact_cur_dist(i)
      drms     = drms     + tmp*tmp
      drms_cur = drms_cur + tmp2*tmp2

    end do
    !$omp end parallel do

    if (dn > EPS) drms     = sqrt(drms/dn)
    if (dn > EPS) drms_cur = sqrt(drms_cur/dn)

    return
  end subroutine compute_drms_two_states

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_zerocopy
  !> @brief        DRMS analysis with true zero-copy (contact data from Python)
  !! @authors      Claude Code
  !! @param[in]    contact_list   : contact atom pairs (2, n_contact)
  !! @param[in]    contact_dist   : reference distances for contacts
  !! @param[in]    n_contact      : number of contacts
  !! @param[in]    trajes_c       : trajectory C structure
  !! @param[in]    ana_period     : analysis period
  !! @param[in]    pbc_correct    : apply PBC correction
  !! @param[out]   dr_results     : DRMS results
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_zerocopy(contact_list, contact_dist, n_contact, &
                              trajes_c, ana_period, pbc_correct, &
                              dr_results)
    use s_trajectories_c_mod

    ! formal arguments
    integer,                 intent(in)    :: contact_list(:,:)
    real(wp),                intent(in)    :: contact_dist(:)
    integer,                 intent(in)    :: n_contact
    type(s_trajectories_c),  intent(in)    :: trajes_c
    integer,                 intent(in)    :: ana_period
    logical,                 intent(in)    :: pbc_correct
    real(wp), pointer,       intent(out)   :: dr_results(:)

    ! local variables
    type(s_trajectory) :: trajectory
    integer            :: nstru, istep
    integer            :: i, i1, i2
    real(wp)           :: d12(3), r12, tmp, dn, drms
    real(wp)           :: box(3)

    ! allocate results array
    allocate(dr_results(trajes_c%nframe / ana_period))

    ! analysis loop
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

    ! Output summary
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze_zerocopy> DRMS analysis completed (zero-copy)'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze_zerocopy

end module drms_impl_mod
