!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  dr_main
!! @brief   analysis of Drms
!! @authors Chigusa Kobayashi (CK), Daisuke Matsuoka (DM), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module drms_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use drms_impl_mod

  use dr_control_mod
  use dr_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use error_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: dr_analysis_c
  public :: deallocate_drms_results_c

  ! Module-level pointer for results (to be deallocated later)
  real(wp), pointer, save :: dr_ptr(:) => null()

contains
  subroutine dr_analysis_c(molecule, s_trajes_c, ana_period, &
                           ctrl_text, ctrl_len, &
                           result_dr, status, msg, msglen) &
        bind(C, name="dr_analysis_c")
    use conv_f_c_util
    implicit none
    type(s_molecule_c), intent(in) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer, intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out) :: result_dr
    integer(c_int),          intent(out) :: status
    character(kind=c_char),  intent(out) :: msg(*)
    integer(c_int),          value       :: msglen

    type(s_molecule) :: f_molecule

    type(s_error) :: err

    ! Deallocate previous results if any
    if (associated(dr_ptr)) then
      deallocate(dr_ptr)
      nullify(dr_ptr)
    end if

    call error_init(err)
    call c2f_s_molecule(molecule, f_molecule)
    call dr_analysis_main( &
        f_molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, dr_ptr, err)
    if (error_has(err)) then
      call error_to_c(err, status, msg, msglen)
      return
    end if
    status = 0
    if (msglen > 0) msg(1) = c_null_char

    result_dr = c_loc(dr_ptr)
  end subroutine dr_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_drms_results_c
  !> @brief        Deallocate DRMS analysis results
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine deallocate_drms_results_c() bind(C, name="deallocate_drms_results_c")
    implicit none

    if (associated(dr_ptr)) then
      deallocate(dr_ptr)
      nullify(dr_ptr)
    end if

  end subroutine deallocate_drms_results_c

  subroutine dr_analysis_main( &
          molecule, s_trajes_c, ana_period, ctrl_text, ctrl_len, dr, err)
    use, intrinsic :: iso_c_binding
    implicit none
    type(s_molecule), intent(inout) :: molecule
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer,                intent(in) :: ana_period
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer,                intent(in) :: ctrl_len
    real(wp), pointer, intent(out) :: dr(:)
    type(s_error),                   intent(inout) :: err

    ! local variables
    type(s_ctrl_data)      :: ctrl_data
    type(s_trajectory)     :: trajectory
    type(s_output)         :: output
    type(s_option)         :: option


    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.


    ! [Step1] Read control parameters from string
    !
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '

    call control_from_string(ctrl_text, ctrl_len, ctrl_data)


    ! [Step2] Set relevant variables and structures 
    !
    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '

    call setup(molecule, ctrl_data, output, option, err)
    if (error_has(err)) return


    ! [Step3] Analyze trajectory
    !
    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '

    call analyze(molecule, s_trajes_c, ana_period, output, option, &
                 dr, err)
    if (error_has(err)) return


    ! [Step4] Deallocate memory
    !
    write(MsgOut,'(A)') '[STEP4] Deallocate memory'
    write(MsgOut,'(A)') ' '

    call dealloc_trajectory(trajectory)
    call dealloc_molecules_all(molecule)
  end subroutine dr_analysis_main

 !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in DRMS_ANALYSIS
  !! @authors      NT
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(molecule, ctrl_data, output, option, err)
    use dr_control_mod
    use dr_option_mod
    use dr_option_str_mod
    use trajectory_mod
    use output_mod
    use input_mod
    use trajectory_str_mod
    use output_str_mod
    use input_str_mod
    use select_mod
    use molecules_mod
    use molecules_str_mod
    use fileio_prmtop_mod
    use fileio_ambcrd_mod
    use fileio_grotop_mod
    use fileio_grocrd_mod
    use fileio_psf_mod
    use fileio_pdb_mod
    use constants_mod
    implicit none

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_output),          intent(inout) :: output
    type(s_option),          intent(inout) :: option
    type(s_error),           intent(inout) :: err


    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)

    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)

    ! setup option
    !
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, option)

    ! setup contact list
    !
    if (option%two_states) then
      call setup_contact_list_two_states(molecule, option, err)
    else
      call setup_contact_list(molecule, option, err)
    end if

    return

  end subroutine setup

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_contact_list
  !> @brief        setup contact list
  !! @authors      CK
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_contact_list(molecule, option, err)

    ! formal arguments
    type(s_molecule), intent(in)    :: molecule
    type(s_option),   intent(inout) :: option
    type(s_error),    intent(inout) :: err

    logical                         :: duplicate
    integer                         :: i, j, k
    integer                         :: imax, imin
    integer                         :: jmax, jmin
    integer                         :: iatm, jatm
    integer                         :: ires, jres
    integer                         :: iexcl, jexcl
    integer                         :: ncount, ncount_d
    real(wp)                        :: di(1:3)
    real(wp)                        :: dj(1:3)
    real(wp)                        :: d(1:3), dr
    real(wp)                        :: box_ref(1:3)

    character(20)                   :: atom_id
    character(4)                    :: seg_i, seg_j
    integer                         :: imol, jmol
    integer                         :: ifst, icount, itmp
    integer                         :: igrp1, igrp2

    integer                         :: alloc_stat, dealloc_stat
    logical                         :: segflag, molflag, contflag
    logical                         :: exclflag

    integer, allocatable            :: temp_contact_list(:,:)
    integer, allocatable            :: temp_conv(:)
    real(wp), allocatable           :: temp_contact_dist(:)

    molflag=(molecule%molecule_no(1) > 0)
    segflag=(len_trim(molecule%segment_name(1)) > 0)

    box_ref(1:3) = option%box_size_ref(1:3)

    if (option%identical_group) then
      igrp1 = 1
      igrp2 = 1
    else
      igrp1 = 1
      igrp2 = 2
    end if
    dealloc_stat = 0
    alloc_stat = 0

    do icount = 1, 2
      ncount = 0
      do i = 1, option%num_atoms_group(igrp1)
        iatm = option%contact_atoms(i,igrp1)
        ires = molecule%residue_no(iatm)
        seg_i = molecule%segment_name(iatm)
        imol  = molecule%molecule_no(iatm)
        iexcl = option%exclude_group_list(iatm)
        di(1:3) = molecule%atom_refcoord(1:3,iatm)
        ifst  = 1
        if (option%identical_group) ifst = i+1
        do j = ifst,option%num_atoms_group(igrp2)
          jatm = option%contact_atoms(j,igrp2)
          jres = molecule%residue_no(jatm)
          seg_j = molecule%segment_name(jatm)
          jmol  = molecule%molecule_no(jatm)
          jexcl = option%exclude_group_list(jatm)
          exclflag= (iexcl > 0 .and. iexcl .eq. jexcl)
          if (exclflag) cycle
          if (iatm == jatm) cycle

          contflag=.true.
          if (segflag) contflag=(seg_i .eq. seg_j)
          if (molflag) contflag=(imol .eq. jmol)
          if (contflag .and. abs(ires-jres) < option%exclude_residues) cycle
          dj(1:3) = molecule%atom_refcoord(1:3,jatm)
          d(1:3) = di(1:3) - dj(1:3)
          if (option%pbc_correct_setup) &
            d(1:3) = d(1:3)-anint(d(1:3)/box_ref(1:3))*box_ref(1:3)
          dr = sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))
          if (dr >= option%minimum_distance .and.  &
              dr < option%maximum_distance) then
            ncount = ncount + 1
            if (icount == 2) then
              temp_contact_dist(ncount)   = dr
              temp_contact_list(1,ncount) = min(iatm,jatm)
              temp_contact_list(2,ncount) = max(iatm,jatm)
            end if
          end if
        end do
      end do

      if (icount == 1) then
        if (ncount > 0) then
          allocate(temp_contact_dist(1:ncount),      &
                   temp_contact_list(1:2, 1:ncount), &
                   temp_conv(1:ncount),              &
                   stat = alloc_stat)
          if (alloc_stat /= 0) then
              call error_set(err, ERROR_CODE, & 
                'Setup_Contact_List> allocate error')
              return
          end if
          call alloc_option(option, DA_Contact, ncount)
        else
          call error_set(err, ERROR_CODE, & 
              'Setup_Contact_List> ERROR : no contact is defined.')
          return
        end if
      end if
    end do

    if (option%avoid_bonding ) then
      if (molecule%num_bonds > 0) then
        call calc_avoid_bonding(molecule, ncount, temp_contact_list, ncount_d, temp_conv)
        option%num_contact = ncount_d
        do i = 1, option%num_contact
          itmp = temp_conv(i)
          option%contact_list(1:2,i) = temp_contact_list(1:2,itmp)
          option%contact_dist(i) = temp_contact_dist(itmp)
        end do
      else
        call error_set(err, ERROR_CODE, & 
             'Setup_Contact_List> ERROR : bond/angle/dihedral information is required in avoid_bonding option')
        return
      end if
    else
      option%num_contact = ncount
      option%contact_list(1:2,1:ncount) = temp_contact_list(1:2,1:ncount)
      option%contact_dist(1:ncount)     = temp_contact_dist(1:ncount)
    end if

    write(MsgOut,'(A)') 'Setup_Contact_List> Contact List'
    write(MsgOut,'(A20,I10)') '  # of contacts   = ', option%num_contact
    do i = 1, option%num_contact
      write(MsgOut,'(A8,I0,A2,$)') ' Contact',i,': '

      do k  = 1, 2
        iatm = option%contact_list(k,i)
        if (molecule%segment_name(iatm) .ne. "") then
          write(atom_id,'(A,A,A,A,I0,A,A,A,A)') trim(molecule%segment_name(iatm)),':', &
                                            trim(molecule%residue_name(iatm)),':', &
                                            molecule%residue_no(iatm),':', &
                                            trim(molecule%atom_name(iatm))
        else
          write(atom_id,'(A,A,I0,A,A)') trim(molecule%residue_name(iatm)),':', &
                                            molecule%residue_no(iatm),':', &
                                            trim(molecule%atom_name(iatm))
        end if
        write(MsgOut,'(I0,A3,A,A3,$)') option%contact_list(k,i), ' ( ',trim(atom_id),' ) '
        if (k == 1) then
          write(MsgOut,'(A3,$)') ' - '
        end if

      end do
      write(MsgOut,'(F10.2)') option%contact_dist(i)

    end do

     deallocate(temp_contact_dist, &
                temp_contact_list, &
                temp_conv,         &
              stat = dealloc_stat)
    if (dealloc_stat /= 0) then
        call error_set(err, ERROR_CODE, & 
            'Setup_Contact_List> deallocation error')
        return
    end if

    return

  end subroutine setup_contact_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_contact_list_two_states
  !> @brief        setup contact list considering two states
  !! @authors      CK
  !! @param[in]    molecule : molecule information
  !! @param[inout] option   : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_contact_list_two_states(molecule, option, err)

    ! formal arguments
    type(s_molecule), intent(in)    :: molecule
    type(s_option),   intent(inout) :: option
    type(s_error),    intent(inout) :: err

    logical                         :: duplicate
    integer                         :: i, j, k
    integer                         :: imax, imin
    integer                         :: jmax, jmin
    integer                         :: iatm, jatm
    integer                         :: ires, jres
    integer                         :: iexcl, jexcl
    integer                         :: ncount, ncount_d
    real(wp)                        :: di(1:3)
    real(wp)                        :: dj(1:3)
    real(wp)                        :: d(1:3), dr
    real(wp)                        :: ei(1:3)
    real(wp)                        :: ej(1:3)
    real(wp)                        :: e(1:3), er
    real(wp)                        :: box_cur(1:3)
    real(wp)                        :: box_ref(1:3)

    character(20)                   :: atom_id
    character(4)                    :: seg_i, seg_j
    integer                         :: imol, jmol
    integer                         :: ifst, icount, itmp
    integer                         :: igrp1, igrp2

    integer                         :: alloc_stat, dealloc_stat
    logical                         :: segflag, molflag, contflag
    logical                         :: exclflag

    integer, allocatable            :: temp_contact_list(:,:)
    integer, allocatable            :: temp_conv(:)
    real(wp), allocatable           :: temp_contact_dist(:)
    real(wp), allocatable           :: temp_contact_dist_other(:)

    molflag=(molecule%molecule_no(1) > 0)
    segflag=(len_trim(molecule%segment_name(1)) > 0)

    box_cur(1:3) = option%box_size_cur(1:3)
    box_ref(1:3) = option%box_size_ref(1:3)

    dealloc_stat = 0
    alloc_stat = 0

    if (option%identical_group) then
      igrp1 = 1
      igrp2 = 1
    else
      igrp1 = 1
      igrp2 = 2
    end if

    do icount = 1, 2
      ncount = 0
      do i = 1, option%num_atoms_group(igrp1)
        iatm = option%contact_atoms(i,igrp1)
        ires = molecule%residue_no(iatm)
        seg_i = molecule%segment_name(iatm)
        imol  = molecule%molecule_no(iatm)
        iexcl = option%exclude_group_list(iatm)
        di(1:3) = molecule%atom_refcoord(1:3,iatm)
        ei(1:3) = molecule%atom_coord(1:3,iatm)
        ifst  = 1
        if (option%identical_group) ifst = i+1
        do j = ifst,option%num_atoms_group(igrp2)
          jatm = option%contact_atoms(j,igrp2)
          jres = molecule%residue_no(jatm)
          seg_j = molecule%segment_name(jatm)
          jmol  = molecule%molecule_no(jatm)
          jexcl = option%exclude_group_list(jatm)
          exclflag= (iexcl > 0 .and. iexcl .eq. jexcl)
          if (exclflag) cycle
          if (iatm == jatm) cycle

          contflag=.true.
          if (segflag) contflag=(seg_i .eq. seg_j)
          if (molflag) contflag=(imol .eq. jmol)
          if (contflag .and. abs(ires-jres) < option%exclude_residues) cycle
          dj(1:3) = molecule%atom_refcoord(1:3,jatm)
          ej(1:3) = molecule%atom_coord(1:3,jatm)
          d(1:3) = di(1:3) - dj(1:3)
          e(1:3) = ei(1:3) - ej(1:3)
          if (option%pbc_correct_setup) &
            d(1:3) = d(1:3)-anint(d(1:3)/box_ref(1:3))*box_ref(1:3)
          if (option%pbc_correct_setup) &
            e(1:3) = e(1:3)-anint(e(1:3)/box_cur(1:3))*box_cur(1:3)
          dr = sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))
          er = sqrt(e(1)*e(1)+e(2)*e(2)+e(3)*e(3))
          if (((dr >= option%minimum_distance .and.  &
              dr < option%maximum_distance) .or.     &
               (er >= option%minimum_distance .and.  &
                er < option%maximum_distance)) .and. &
              abs(dr-er) >= option%minimum_difference) then
            ncount = ncount + 1
            if (icount == 2) then
              temp_conv(ncount)               = ncount
              temp_contact_dist(ncount)       = dr
              temp_contact_dist_other(ncount) = er
              temp_contact_list(1,ncount)     = min(iatm,jatm)
              temp_contact_list(2,ncount)     = max(iatm,jatm)
            end if
          end if
        end do
      end do

      if (icount == 1) then
        if (ncount > 0) then
          allocate(temp_contact_dist(1:ncount), &
                   temp_contact_dist_other(1:ncount), &
                   temp_contact_list(1:2, 1:ncount),  &
                   temp_conv(1:ncount),               &
                   stat = alloc_stat)
          if (alloc_stat /= 0) then
            call error_set(err, ERROR_CODE, & 
                 'Setup_Contact_List> allocation error')
            return
          end if
          call alloc_option(option, DA_Contact, ncount)
        else
          call error_set(err, ERROR_CODE, & 
                 'Setup_Contact_List> ERROR : no contact is defined.')
          return
        end if
      end if
    end do

    if (option%avoid_bonding ) then
      if (molecule%num_bonds > 0) then
        call calc_avoid_bonding(molecule, ncount, temp_contact_list, ncount_d, temp_conv)
        option%num_contact = ncount_d
        do i = 1, option%num_contact
          itmp = temp_conv(i)
          option%contact_list(1:2,i) = temp_contact_list(1:2,itmp)
          option%contact_dist(i)     = temp_contact_dist(itmp)
          option%contact_cur_dist(i) = temp_contact_dist_other(itmp)
        end do
      else
          call error_set(err, ERROR_CODE, & 
             'Setup_Contact_List> ERROR : bond/angle/dihedral information is required in avoid_bonding option')
          return
      end if
    else
      option%num_contact = ncount
      option%contact_list(1:2,1:ncount) = temp_contact_list(1:2,1:ncount)
      option%contact_dist(1:ncount)     = temp_contact_dist(1:ncount)
      option%contact_cur_dist(1:ncount) = temp_contact_dist_other(1:ncount)
    end if

    write(MsgOut,'(A)') 'Setup_Contact_List> Contact List'
    write(MsgOut,'(A20,I10)') '  # of contacts   = ', option%num_contact
    do i = 1, option%num_contact
      write(MsgOut,'(A8,I0,A2,$)') ' Contact',i,': '

      do k  = 1, 2
        iatm = option%contact_list(k,i)
        if (molecule%segment_name(iatm) .ne. "") then
          write(atom_id,'(A,A,A,A,I0,A,A,A,A)') trim(molecule%segment_name(iatm)),':', &
                                            trim(molecule%residue_name(iatm)),':', &
                                            molecule%residue_no(iatm),':', &
                                            trim(molecule%atom_name(iatm))
        else
          write(atom_id,'(A,A,I0,A,A)') trim(molecule%residue_name(iatm)),':', &
                                            molecule%residue_no(iatm),':', &
                                            trim(molecule%atom_name(iatm))
        end if
        write(MsgOut,'(I0,A3,A,A3,$)') option%contact_list(k,i), ' ( ',trim(atom_id),' ) '
        if (k == 1) then
          write(MsgOut,'(A3,$)') ' - '
        end if

      end do
      write(MsgOut,'(2F10.2)') option%contact_dist(i), option%contact_cur_dist(i)

    end do

     deallocate(temp_contact_dist,       &
                temp_contact_list,       &
                temp_contact_dist_other, &
                temp_conv,               &
              stat = dealloc_stat)
     if (dealloc_stat /= 0) then
       call error_set(err, ERROR_CODE, & 
            'Setup_Contact_List_Two_states> deallocation error')
       return
     end if

    return

  end subroutine setup_contact_list_two_states

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calc_avoid_bonding
  !> @brief        calc_avoid_bonding
  !! @authors      CK
  !! @param[in]    molecule          : molecule information
  !! @param[in]    ncount            : number of members in temporary contact lists
  !! @param[in]    temp_contact_list : temporary contact lists
  !! @param[out]   ncount_d          : number of member in real contact lists
  !! @param[out]   temp_conv         : temporary convert lists
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calc_avoid_bonding(molecule, ncount, temp_contact_list, ncount_d, temp_conv)

    ! formal arguments
    type(s_molecule), intent(in)    :: molecule
    integer,          intent(in)    :: ncount
    integer,          intent(in)    :: temp_contact_list(:,:)
    integer,          intent(out)   :: ncount_d
    integer,          intent(out)   :: temp_conv(:)

    logical                         :: duplicate
    integer                         :: i, j, k
    integer                         :: imax, imin
    integer                         :: jmax, jmin
    integer                         :: iatm, jatm

!
! skip contacts with bonding terms
!
    ncount_d = 0
    do i = 1, ncount
      imin = temp_contact_list(1,i)
      imax = temp_contact_list(2,i)
      duplicate = .false.

      do j = 1, molecule%num_bonds
        jmax = max(molecule%bond_list(1,j), molecule%bond_list(2,j))
        jmin = min(molecule%bond_list(1,j), molecule%bond_list(2,j))
        if (imax == jmax .and. imin == jmin) then
          duplicate = .true.
          exit
        end if
      end do
      if (.not. duplicate) then
        do j = 1, molecule%num_angles
          jmax = max(molecule%angl_list(1,j), molecule%angl_list(3,j))
          jmin = min(molecule%angl_list(1,j), molecule%angl_list(3,j))
          if (imax == jmax .and. imin == jmin) then
            duplicate = .true.
            exit
          end if
        end do
      end if
      if (.not. duplicate) then
        do j = 1, molecule%num_dihedrals
          jmax = max(molecule%dihe_list(1,j), molecule%dihe_list(4,j))
          jmin = min(molecule%dihe_list(1,j), molecule%dihe_list(4,j))
          if (imax == jmax .and. imin == jmin) then
            duplicate = .true.
            exit
          end if
        end do
      end if
      if (.not. duplicate) then
        do j = 1, molecule%num_impropers
          jmax = max(molecule%impr_list(1,j), molecule%impr_list(4,j))
          jmin = min(molecule%impr_list(1,j), molecule%impr_list(4,j))
          if (imax == jmax .and. imin == jmin) then
            duplicate = .true.
            exit
          end if
        end do
      end if
      if (.not. duplicate) then
        ncount_d = ncount_d + 1
        temp_conv(ncount_d) = i
      end if
    end do

    return

  end subroutine calc_avoid_bonding

end module drms_c_mod

