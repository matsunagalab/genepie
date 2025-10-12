!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   crd_convert_convert
!> @brief   convert trajectory files
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module crd_convert_convert

  use cc_option_mod
  use cc_option_str_mod
  use pbc_correct_mod
  use fitting_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use select_atoms_mod
  use molecules_mod
  use molecules_str_mod
  use constants_mod
  use measure_mod
  use fileio_trj_mod
  use fileio_mod
  use messages_mod
  use string_mod
  use fileio_pdb_mod
  use s_trajectories_c_mod

  implicit none
  private

  ! subroutines
  public  :: convert
  private :: centering
  private :: output_split_trjpdb
  private :: get_filename

  real(wp), allocatable, save   :: tmp_coord(:,:)

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert
  !> @brief        convert trajectory files
  !! @authors      NT
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !! @param[inout] s_trajs_c_array : output multi flame trajectories
  !! @param[out]   num_trajs       : length of s_trajs_c_array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

 subroutine convert(molecule,   &
                    trj_list,   &
                    trajectory, &
                    fitting,    &
                    option,     &
                    output,     &
                    s_trajs_c_array, &
                    num_trajs)
    use, intrinsic :: iso_c_binding

    ! formal arguments
    type(s_molecule),   intent(inout) :: molecule
    type(s_trj_list),   intent(inout) :: trj_list
    type(s_trajectory), intent(inout) :: trajectory
    type(s_fitting),    intent(inout) :: fitting
    type(s_option),     intent(inout) :: option
    type(s_output),     intent(inout) :: output
    type(c_ptr), intent(inout) :: s_trajs_c_array
    integer(c_int), intent(out) :: num_trajs

    ! local variables
    type(s_trj_file)         :: trj_in, trj_out
    integer                  :: nstru, irun, itrj
    integer                  :: rms_out, trr_out
    type(s_trajectories_c), pointer :: s_trajs_c_buf(:)


    ! check-only
    if (option%check_only) &
      return

    ! open output files
    ! if (output%trjfile /= '' .and. .not. option%split_trjpdb) &
    !   call open_trj (trj_out,              &
    !                  output%trjfile,       &
    !                  option%trjout_format, &
    !                  option%trjout_type,   &
    !                  IOFileOutputNew)

    if (output%rmsfile /= '') &
      call open_file(rms_out, output%rmsfile, IOFileOutputNew)

    if (output%trrfile /= '') &
      call open_file(trr_out, output%trrfile, IOFileOutputNew)

    nstru = 0

    ! init s_trajs_c_array
    num_trajs = size(trj_list%md_steps)
    allocate(s_trajs_c_buf(num_trajs))

    do irun = 1, size(trj_list%md_steps)

      call open_trj(trj_in,                   &
                    trj_list%filenames(irun), &
                    trj_list%trj_format,      &
                    trj_list%trj_type,        &
                    IOFileInput)

      do itrj = 1, trj_list%md_steps(irun)

        ! input trj
        !
        call read_trj(trj_in, trajectory)

        if (itrj == 1) then
            ! Initialize trajectory with selected atom count
            if (allocated(option%trjout_atom%idx)) then
              call init_empty_s_trajectories_c( &
                  s_trajs_c_buf(irun), size(option%trjout_atom%idx), sum(trj_list%md_steps))
            else
              call init_empty_s_trajectories_c( &
                  s_trajs_c_buf(irun), size(trajectory%coord, dim=2), sum(trj_list%md_steps))
            end if
        end if


        if (mod(itrj, trj_list%ana_periods(irun)) == 0) then

          nstru = nstru + 1
          write(MsgOut,*) '      number of structures = ', nstru

          ! selection
          !
          call reselect_atom(molecule, &
                             option%trjout_atom_exp, &
                             trajectory%coord, &
                             option%trjout_atom, &
                             option%trjout_atom_trj)

          ! centering
          !
          call centering(molecule, trajectory%coord, option)

          ! pbc-correct
          !
          call run_pbc_correct(option%pbcc_mode, &
                               molecule, &
                               trajectory)

          ! fitting
          !
          if (fitting%mass_weight) then
            call run_fitting(fitting, &
                             molecule%atom_coord, &
                             trajectory%coord, &
                             trajectory%coord, &
                             molecule%mass)

          else
            call run_fitting(fitting, &
                             molecule%atom_coord, &
                             trajectory%coord, &
                             trajectory%coord)

          end if

          ! write data
          !
          if (output%rmsfile /= '') &
            call out_rmsd (rms_out, nstru, fitting)

          if (output%trrfile /= '') &
            call out_trrot(trr_out, nstru, fitting)

          if (output%trjfile /= '') then
            if (option%split_trjpdb) then
              call output_split_trjpdb(nstru, molecule, trajectory, option, output)
            else
              ! call write_trj(trj_out, trajectory, option%trjout_atom_trj,molecule)
            end if
          end if

          ! Filter trajectory coordinates based on selected atoms
          if (allocated(option%trjout_atom%idx)) then
            call set_frame_filtered(s_trajs_c_buf(irun), trajectory, itrj, option%trjout_atom)
          else
            call set_frame(s_trajs_c_buf(irun), trajectory, itrj)
          end if
        end if

      end do

      call close_trj(trj_in)

    end do

    if (output%trrfile /= '') call close_file(trr_out)
    if (output%rmsfile /= '') call close_file(rms_out)
    ! if (output%trjfile /= '' .and. .not. option%split_trjpdb) &
    !                           call close_trj (trj_out)

    s_trajs_c_array = c_loc(s_trajs_c_buf)
    return

  end subroutine convert

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    centering
  !> @brief        move the COM of the fitting target to the origin
  !! @authors      DM
  !! @param[in]    molecule   : molecule information
  !! @param[inout] coord      : atom coordinates
  !! @param[in]    fitting    : fitting information
  !! @param[in]    option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine centering(molecule, coord, option)

  ! formal arguments
  type(s_molecule), intent(in)    :: molecule
  real(wp),         intent(inout) :: coord(:,:)
  type(s_option),   intent(in)    :: option

  ! local variables
  integer  :: iatm, natm
  real(wp) :: com(3)

  if (.not. option%centering) return

  natm = size(molecule%atom_no)
  com  = compute_com(coord, molecule%mass, option%centering_atom%idx)

  do iatm = 1, natm
    coord(:,iatm) = coord(:,iatm) - com(:) + option%center_coord(:)
  end do

  return

  end subroutine centering

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_split_trjpdb
  !> @brief        output split PDB files as trajectory
  !! @authors      TM
  !! @param[inout] nstru      : structure index
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_split_trjpdb(nstru, molecule, trajectory, option, output)

    ! formal arguments
    integer,                 intent(inout) :: nstru
    type(s_molecule),        intent(inout) :: molecule
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_pdb)              :: tmp_pdb


    ! allocate tmp_coord to save molecule%atom_coord
    !
    if(.not. allocated(tmp_coord)) then
      allocate(tmp_coord(3,molecule%num_atoms))
    end if

    tmp_coord(:,:) = molecule%atom_coord(:,:)
    molecule%atom_coord(:,:) = trajectory%coord(:,:)

    call export_molecules(molecule, option%trjout_atom, tmp_pdb)

    if (option%trjout_type == TrjTypeCoorBox) then
      tmp_pdb%cryst_rec = .true.
      tmp_pdb%pbc_box(1:3,1:3) = trajectory%pbc_box(1:3,1:3)
    else
      tmp_pdb%cryst_rec = .false.
    end if

    tmp_pdb%model_rec = .false.

    call output_pdb(get_filename(output%trjfile,nstru), tmp_pdb)

    molecule%atom_coord(:,:) = tmp_coord(:,:)

    return

  end subroutine output_split_trjpdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_filename
  !> @brief        insert snapshot index into {} in the filename
  !! @authors      TM
  !! @param[in]    filename      : filename
  !! @param[in]    no            : index
  !! @note         this subroutine was originally made by NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_filename(filename, no)

    ! return
    character(Maxfilename)   :: get_filename

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: no

    ! local variables
    integer                  :: bl, br
    character(100)           :: fid


    bl = index(filename, '{', back=.true.)
    br = index(filename, '}', back=.true.)

    if (bl == 0 .or. br == 0 .or. bl > br) &
      call error_msg('Get_Filename> {} is not found in the output trjfile name')

    write(fid,'(i0)') no
    get_filename = filename(:bl-1) // trim(fid) // filename(br+1:)

    return

  end function get_filename

end module crd_convert_convert
