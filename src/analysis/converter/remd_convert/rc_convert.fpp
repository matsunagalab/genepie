!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rc_convert_mod
!> @brief   convert remd trajectory files
!! @authors Takaharu Mori (TM), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rc_convert_mod

  use constants_mod
  use rc_option_str_mod
  use pbc_correct_mod
  use fitting_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use input_str_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use, intrinsic :: iso_c_binding
  use fileio_trj_mod
  use fileio_mod
  use string_mod
  use measure_mod
  use messages_mod
 
  implicit none
  private

  ! subroutines
  public  :: convert
  public  :: remd_count_frames
  public  :: remd_fill_frames
  private :: setup_convert
  private :: centering
  private :: get_replicate_name
  private :: transfer_dcd_frame
  private :: expr_has_around

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert
  !> @brief        convert remd trajectory files
  !! @authors      TM
  !! @param[in]    input      : input information
  !! @param[in]    output     : output information
  !! @param[in]    option     : option information
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine convert(input, output, option, molecule, trajectory, fitting)

    ! formal arguments
    type(s_input),           intent(in)    :: input
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(in)    :: option
    type(s_molecule),        intent(inout) :: molecule
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting

    ! local variables
    integer                       :: i, j, k, m, n
    integer                       :: num_conv_id, nlines
    integer                       :: num_open_files
    integer                       :: istep, parmid, unit_no, imdstep
    logical                       :: conv_crd, conv_log, conv_ene, do_skip
    logical                       :: done, do_write
    logical                       :: static_sel, sel_all, fast_dcd
    integer                       :: ierr, natm
    character(MaxFilename)        :: filename
    character(2000)               :: line
    type(s_selatoms)              :: trjout_atom_trj
    integer,          allocatable :: conv_id(:), open_id(:), rem_unit_in(:)
    integer,          allocatable :: log_unit_in(:), log_unit_out(:)
    integer,          allocatable :: ene_unit_in(:), ene_unit_out(:)
    integer,          allocatable :: param_index(:)
    integer,          allocatable :: dest(:)
    logical,          allocatable :: do_open_file(:)
    logical,          allocatable :: in_ready(:), out_ready(:)
    real(sp),         allocatable :: xbuf(:), ybuf(:), zbuf(:)
    type(s_trj_file), allocatable :: trj_in(:), trj_out(:)


    ! check check only
    !
    if (option%check_only) &
      return


    ! check remfile, dcdfile, logfile, enefile existing
    !
    call setup_convert(input, output, option, conv_crd, conv_log, conv_ene)


    ! make convert id lists
    !
    if (.not. allocated(option%convert_ids)) then
      allocate(conv_id(option%num_replicas))
      do i = 1, option%num_replicas
        conv_id(i) = i
      end do
    else
      allocate(conv_id(size(option%convert_ids)))
      conv_id(:) = option%convert_ids(:)
      do i = 1, size(conv_id)
        if (conv_id(i) > option%num_replicas) &
          call error_msg('Convert> Convert IDs is out-of-range.')
      end do
    end if
    num_conv_id = size(conv_id)


    ! other setup
    !
    allocate(open_id     (option%num_replicas), &
             do_open_file(option%num_replicas), &
             trj_in      (option%num_replicas), &
             trj_out     (option%num_replicas), &
             rem_unit_in (option%num_replicas), &
             log_unit_in (option%num_replicas), &
             log_unit_out(option%num_replicas), &
             ene_unit_in (option%num_replicas), &
             ene_unit_out(option%num_replicas), &
             param_index (option%num_replicas))


    ! make lists to be opened
    !
    if (option%convert_type == ConvertTypeParameter) then

      do_open_file(:) = .false.
      do i = 1, option%num_replicas

        nlines = option%nsteps/option%exchange_period
        filename = get_replicate_name(input%remfile,i)
        call open_file(unit_no, filename, IOFileInput)

        do j = 1, nlines
          read(unit_no,*) istep, parmid
          do k = 1, num_conv_id
            m = conv_id(k)
            if (parmid == m) then
              do_open_file(i) = .true.
            end if
          end do
          if (do_open_file(i)) exit
        end do

        call close_file(unit_no)
      end do

    else if (option%convert_type == ConvertTypeReplica) then

      do_open_file(:) = .false.
      do i = 1, option%num_replicas

        do k = 1, num_conv_id
          m = conv_id(k)
          if (m == i) then
            do_open_file(i) = .true.
          end if
        end do

      end do

    end if


    ! print open file index lists
    !
    write(MsgOut,'(A)')   'Convert> Replicas that contain "convert_ids"'
    write(MsgOut,'(A,$)') '  Replicas: '
    do i = 1, option%num_replicas
      if (do_open_file(i)) then
        write(MsgOut,'(I5,$)') i
      end if
    end do
    write(MsgOut,'(A)') ''


    ! open input files
    !
    trj_in(:)%unit_no = 0
    rem_unit_in(:) = 0
    log_unit_in(:) = 0
    ene_unit_in(:) = 0

    do i = 1, option%num_replicas
       if (do_open_file(i)) then

         ! assign file unit no for remfiles
         if (option%convert_type == ConvertTypeParameter) then
           filename = get_replicate_name(input%remfile,i)
           call open_file(unit_no, filename, IOFileInput)
           rem_unit_in(i) = unit_no
         end if

         ! assign file unit no for dcdfiles
         if (conv_crd) then
           filename = get_replicate_name(input%dcdfile,i)
           call open_trj(trj_in(i), filename, TrjFormatDCD, &
                         TrjTypeCoorBox, IOFileInput)
         end if

         ! assign file unit no for logfiles
         if (conv_log) then
           filename = get_replicate_name(input%logfile,i)
           call open_file(unit_no, filename, IOFileInput)
           log_unit_in(i) = unit_no
         end if

         ! assign file unit no for enefiles of gREST
         if (conv_ene) then
           filename = get_replicate_name(input%enefile,i)
           call open_file(unit_no, filename, IOFileInput)
           ene_unit_in(i) = unit_no
         end if

       end if
    end do


    ! open output files
    !
    trj_out(:)%unit_no = 0
    log_unit_out(:)    = 0
    ene_unit_out(:)    = 0

    do k = 1, num_conv_id
      m = conv_id(k)

      if (conv_crd) then
        filename = get_replicate_name(output%trjfile,m)
        call open_trj(trj_out(k), filename, option%trjout_format, &
                      option%trjout_type, IOFileOutputNew)
      end if

      if (conv_log) then
        filename = get_replicate_name(output%logfile,m)
        call open_file(unit_no, filename, IOFileOutputNew)
        log_unit_out(k) = unit_no
      end if

      if (conv_ene) then
        filename = get_replicate_name(output%enefile,m)
        call open_file(unit_no, filename, IOFileOutputNew)
        ene_unit_out(k) = unit_no
      end if

    end do


    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Convert> Assign Input/Output file unit no'
    write(MsgOut,*)     ' Input  remfiles :', rem_unit_in(:)
    if (conv_crd) &
      write(MsgOut,*)   ' Input  dcdfiles :', trj_in(:)%unit_no
    if (conv_log) &
      write(MsgOut,*)   ' Input  logfiles :', log_unit_in(:)
    if (conv_ene) &
      write(MsgOut,*)   ' Input  logfiles :', ene_unit_in(:)
    if (conv_crd) &
      write(MsgOut,*)   ' Output trjfiles :', trj_out(:)%unit_no
    if (conv_log) &
      write(MsgOut,*)   ' Output logfiles :', log_unit_out(:)
    if (conv_ene) &
      write(MsgOut,*)   ' Output logfiles :', ene_unit_out(:)
    write(MsgOut,'(A)') ''


    ! read 0 step parameter index
    !
    if (option%convert_type == ConvertTypeParameter) then

      param_index(:) = 0
      do i = 1, option%num_replicas
        if (do_open_file(i)) then
          read(rem_unit_in(i),*) istep, param_index(i)
        end if
      end do

    else if (option%convert_type == ConvertTypeReplica) then

      param_index(:) = 0
      do i = 1, option%num_replicas
        if (do_open_file(i)) then
          param_index(i) = i
        end if
      end do

    end if


    ! read(i)/write(j) logfile header
    !
    if (conv_log) then

      do i = 1, option%num_replicas
        if (do_open_file(i)) then

          ! check length of the header
          read (log_unit_in(i),'(A)') line
          if (line(1:17) == 'Generate_Velocity') then
            n = 7
          else
            n = 1
          end if

          do_skip = .true.
          do j = 1, num_conv_id
            m = conv_id(j)
            if (m == param_index(i)) then
              write(log_unit_out(j),'(A)') trim(line)
              do k = 1, n
                read (log_unit_in (i),'(A)') line
                write(log_unit_out(j),'(A)') trim(line)
              end do
              do_skip = .false.
            end if
          end do

          if (do_skip) then
            do k = 1, n
              read (log_unit_in(i),'(A)') line
            end do
          end if

        end if
      end do

    end if


    ! parameter id -> output slot. Fixed for the whole run.
    !
    allocate(dest(option%num_replicas))
    dest(:) = 0
    do j = 1, num_conv_id
      m = conv_id(j)
      if (m >= 1 .and. m <= option%num_replicas) then
        if (dest(m) == 0) dest(m) = j
      end if
    end do


    ! A static selection does not depend on the frame. "around" does.
    ! With no centering, fitting, or PBC correction and every atom kept,
    ! DCD frames can be copied as single precision.
    !
    static_sel = .true.
    sel_all    = .false.
    fast_dcd   = .false.
    if (conv_crd) then
      natm = size(trajectory%coord, 2)
      static_sel = .not. expr_has_around(option%trjout_atom_exp)
      if (static_sel .and. allocated(option%trjout_atom%idx)) then
        sel_all = (size(option%trjout_atom%idx) == natm)
        if (sel_all) then
          do i = 1, natm
            if (option%trjout_atom%idx(i) /= i) then
              sel_all = .false.
              exit
            end if
          end do
        end if
      end if
      fast_dcd = static_sel .and. sel_all .and. &
                 (.not. option%centering) .and. &
                 (option%pbcc_mode == PBCCModeNo) .and. &
                 (fitting%fitting_method == FittingMethodNO) .and. &
                 (option%trjout_format == TrjFormatDCD) .and. &
                 (natm > 0)
      if (fast_dcd) then
        allocate(xbuf(natm), ybuf(natm), zbuf(natm))
        allocate(in_ready(option%num_replicas), out_ready(num_conv_id))
        in_ready(:)  = .false.
        out_ready(:) = .false.
      end if
    end if


    ! main loop
    !
    DO imdstep = 1, option%nsteps

      ! read dcdfile(i) / write trjfile(j)
      !
      if (conv_crd) then

        if (mod(imdstep,option%crdout_period) == 0) then

          write(MsgOut,'(A,I10,A)') 'Convert> ', imdstep, ' step  read dcdfile > write trjfile'
          write(MsgOut,'(A,$)')     '         '

          do i = 1, option%num_replicas
            if (do_open_file(i)) then

              j = 0
              if (param_index(i) >= 1 .and. &
                  param_index(i) <= option%num_replicas) &
                j = dest(param_index(i))

              do_write = .false.
              if (j > 0) &
                do_write = (mod(imdstep, option%trjout_period) == 0)

              if (j > 0 .and. fast_dcd) then

                ! First frame of each file goes through read_trj / write_trj
                ! so the DCD header is consumed and written. Later frames are
                ! copied as single precision.
                if (in_ready(i) .and. &
                    .not. (do_write .and. .not. out_ready(j)) .and. &
                    .not. (do_write .and. &
                           trj_in(i)%trj_type == TrjTypeCoor .and. &
                           trj_out(j)%trj_type == TrjTypeCoorBox)) then
                  call transfer_dcd_frame(trj_in(i), xbuf, ybuf, zbuf, &
                                          trj_out(j), do_write)
                else
                  call read_trj(trj_in(i), trajectory)
                  in_ready(i) = .true.
                  if (do_write) then
                    call write_trj(trj_out(j), trajectory)
                    out_ready(j) = .true.
                  end if
                end if

                if (do_write) &
                  write(MsgOut,'(I5,A,I5,$)') &
                    trj_in(i)%unit_no, '>', trj_out(j)%unit_no

              else if (j > 0) then

                call read_trj(trj_in(i), trajectory)

                if (.not. static_sel) then
                  call reselect_atom(molecule,               &
                                     option%trjout_atom_exp, &
                                     trajectory%coord,       &
                                     option%trjout_atom,     &
                                     trjout_atom_trj)
                end if

                call centering(molecule, trajectory%coord, option)
                call run_pbc_correct(option%pbcc_mode, molecule, trajectory)

                if (fitting%mass_weight) then
                  call run_fitting(fitting,             &
                                   molecule%atom_coord, &
                                   trajectory%coord,    &
                                   trajectory%coord,    &
                                   molecule%mass)
                else
                  call run_fitting(fitting,             &
                                   molecule%atom_coord, &
                                   trajectory%coord,    &
                                   trajectory%coord)
                end if

                if (do_write) then
                  if (static_sel .and. sel_all .and. &
                      option%trjout_format == TrjFormatDCD) then
                    call write_trj(trj_out(j), trajectory)
                  else if (static_sel) then
                    call write_trj(trj_out(j), trajectory, &
                                   option%trjout_atom, molecule)
                  else
                    call write_trj(trj_out(j), trajectory, &
                                   trjout_atom_trj, molecule)
                  end if
                  write(MsgOut,'(I5,A,I5,$)') &
                    trj_in(i)%unit_no, '>', trj_out(j)%unit_no
                end if

              else

                if (fast_dcd .and. in_ready(i)) then
                  call transfer_dcd_frame(trj_in(i), xbuf, ybuf, zbuf)
                else
                  call read_trj(trj_in(i), trajectory)
                  if (fast_dcd) in_ready(i) = .true.
                end if

              end if

            end if
          end do

          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') ''

        end if

      end if


      ! read logfile(i) / write logfile(j)
      !
      if (conv_log) then

        if (mod(imdstep,option%eneout_period) == 0) then

          write(MsgOut,'(A,I10,A)') 'Convert> ', imdstep, ' step  read logfile > write logfile'
          write(MsgOut,'(A,$)')     '         '

          do i = 1, option%num_replicas
            if (do_open_file(i)) then

              j = 0
              if (param_index(i) >= 1 .and. &
                  param_index(i) <= option%num_replicas) &
                j = dest(param_index(i))

              if (j > 0) then

                  done = .false.
                  do while (.not. done)

                    read(log_unit_in(i),'(A)',iostat=ierr) line

                    if (ierr < 0) then
                      write(0,*) trim(line)
                      write(MsgOut,'(A,I10,A)') 'Convert> Something is wrong around ', imdstep, ' step in logfile'
                      write(MsgOut,'(A,I10,A,I10)') 'Convert> error: ', ierr, &
                                                    ' unit: ', log_unit_in(i)
                      call error_msg('Convert> Convert is failed due to unrecognized lines in logfile')
                    end if

                    if (line(1:5) == 'INFO:') then
                      if (mod(imdstep,option%logout_period) == 0) then
                        write(log_unit_out(j),'(A)') trim(line)
                        write(log_unit_out(j),'(A)') ''
                      endif
                      done = .true.
                    end if

                  end do

                  write(MsgOut,'(I5,A,I5,$)') log_unit_in(i),'>',log_unit_out(j)

              else

                done = .false.
                do while (.not. done)
                  read (log_unit_in(i),'(A)') line
                  if (line(1:5) == 'INFO:') then
                    read (log_unit_in(i),'(A)') line
                    done = .true.
                  end if
                end do

              end if

            end if
          end do

          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') ''

        end if

      end if


      ! read enefile(i) / write enefile(j) for gREST
      !
      if (conv_ene) then

        if (mod(imdstep,option%crdout_period) == 0) then

          write(MsgOut,'(A,I10,A)') 'Convert> ', imdstep, ' step  read enefile > write enefile'
          write(MsgOut,'(A,$)')     '         '

          do i = 1, option%num_replicas
            if (do_open_file(i)) then

              j = 0
              if (param_index(i) >= 1 .and. &
                  param_index(i) <= option%num_replicas) &
                j = dest(param_index(i))

              if (j > 0) then

                  done = .false.
                  do while (.not. done)

                    read(ene_unit_in(i),'(A)',iostat=ierr) line

                    if (ierr < 0) then
                      write(0,*) trim(line)
                      write(MsgOut,'(A,I10,A)') 'Convert> Something is wrong around ', imdstep, ' step in enefile'
                      write(MsgOut,'(A,I10,A,I10)') 'Convert> error: ', ierr, &
                                                    ' unit: ', ene_unit_in(i)
                      call error_msg('Convert> Convert is failed due to unrecognized lines in enefile')
                    end if


                    if (line(1:1) .ne. '#' .and. line(1:1) .ne. '@') then
                      write(ene_unit_out(j),'(A)') trim(line)
                      done = .true.
                    end if

                  end do

                  write(MsgOut,'(I5,A,I5,$)') ene_unit_in(i),'>',ene_unit_out(j)

              else

                done = .false.
                do while (.not. done)
                  read (ene_unit_in(i),'(A)') line
                  if (line(1:1) .ne. '#' .and. line(1:1) .ne. '@') then
                    done = .true.
                  end if
                end do

              end if

            end if
          end do

          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') ''

        end if

      end if


      ! read remfile(i)
      !
      if (mod(imdstep,option%exchange_period) == 0 .and. &
          imdstep /= option%nsteps) then

        if (option%convert_type == ConvertTypeParameter) then

          do i = 1, option%num_replicas
            if (do_open_file(i)) then
              read(rem_unit_in(i),*) istep, param_index(i)
            end if
          end do

        else if (option%convert_type == ConvertTypeReplica) then

          do i = 1, option%num_replicas
            if (do_open_file(i)) then
              param_index(i) = i
            end if
          end do

        end if

      end if

    END DO


    ! close files
    !
    do i = 1, option%num_replicas
      if (do_open_file(i)) then
        if (option%convert_type == ConvertTypeParameter) &
          call close_file(rem_unit_in(i))
        if (conv_crd) call close_trj(trj_in(i))
        if (conv_log) call close_file(log_unit_in(i))
        if (conv_ene) call close_file(ene_unit_in(i))
      end if
    end do

    do j = 1, num_conv_id
      if (conv_crd) call close_trj(trj_out(j))
      if (conv_log) call close_file(log_unit_out(j))
      if (conv_ene) call close_file(ene_unit_out(j))
    end do


    return

  end subroutine convert

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_convert
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_convert(input, output, option, conv_crd, conv_log, conv_ene)

    ! formal arguments
    type(s_input),           intent(in)    :: input
    type(s_output),          intent(in)    :: output
    type(s_option),          intent(in)    :: option
    logical,                 intent(out)   :: conv_crd
    logical,                 intent(out)   :: conv_log
    logical,                 intent(out)   :: conv_ene


    conv_crd = .false.
    conv_log = .false.
    conv_ene = .false.

    ! check remfile, dcdfile, logfile existing
    !
    if (option%convert_type == ConvertTypeParameter) then
      if (input%remfile == '') then
        call error_msg('Setup_Convert> Error: remfile is not specified in [INPUT]')
      end if
    end if

    if ((input%dcdfile == '' .and. output%trjfile /= '') .or. &
        (input%dcdfile /= '' .and. output%trjfile == '')) then
      call error_msg('Setup_Convert> dcdfile in [INPUT] or trjfile in [OUTPUT] is not specified')
    else if (input%dcdfile == '' .and. output%trjfile == '') then
      conv_crd = .false.
    else
      conv_crd = .true.
    end if

    if ((input%logfile == '' .and. output%logfile /= '') .or. &
        (input%logfile /= '' .and. output%logfile == '')) then
      call error_msg('Setup_Convert> logfile in [INPUT] or logfile in [OUTPUT] is not specified')
    else if (input%logfile == '' .and. output%logfile == '') then
      conv_log = .false.
    else
      conv_log = .true.
    end if

    if ((input%enefile == '' .and. output%enefile /= '') .or. &
        (input%enefile /= '' .and. output%enefile == '')) then
      call error_msg('Setup_Convert> enefile in [INPUT] or enefile in [OUTPUT] is not specified')
    else if (input%enefile == '' .and. output%enefile == '') then
      conv_ene = .false.
    else
      conv_ene = .true.
    end if

    if (conv_crd .and. option%crdout_period == 0) &
      call error_msg('Setup_Convert> error: crdout_period in [OPTION] = 0')

    if (conv_log .and. option%eneout_period == 0) &
      call error_msg('Setup_Convert> error: eneout_period in [OPTION] = 0')

    if (conv_ene .and. option%crdout_period == 0) &
      call error_msg('Setup_Convert> error: crdout_period in [OPTION] = 0')

    if (option%nsteps == 0) &
      call error_msg('Setup_Convert> error: nsteps in [OPTION] = 0')

    if (option%num_replicas == 0) &
      call error_msg('Setup_Convert> error: num_replicas in [OPTION] = 0')

    if (.not. conv_crd .and. .not. conv_log .and. .not. conv_ene) then
      write(MsgOut,'(A)') 'Setup_Convert> Nothing was done, because dcdfile, trjfile, logfile, and enefile are empty'
      write(MsgOut,'(A)') ''
      return
    end if

    return

  end subroutine setup_convert

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
  !  Function      get_replicate_name
  !> @brief        get replicate name
  !! @authors      NT
  !! @param[in]    filename   : filename
  !! @param[in]    replica_no : replica number
  !! @return       replicate name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_replicate_name(filename, replica_no)

    ! return
    character(Maxfilename)   :: get_replicate_name

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: replica_no

    ! local variables
    integer                  :: bl, br
    character(MaxLine)       :: bl_str, br_str
    character(20)            :: fmt
    integer                  :: i, ndigit, comp


    bl = index(filename, '{', back=.true.)
    br = index(filename, '}', back=.true.)

    if (bl == 0 .or. br == 0) then
      get_replicate_name = filename
      call error_msg('Get_Replicate_Name> Error: "{}" is not found in the file name')
    end if

    if (bl > br) &
      call error_msg('Get_Replicate_Name> Syntax error.')

    if (bl > 1) then
      bl_str = filename(1:bl-1)
    else
      bl_str = ''
    end if

    if (br < len_trim(filename)) then
      br_str = filename(br+1:)
    else
      br_str = ''
    end if
    
    do i = 1, 100
      comp = 10**i
      if (replica_no < comp) then
        ndigit = i
        exit
      end if
    end do
    
    write(fmt,'(a,i0,a,i0,a)') '(a,i', ndigit, '.', ndigit, ',a)'
    write(get_replicate_name,fmt) trim(bl_str), replica_no, trim(br_str)

    return

  end function get_replicate_name

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      expr_has_around
  !> @brief        true when a selection must be recomputed every frame
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function expr_has_around(expr) result(yes)

    ! formal arguments
    character(*), intent(in) :: expr

    ! return
    logical :: yes

    ! local variables
    character(len=len(expr)) :: tmp
    integer :: i, c


    tmp = expr
    do i = 1, len_trim(tmp)
      c = iachar(tmp(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) &
        tmp(i:i) = achar(c + 32)
    end do
    yes = (index(tmp, 'around') > 0)

    return

  end function expr_has_around

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    transfer_dcd_frame
  !> @brief        move one DCD frame without expanding it to double precision
  !! @param[inout] src      : input trajectory, header already read
  !! @param[inout] xbuf     : single-precision X buffer
  !! @param[inout] ybuf     : single-precision Y buffer
  !! @param[inout] zbuf     : single-precision Z buffer
  !! @param[inout] dst      : output trajectory (required when do_write)
  !! @param[in]    do_write : write the frame, or only advance the input
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine transfer_dcd_frame(src, xbuf, ybuf, zbuf, dst, do_write)

    ! formal arguments
    type(s_trj_file), intent(inout)           :: src
    real(sp),         intent(inout)           :: xbuf(:)
    real(sp),         intent(inout)           :: ybuf(:)
    real(sp),         intent(inout)           :: zbuf(:)
    type(s_trj_file), intent(inout), optional :: dst
    logical,          intent(in),    optional :: do_write

    ! local variables
    real(dp) :: symmat(6), pbcbox(6), symout(6)
    integer  :: ierr
    logical  :: write_frame, in_box, out_box


    write_frame = .false.
    if (present(do_write)) write_frame = do_write

    in_box = (src%trj_type == TrjTypeCoorBox)
    out_box = .false.
    if (write_frame) then
      if (.not. present(dst)) &
        call error_msg('Convert> internal error: missing output trajectory')
      out_box = (dst%trj_type == TrjTypeCoorBox)
    end if

    if (in_box) then
      read(src%unit_no, iostat=ierr) symmat
      if (ierr /= 0) &
        call error_msg('Convert> failed to read a DCD unit cell')
      if (write_frame .and. out_box) then
        call to_pbcbox(symmat, pbcbox)
        ! Same orthogonal cell written by write_dcd_trj / to_symmat.
        if (pbcbox(4) - EPS > 90.0_wp .or. &
            pbcbox(4) + EPS < 90.0_wp .or. &
            pbcbox(5) - EPS > 90.0_wp .or. &
            pbcbox(5) + EPS < 90.0_wp .or. &
            pbcbox(6) - EPS > 90.0_wp .or. &
            pbcbox(6) + EPS < 90.0_wp) &
          write(MsgOut,*) &
            'To_Symmat> WARNING: Non rectangular unit cell is not supported.'
        symout(1) = real(real(pbcbox(1), wp), dp)
        symout(2) = 0.0_dp
        symout(3) = real(real(pbcbox(2), wp), dp)
        symout(4) = 0.0_dp
        symout(5) = 0.0_dp
        symout(6) = real(real(pbcbox(3), wp), dp)
        write(dst%unit_no) symout
      end if
    end if

    read(src%unit_no, iostat=ierr) xbuf
    if (ierr /= 0) &
      call error_msg('Convert> failed to read DCD coordinates')
    read(src%unit_no, iostat=ierr) ybuf
    if (ierr /= 0) &
      call error_msg('Convert> failed to read DCD coordinates')
    read(src%unit_no, iostat=ierr) zbuf
    if (ierr /= 0) &
      call error_msg('Convert> failed to read DCD coordinates')

    if (write_frame) then
      write(dst%unit_no) xbuf
      write(dst%unit_no) ybuf
      write(dst%unit_no) zbuf
      dst%step_no = dst%step_no + 1
    end if
    src%step_no = src%step_no + 1

    return

  end subroutine transfer_dcd_frame

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remd_count_frames
  !> @brief        frames each parameter (or replica) id will receive
  !!               The walk matches convert(): one REMD step clock, a frame
  !!               read every crdout_period, a parameter update every
  !!               exchange_period after that frame is placed.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remd_count_frames(rem_packed, nreplica, name_len, &
                               convert_type, conv_ids, n_ids, &
                               nsteps, exchange_period, &
                               crdout_period, trjout_period, &
                               frame_counts, n_out)
    implicit none
    character(kind=c_char), intent(in)  :: rem_packed(*)
    integer,                intent(in)  :: nreplica
    integer,                intent(in)  :: name_len
    integer,                intent(in)  :: convert_type
    integer,                intent(in)  :: conv_ids(*)
    integer,                intent(in)  :: n_ids
    integer,                intent(in)  :: nsteps
    integer,                intent(in)  :: exchange_period
    integer,                intent(in)  :: crdout_period
    integer,                intent(in)  :: trjout_period
    integer,                intent(inout) :: frame_counts(n_out)
    integer,                intent(in)    :: n_out

    call remd_drive(rem_packed, rem_packed, nreplica, name_len, &
                    convert_type, conv_ids, n_ids, &
                    nsteps, exchange_period, crdout_period, trjout_period, &
                    conv_ids, 0, 0, &
                    frame_counts, n_out, do_fill=.false.)
  end subroutine remd_count_frames

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remd_fill_frames
  !> @brief        sort replica DCDs into caller-owned coordinate buffers
  !! @param        coords_ptrs  : n_out pointers, each (3, n_selected, nframe)
  !! @param        box_ptrs     : n_out pointers, each (3, 3, nframe)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remd_fill_frames(dcd_packed, rem_packed, nreplica, name_len, &
                              convert_type, conv_ids, n_ids, &
                              nsteps, exchange_period, &
                              crdout_period, trjout_period, &
                              selected, n_selected, natom, &
                              frame_counts, n_out, coords_ptrs, box_ptrs)
    implicit none
    character(kind=c_char), intent(in) :: dcd_packed(*)
    character(kind=c_char), intent(in) :: rem_packed(*)
    integer,                intent(in) :: nreplica
    integer,                intent(in) :: name_len
    integer,                intent(in) :: convert_type
    integer,                intent(in) :: conv_ids(*)
    integer,                intent(in) :: n_ids
    integer,                intent(in) :: nsteps
    integer,                intent(in) :: exchange_period
    integer,                intent(in) :: crdout_period
    integer,                intent(in) :: trjout_period
    integer,                intent(in) :: selected(*)
    integer,                intent(in) :: n_selected
    integer,                intent(in) :: natom
    integer,                intent(inout) :: frame_counts(n_out)
    integer,                intent(in) :: n_out
    type(c_ptr),            intent(in) :: coords_ptrs(n_out)
    type(c_ptr),            intent(in) :: box_ptrs(n_out)

    call remd_drive(dcd_packed, rem_packed, nreplica, name_len, &
                    convert_type, conv_ids, n_ids, &
                    nsteps, exchange_period, crdout_period, trjout_period, &
                    selected, n_selected, natom, &
                    frame_counts, n_out, coords_ptrs, box_ptrs, .true.)
  end subroutine remd_fill_frames

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remd_drive(dcd_packed, rem_packed, nreplica, name_len, &
                        convert_type, conv_ids, n_ids, &
                        nsteps, exchange_period, &
                        crdout_period, trjout_in, &
                        selected, n_selected, natom, &
                        frame_counts, n_out, coords_ptrs, box_ptrs, do_fill)
    implicit none
    character(kind=c_char), intent(in) :: dcd_packed(*)
    character(kind=c_char), intent(in) :: rem_packed(*)
    integer,                intent(in) :: nreplica
    integer,                intent(in) :: name_len
    integer,                intent(in) :: convert_type
    integer,                intent(in) :: conv_ids(*)
    integer,                intent(in) :: n_ids
    integer,                intent(in) :: nsteps
    integer,                intent(in) :: exchange_period
    integer,                intent(in) :: crdout_period
    integer,                intent(in) :: trjout_in
    integer,                intent(in) :: selected(*)
    integer,                intent(in) :: n_selected
    integer,                intent(in) :: natom
    integer,                intent(inout) :: frame_counts(n_out)
    integer,                intent(in)    :: n_out
    type(c_ptr),            intent(in), optional :: coords_ptrs(:)
    type(c_ptr),            intent(in), optional :: box_ptrs(:)
    logical,                intent(in)    :: do_fill

    integer :: i, j, m, imdstep, istep, parmid, unit_no, nlines
    integer :: trjout_period, n_conv, cursor
    integer, allocatable :: conv_id(:), dest(:), param_index(:)
    integer, allocatable :: filled(:), rem_unit(:)
    logical, allocatable :: do_open(:), in_ready(:)
    logical :: use_identity, write_this, got_box
    type(s_trj_file), allocatable :: trj_in(:)
    type(s_trajectory) :: trajectory
    real(sp), allocatable :: xbuf(:), ybuf(:), zbuf(:)
    real(wp), pointer :: coords(:,:,:), pbc(:,:,:)
    real(dp) :: llen(3)
    character(len=:), allocatable :: filename
    integer :: ierr


    if (nreplica <= 0) &
      call error_msg('Remd_Convert> num_replicas must be positive')
    if (nsteps <= 0) &
      call error_msg('Remd_Convert> nsteps must be positive')
    if (exchange_period <= 0 .or. crdout_period <= 0) &
      call error_msg('Remd_Convert> periods must be positive')
    if (convert_type /= ConvertTypeParameter .and. &
        convert_type /= ConvertTypeReplica) &
      call error_msg('Remd_Convert> convert_type must be PARAMETER or REPLICA')
    if (name_len <= 0) &
      call error_msg('Remd_Convert> empty filename')

    if (trjout_in <= 0) then
      trjout_period = crdout_period
    else
      trjout_period = trjout_in
    end if

    if (n_ids <= 0) then
      n_conv = nreplica
    else
      n_conv = n_ids
    end if
    if (n_out /= n_conv) &
      call error_msg('Remd_Convert> output count does not match convert_ids')

    allocate(conv_id(n_conv), dest(nreplica), param_index(nreplica))
    allocate(do_open(nreplica), filled(n_out))
    dest(:) = 0
    param_index(:) = 0
    do_open(:) = .false.
    filled(:) = 0
    use_identity = .false.
    if (.not. do_fill) frame_counts(:) = 0

    do j = 1, n_conv
      if (n_ids <= 0) then
        m = j
      else
        m = conv_ids(j)
      end if
      if (m < 1 .or. m > nreplica) &
        call error_msg('Remd_Convert> convert_ids is out of range')
      if (dest(m) /= 0) &
        call error_msg('Remd_Convert> convert_ids contains a duplicate')
      conv_id(j) = m
      dest(m) = j
    end do

    allocate(character(len=name_len) :: filename)

    ! Which replica files actually visit a requested id. Same early-exit
    ! scan as convert(), so an id that never appears is not read.
    if (convert_type == ConvertTypeReplica) then
      do i = 1, nreplica
        do_open(i) = (dest(i) /= 0)
      end do
    else
      nlines = nsteps / exchange_period
      do i = 1, nreplica
        call packed_name(rem_packed, i, name_len, filename)
        call open_file(unit_no, filename, IOFileInput)
        do j = 1, nlines
          read(unit_no, *, iostat=ierr) istep, parmid
          if (ierr /= 0) &
            call error_msg('Remd_Convert> remfile ended early')
          if (parmid >= 1 .and. parmid <= nreplica) then
            if (dest(parmid) /= 0) do_open(i) = .true.
          end if
          if (do_open(i)) exit
        end do
        call close_file(unit_no)
      end do
    end if

    allocate(rem_unit(nreplica), in_ready(nreplica))
    rem_unit(:) = 0
    in_ready(:) = .false.
    allocate(trj_in(nreplica))
    trj_in(:)%unit_no = 0

    if (convert_type == ConvertTypeParameter) then
      do i = 1, nreplica
        if (do_open(i)) then
          call packed_name(rem_packed, i, name_len, filename)
          call open_file(unit_no, filename, IOFileInput)
          rem_unit(i) = unit_no
          read(unit_no, *, iostat=ierr) istep, param_index(i)
          if (ierr /= 0) &
            call error_msg('Remd_Convert> remfile ended early')
        end if
      end do
    else
      do i = 1, nreplica
        if (do_open(i)) param_index(i) = i
      end do
    end if

    if (do_fill) then
      if (natom <= 0) &
        call error_msg('Remd_Convert> natom must be positive')
      if (n_selected <= 0) &
        call error_msg('Remd_Convert> selection is empty')
      use_identity = (n_selected == natom)
      do i = 1, n_selected
        if (selected(i) < 1 .or. selected(i) > natom) &
          call error_msg('Remd_Convert> selection index is out of range')
        if (selected(i) /= i) use_identity = .false.
      end do
      allocate(xbuf(natom), ybuf(natom), zbuf(natom))
      call alloc_trajectory(trajectory, natom)
      do i = 1, nreplica
        if (do_open(i)) then
          call packed_name(dcd_packed, i, name_len, filename)
          call open_trj(trj_in(i), filename, TrjFormatDCD, &
                        TrjTypeCoorBox, IOFileInput)
        end if
      end do
    end if

    do imdstep = 1, nsteps

      if (mod(imdstep, crdout_period) == 0) then
        do i = 1, nreplica
          if (.not. do_open(i)) cycle
          j = 0
          if (param_index(i) >= 1 .and. param_index(i) <= nreplica) &
            j = dest(param_index(i))
          write_this = (j > 0 .and. mod(imdstep, trjout_period) == 0)

          if (.not. do_fill) then
            if (write_this) frame_counts(j) = frame_counts(j) + 1
          else
            if (.not. in_ready(i)) then
              call read_trj(trj_in(i), trajectory)
              in_ready(i) = .true.
              if (write_this) then
                filled(j) = filled(j) + 1
                cursor = filled(j)
                if (cursor > frame_counts(j)) &
                  call error_msg('Remd_Convert> frame buffer is too small')
                call c_f_pointer(coords_ptrs(j), coords, &
                                 [3, n_selected, frame_counts(j)])
                call c_f_pointer(box_ptrs(j), pbc, &
                                 [3, 3, frame_counts(j)])
                call store_wp_frame(coords, pbc, cursor, &
                                    trajectory%coord, trajectory%pbc_box, &
                                    selected, n_selected, use_identity)
              end if
            else
              got_box = (trj_in(i)%trj_type == TrjTypeCoorBox)
              call read_sp_dcd_frame(trj_in(i)%unit_no, got_box, &
                                     xbuf, ybuf, zbuf, llen)
              if (write_this) then
                filled(j) = filled(j) + 1
                cursor = filled(j)
                if (cursor > frame_counts(j)) &
                  call error_msg('Remd_Convert> frame buffer is too small')
                call c_f_pointer(coords_ptrs(j), coords, &
                                 [3, n_selected, frame_counts(j)])
                call c_f_pointer(box_ptrs(j), pbc, &
                                 [3, 3, frame_counts(j)])
                call store_sp_frame(coords, pbc, cursor, &
                                    xbuf, ybuf, zbuf, llen, got_box, &
                                    selected, n_selected, use_identity)
              end if
            end if
          end if
        end do
      end if

      if (mod(imdstep, exchange_period) == 0 .and. imdstep /= nsteps) then
        if (convert_type == ConvertTypeParameter) then
          do i = 1, nreplica
            if (do_open(i)) then
              read(rem_unit(i), *, iostat=ierr) istep, param_index(i)
              if (ierr /= 0) &
                call error_msg('Remd_Convert> remfile ended early')
            end if
          end do
        else
          do i = 1, nreplica
            if (do_open(i)) param_index(i) = i
          end do
        end if
      end if

    end do

    do i = 1, nreplica
      if (rem_unit(i) /= 0) call close_file(rem_unit(i))
      if (do_fill .and. do_open(i)) call close_trj(trj_in(i))
    end do
    if (do_fill) call dealloc_trajectory(trajectory)

    return
  end subroutine remd_drive

  subroutine packed_name(packed, irep, name_len, filename)
    implicit none
    character(kind=c_char), intent(in)  :: packed(*)
    integer,                intent(in)  :: irep
    integer,                intent(in)  :: name_len
    character(*),           intent(out) :: filename
    integer :: i, offset

    filename = ''
    offset = (irep - 1) * name_len
    do i = 1, min(len(filename), name_len)
      if (packed(offset + i) == c_null_char) exit
      filename(i:i) = packed(offset + i)
    end do
  end subroutine packed_name

  subroutine read_sp_dcd_frame(unit_no, got_box, xbuf, ybuf, zbuf, llen)
    implicit none
    integer,  intent(in)    :: unit_no
    logical,  intent(in)    :: got_box
    real(sp), intent(inout) :: xbuf(:), ybuf(:), zbuf(:)
    real(dp), intent(out)   :: llen(3)
    real(dp) :: symmat(6), pbcbox(6)
    integer :: ierr

    llen(:) = 0.0_dp
    if (got_box) then
      read(unit_no, iostat=ierr) symmat
      if (ierr /= 0) call error_msg('Remd_Convert> failed to read a DCD unit cell')
      call to_pbcbox(symmat, pbcbox)
      llen(1) = pbcbox(1)
      llen(2) = pbcbox(2)
      llen(3) = pbcbox(3)
    end if
    read(unit_no, iostat=ierr) xbuf
    if (ierr /= 0) call error_msg('Remd_Convert> failed to read DCD coordinates')
    read(unit_no, iostat=ierr) ybuf
    if (ierr /= 0) call error_msg('Remd_Convert> failed to read DCD coordinates')
    read(unit_no, iostat=ierr) zbuf
    if (ierr /= 0) call error_msg('Remd_Convert> failed to read DCD coordinates')
  end subroutine read_sp_dcd_frame

  subroutine store_wp_frame(coords, pbc, iframe, coord, box, &
                            selected, n_selected, identity)
    implicit none
    real(wp), intent(inout) :: coords(:,:,:)
    real(wp), intent(inout) :: pbc(:,:,:)
    integer,  intent(in)    :: iframe
    real(wp), intent(in)    :: coord(:,:)
    real(wp), intent(in)    :: box(:,:)
    integer,  intent(in)    :: selected(*)
    integer,  intent(in)    :: n_selected
    logical,  intent(in)    :: identity
    integer :: ia, k

    if (identity) then
      coords(1, :, iframe) = coord(1, :)
      coords(2, :, iframe) = coord(2, :)
      coords(3, :, iframe) = coord(3, :)
    else
      do ia = 1, n_selected
        k = selected(ia)
        coords(:, ia, iframe) = coord(:, k)
      end do
    end if
    pbc(:, :, iframe) = box(:, :)
  end subroutine store_wp_frame

  subroutine store_sp_frame(coords, pbc, iframe, xbuf, ybuf, zbuf, llen, &
                            got_box, selected, n_selected, identity)
    implicit none
    real(wp), intent(inout) :: coords(:,:,:)
    real(wp), intent(inout) :: pbc(:,:,:)
    integer,  intent(in)    :: iframe
    real(sp), intent(in)    :: xbuf(:), ybuf(:), zbuf(:)
    real(dp), intent(in)    :: llen(3)
    logical,  intent(in)    :: got_box
    integer,  intent(in)    :: selected(*)
    integer,  intent(in)    :: n_selected
    logical,  intent(in)    :: identity
    integer :: ia, k

    if (identity) then
      coords(1, :, iframe) = real(xbuf(:), wp)
      coords(2, :, iframe) = real(ybuf(:), wp)
      coords(3, :, iframe) = real(zbuf(:), wp)
    else
      do ia = 1, n_selected
        k = selected(ia)
        coords(1, ia, iframe) = real(xbuf(k), wp)
        coords(2, ia, iframe) = real(ybuf(k), wp)
        coords(3, ia, iframe) = real(zbuf(k), wp)
      end do
    end if
    pbc(:, :, iframe) = 0.0_wp
    if (got_box) then
      pbc(1, 1, iframe) = real(llen(1), wp)
      pbc(2, 2, iframe) = real(llen(2), wp)
      pbc(3, 3, iframe) = real(llen(3), wp)
    end if
  end subroutine store_sp_frame

end module rc_convert_mod
