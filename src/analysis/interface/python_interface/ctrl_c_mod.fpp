module ctrl_c_mod
  use iso_c_binding
  implicit none

  type, bind(C) :: c_s_inp_info
    type(c_ptr) :: ambcrdfile  = C_NULL_PTR ! char*
    type(c_ptr) :: ambreffile  = C_NULL_PTR ! char*
    type(c_ptr) :: atpfile     = C_NULL_PTR ! char*
    type(c_ptr) :: coorfile    = C_NULL_PTR ! char*
    type(c_ptr) :: cvfile      = C_NULL_PTR ! char*
    type(c_ptr) :: targetfile  = C_NULL_PTR ! char*
    type(c_ptr) :: dcdfile     = C_NULL_PTR ! char*
    type(c_ptr) :: dcdvelfile  = C_NULL_PTR ! char*
    type(c_ptr) :: excfile     = C_NULL_PTR ! char*
    type(c_ptr) :: gprfile     = C_NULL_PTR ! char*
    type(c_ptr) :: grocrdfile  = C_NULL_PTR ! char*
    type(c_ptr) :: groreffile  = C_NULL_PTR ! char*
    type(c_ptr) :: grotopfile  = C_NULL_PTR ! char*
    type(c_ptr) :: mtfile      = C_NULL_PTR ! char*
    type(c_ptr) :: indexfile   = C_NULL_PTR ! char*
    type(c_ptr) :: logfile     = C_NULL_PTR ! char*
    type(c_ptr) :: enefile     = C_NULL_PTR ! char*
    type(c_ptr) :: msdfile     = C_NULL_PTR ! char*
    type(c_ptr) :: pathfile    = C_NULL_PTR ! char*
    type(c_ptr) :: pathcvfile  = C_NULL_PTR ! char*
    type(c_ptr) :: pcafile     = C_NULL_PTR ! char*
    type(c_ptr) :: pdbfile     = C_NULL_PTR ! char*
    type(c_ptr) :: pdb_tgtfile = C_NULL_PTR ! char*
    type(c_ptr) :: pdb_avefile = C_NULL_PTR ! char*
    type(c_ptr) :: pdb_aftfile = C_NULL_PTR ! char*
    type(c_ptr) :: pdb_sphfile = C_NULL_PTR ! char*
    type(c_ptr) :: pdb_wbxfile = C_NULL_PTR ! char*
    type(c_ptr) :: prmtopfile  = C_NULL_PTR ! char*
    type(c_ptr) :: psffile     = C_NULL_PTR ! char*
    type(c_ptr) :: radfile     = C_NULL_PTR ! char*
    type(c_ptr) :: refenefile  = C_NULL_PTR ! char*
    type(c_ptr) :: reffile     = C_NULL_PTR ! char*
    type(c_ptr) :: fitfile     = C_NULL_PTR ! char*
    type(c_ptr) :: remfile     = C_NULL_PTR ! char*
    type(c_ptr) :: rstfile     = C_NULL_PTR ! char*
    type(c_ptr) :: rtpfile     = C_NULL_PTR ! char*
    type(c_ptr) :: topfile     = C_NULL_PTR ! char*
    type(c_ptr) :: valfile     = C_NULL_PTR ! char*
    type(c_ptr) :: vecfile     = C_NULL_PTR ! char*
    type(c_ptr) :: velfile     = C_NULL_PTR ! char*
    type(c_ptr) :: weightfile  = C_NULL_PTR ! char*
    type(c_ptr) :: xscfile     = C_NULL_PTR ! char*
    type(c_ptr) :: distfile    = C_NULL_PTR ! char*
  end type c_s_inp_info

  type, bind(C) :: c_s_out_info
    type(c_ptr) :: ambcrdfile   = C_NULL_PTR ! char*
    type(c_ptr) :: angfile      = C_NULL_PTR ! char*
    type(c_ptr) :: cntfile      = C_NULL_PTR ! char*
    type(c_ptr) :: comangfile   = C_NULL_PTR ! char*
    type(c_ptr) :: comdisfile   = C_NULL_PTR ! char*
    type(c_ptr) :: comtorfile   = C_NULL_PTR ! char*
    type(c_ptr) :: coorfile     = C_NULL_PTR ! char*
    type(c_ptr) :: crdfile      = C_NULL_PTR ! char*
    type(c_ptr) :: crsfile      = C_NULL_PTR ! char*
    type(c_ptr) :: disfile      = C_NULL_PTR ! char*
    type(c_ptr) :: enefile      = C_NULL_PTR ! char*
    type(c_ptr) :: excfile      = C_NULL_PTR ! char*
    type(c_ptr) :: fenefile     = C_NULL_PTR ! char*
    type(c_ptr) :: gprfile      = C_NULL_PTR ! char*
    type(c_ptr) :: grotopfile   = C_NULL_PTR ! char*
    type(c_ptr) :: grocrdfile   = C_NULL_PTR ! char*
    type(c_ptr) :: grocrd_tgtfile = C_NULL_PTR ! char*
    type(c_ptr) :: hb_listfile  = C_NULL_PTR ! char*
    type(c_ptr) :: indexfile    = C_NULL_PTR ! char*
    type(c_ptr) :: logfile      = C_NULL_PTR ! char*
    type(c_ptr) :: mapfile      = C_NULL_PTR ! char*
    type(c_ptr) :: morphfile    = C_NULL_PTR ! char*
    type(c_ptr) :: msdfile      = C_NULL_PTR ! char*
    type(c_ptr) :: outfile      = C_NULL_PTR ! char*
    type(c_ptr) :: parfile      = C_NULL_PTR ! char*
    type(c_ptr) :: pathcvfile   = C_NULL_PTR ! char*
    type(c_ptr) :: pcafile      = C_NULL_PTR ! char*
    type(c_ptr) :: pdbfile      = C_NULL_PTR ! char*
    type(c_ptr) :: pdb_avefile  = C_NULL_PTR ! char*
    type(c_ptr) :: pdb_aftfile  = C_NULL_PTR ! char*
    type(c_ptr) :: pdb_tgtfile  = C_NULL_PTR ! char*
    type(c_ptr) :: pmffile      = C_NULL_PTR ! char*
    type(c_ptr) :: pmlfile      = C_NULL_PTR ! char*
    type(c_ptr) :: prjfile      = C_NULL_PTR ! char*
    type(c_ptr) :: probfile     = C_NULL_PTR ! char*
    type(c_ptr) :: qmmm_crdfile = C_NULL_PTR ! char*
    type(c_ptr) :: qmmm_psffile = C_NULL_PTR ! char*
    type(c_ptr) :: qmmm_pdbfile = C_NULL_PTR ! char*
    type(c_ptr) :: qntfile      = C_NULL_PTR ! char*
    type(c_ptr) :: rdffile      = C_NULL_PTR ! char*
    type(c_ptr) :: rmsfile      = C_NULL_PTR ! char*
    type(c_ptr) :: rgfile       = C_NULL_PTR ! char*
    type(c_ptr) :: rstfile      = C_NULL_PTR ! char*
    type(c_ptr) :: tblfile      = C_NULL_PTR ! char*
    type(c_ptr) :: txtfile      = C_NULL_PTR ! char*
    type(c_ptr) :: topfile      = C_NULL_PTR ! char*
    type(c_ptr) :: torfile      = C_NULL_PTR ! char*
    type(c_ptr) :: trjfile      = C_NULL_PTR ! char*
    type(c_ptr) :: trrfile      = C_NULL_PTR ! char*
    type(c_ptr) :: valfile      = C_NULL_PTR ! char*
    type(c_ptr) :: vcvfile      = C_NULL_PTR ! char*
    type(c_ptr) :: vecfile      = C_NULL_PTR ! char*
    type(c_ptr) :: velfile      = C_NULL_PTR ! char*
    type(c_ptr) :: voronoifile  = C_NULL_PTR ! char*
    type(c_ptr) :: vmdfile      = C_NULL_PTR ! char*
    type(c_ptr) :: weightfile   = C_NULL_PTR ! char*
    type(c_ptr) :: xscfile      = C_NULL_PTR ! char*
  end type c_s_out_info

  type, bind(C) :: c_s_trj_info
    type(c_ptr) :: trj_files      = C_NULL_PTR ! char**
    integer(c_int) :: n_trj_files = 0          ! size of trj_files array
    type(c_ptr) :: md_steps       = C_NULL_PTR ! int*
    integer(c_int) :: n_md_steps  = 0          ! size of md_steps array
    type(c_ptr) :: mdout_periods  = C_NULL_PTR ! int*
    integer(c_int) :: n_mdout_periods = 0      ! size of mdout_periods array
    type(c_ptr) :: ana_periods    = C_NULL_PTR ! int*
    integer(c_int) :: n_ana_periods = 0        ! size of ana_periods array
    type(c_ptr) :: start_steps    = C_NULL_PTR ! int*
    integer(c_int) :: n_start_steps = 0        ! size of start_steps array
    integer(c_int) :: trj_format   = 0
    integer(c_int) :: trj_type     = 0
    integer(c_int) :: trj_natom    = 0
  end type c_s_trj_info

  type, bind(C) :: c_s_sel_info
    type(c_ptr) :: groups      = C_NULL_PTR  ! char**
    integer(c_int) :: n_groups = 0           ! size of groups array
    type(c_ptr) :: mole_names  = C_NULL_PTR  ! char**
    integer(c_int) :: n_mole_names = 0       ! size of mole_names array
  end type c_s_sel_info

contains

  subroutine c2f_inp_info(c_inp, f_inp)
    use input_mod
    use conv_f_c_util
    type(c_s_inp_info), intent(in) :: c_inp
    type(s_inp_info), intent(out) :: f_inp

    call cptr_to_fstring(c_inp%ambcrdfile, f_inp%ambcrdfile)
    call cptr_to_fstring(c_inp%ambreffile, f_inp%ambreffile)
    call cptr_to_fstring(c_inp%atpfile, f_inp%atpfile)
    call cptr_to_fstring(c_inp%coorfile, f_inp%coorfile)
    call cptr_to_fstring(c_inp%cvfile, f_inp%cvfile)
    call cptr_to_fstring(c_inp%targetfile, f_inp%targetfile)
    call cptr_to_fstring(c_inp%dcdfile, f_inp%dcdfile)
    call cptr_to_fstring(c_inp%dcdvelfile, f_inp%dcdvelfile)
    call cptr_to_fstring(c_inp%excfile, f_inp%excfile)
    call cptr_to_fstring(c_inp%gprfile, f_inp%gprfile)
    call cptr_to_fstring(c_inp%grocrdfile, f_inp%grocrdfile)
    call cptr_to_fstring(c_inp%groreffile, f_inp%groreffile)
    call cptr_to_fstring(c_inp%grotopfile, f_inp%grotopfile)
    call cptr_to_fstring(c_inp%mtfile, f_inp%mtfile)
    call cptr_to_fstring(c_inp%indexfile, f_inp%indexfile)
    call cptr_to_fstring(c_inp%logfile, f_inp%logfile)
    call cptr_to_fstring(c_inp%enefile, f_inp%enefile)
    call cptr_to_fstring(c_inp%msdfile, f_inp%msdfile)
    call cptr_to_fstring(c_inp%pathfile, f_inp%pathfile)
    call cptr_to_fstring(c_inp%pathcvfile, f_inp%pathcvfile)
    call cptr_to_fstring(c_inp%pcafile, f_inp%pcafile)
    call cptr_to_fstring(c_inp%pdbfile, f_inp%pdbfile)
    call cptr_to_fstring(c_inp%pdb_tgtfile, f_inp%pdb_tgtfile)
    call cptr_to_fstring(c_inp%pdb_avefile, f_inp%pdb_avefile)
    call cptr_to_fstring(c_inp%pdb_aftfile, f_inp%pdb_aftfile)
    call cptr_to_fstring(c_inp%pdb_sphfile, f_inp%pdb_sphfile)
    call cptr_to_fstring(c_inp%pdb_wbxfile, f_inp%pdb_wbxfile)
    call cptr_to_fstring(c_inp%prmtopfile, f_inp%prmtopfile)
    call cptr_to_fstring(c_inp%psffile, f_inp%psffile)
    call cptr_to_fstring(c_inp%radfile, f_inp%radfile)
    call cptr_to_fstring(c_inp%refenefile, f_inp%refenefile)
    call cptr_to_fstring(c_inp%reffile, f_inp%reffile)
    call cptr_to_fstring(c_inp%fitfile, f_inp%fitfile)
    call cptr_to_fstring(c_inp%remfile, f_inp%remfile)
    call cptr_to_fstring(c_inp%rstfile, f_inp%rstfile)
    call cptr_to_fstring(c_inp%rtpfile, f_inp%rtpfile)
    call cptr_to_fstring(c_inp%topfile, f_inp%topfile)
    call cptr_to_fstring(c_inp%valfile, f_inp%valfile)
    call cptr_to_fstring(c_inp%vecfile, f_inp%vecfile)
    call cptr_to_fstring(c_inp%velfile, f_inp%velfile)
    call cptr_to_fstring(c_inp%weightfile, f_inp%weightfile)
    call cptr_to_fstring(c_inp%xscfile, f_inp%xscfile)
    call cptr_to_fstring(c_inp%distfile, f_inp%distfile)
  end subroutine c2f_inp_info

  subroutine c2f_out_info(c_out, f_out)
    use output_mod
    use conv_f_c_util
    type(c_s_out_info), intent(in) :: c_out
    type(s_out_info), intent(out) :: f_out

    call cptr_to_fstring(c_out%ambcrdfile, f_out%ambcrdfile)
    call cptr_to_fstring(c_out%angfile, f_out%angfile)
    call cptr_to_fstring(c_out%cntfile, f_out%cntfile)
    call cptr_to_fstring(c_out%comangfile, f_out%comangfile)
    call cptr_to_fstring(c_out%comdisfile, f_out%comdisfile)
    call cptr_to_fstring(c_out%comtorfile, f_out%comtorfile)
    call cptr_to_fstring(c_out%coorfile, f_out%coorfile)
    call cptr_to_fstring(c_out%crdfile, f_out%crdfile)
    call cptr_to_fstring(c_out%crsfile, f_out%crsfile)
    call cptr_to_fstring(c_out%disfile, f_out%disfile)
    call cptr_to_fstring(c_out%enefile, f_out%enefile)
    call cptr_to_fstring(c_out%excfile, f_out%excfile)
    call cptr_to_fstring(c_out%fenefile, f_out%fenefile)
    call cptr_to_fstring(c_out%gprfile, f_out%gprfile)
    call cptr_to_fstring(c_out%grotopfile, f_out%grotopfile)
    call cptr_to_fstring(c_out%grocrdfile, f_out%grocrdfile)
    call cptr_to_fstring(c_out%grocrd_tgtfile, f_out%grocrd_tgtfile)
    call cptr_to_fstring(c_out%hb_listfile, f_out%hb_listfile)
    call cptr_to_fstring(c_out%indexfile, f_out%indexfile)
    call cptr_to_fstring(c_out%logfile, f_out%logfile)
    call cptr_to_fstring(c_out%mapfile, f_out%mapfile)
    call cptr_to_fstring(c_out%morphfile, f_out%morphfile)
    call cptr_to_fstring(c_out%msdfile, f_out%msdfile)
    call cptr_to_fstring(c_out%outfile, f_out%outfile)
    call cptr_to_fstring(c_out%parfile, f_out%parfile)
    call cptr_to_fstring(c_out%pathcvfile, f_out%pathcvfile)
    call cptr_to_fstring(c_out%pcafile, f_out%pcafile)
    call cptr_to_fstring(c_out%pdbfile, f_out%pdbfile)
    call cptr_to_fstring(c_out%pdb_avefile, f_out%pdb_avefile)
    call cptr_to_fstring(c_out%pdb_aftfile, f_out%pdb_aftfile)
    call cptr_to_fstring(c_out%pdb_tgtfile, f_out%pdb_tgtfile)
    call cptr_to_fstring(c_out%pmffile, f_out%pmffile)
    call cptr_to_fstring(c_out%pmlfile, f_out%pmlfile)
    call cptr_to_fstring(c_out%prjfile, f_out%prjfile)
    call cptr_to_fstring(c_out%probfile, f_out%probfile)
    call cptr_to_fstring(c_out%qmmm_crdfile, f_out%qmmm_crdfile)
    call cptr_to_fstring(c_out%qmmm_psffile, f_out%qmmm_psffile)
    call cptr_to_fstring(c_out%qmmm_pdbfile, f_out%qmmm_pdbfile)
    call cptr_to_fstring(c_out%qntfile, f_out%qntfile)
    call cptr_to_fstring(c_out%rdffile, f_out%rdffile)
    call cptr_to_fstring(c_out%rmsfile, f_out%rmsfile)
    call cptr_to_fstring(c_out%rgfile, f_out%rgfile)
    call cptr_to_fstring(c_out%rstfile, f_out%rstfile)
    call cptr_to_fstring(c_out%tblfile, f_out%tblfile)
    call cptr_to_fstring(c_out%txtfile, f_out%txtfile)
    call cptr_to_fstring(c_out%topfile, f_out%topfile)
    call cptr_to_fstring(c_out%torfile, f_out%torfile)
    call cptr_to_fstring(c_out%trjfile, f_out%trjfile)
    call cptr_to_fstring(c_out%trrfile, f_out%trrfile)
    call cptr_to_fstring(c_out%valfile, f_out%valfile)
    call cptr_to_fstring(c_out%vcvfile, f_out%vcvfile)
    call cptr_to_fstring(c_out%vecfile, f_out%vecfile)
    call cptr_to_fstring(c_out%velfile, f_out%velfile)
    call cptr_to_fstring(c_out%voronoifile, f_out%voronoifile)
    call cptr_to_fstring(c_out%vmdfile, f_out%vmdfile)
    call cptr_to_fstring(c_out%weightfile, f_out%weightfile)
    call cptr_to_fstring(c_out%xscfile, f_out%xscfile)
  end subroutine c2f_out_info

  subroutine c2f_trj_info(c_trj, f_trj)
    use string_mod
    use trajectory_mod
    use conv_f_c_util
    type(c_s_trj_info), intent(in) :: c_trj
    type(s_trj_info), intent(out) :: f_trj
    character(kind=c_char), pointer :: c_trj_files_array(:,:) => null()
    integer :: i

    if (c_associated(c_trj%trj_files)) then
      call c_f_pointer(c_trj%trj_files, c_trj_files_array, [c_trj%n_trj_files, MaxFilename])
      allocate(character(len=MaxFilename) :: f_trj%trj_files(c_trj%n_trj_files))
      do i = 1, c_trj%n_trj_files
        call cptr_to_fstring(c_loc(c_trj_files_array(i, :)), f_trj%trj_files(i))
      end do
    else
      allocate(character(len=MaxFilename) :: f_trj%trj_files(0))
    end if

    if (c_associated(c_trj%md_steps)) then
      call c2f_int_array(f_trj%md_steps, c_trj%md_steps, [c_trj%n_md_steps])
    else
      allocate(f_trj%md_steps(0))
    end if
    if (c_associated(c_trj%mdout_periods)) then
      call c2f_int_array(f_trj%mdout_periods, c_trj%mdout_periods, [c_trj%n_mdout_periods])
    else
      allocate(f_trj%mdout_periods(0))
    end if
    if (c_associated(c_trj%ana_periods)) then
      call c2f_int_array(f_trj%ana_periods, c_trj%ana_periods, [c_trj%n_ana_periods])
    else
      allocate(f_trj%ana_periods(0))
    end if
    if (c_associated(c_trj%start_steps)) then
      call c2f_int_array(f_trj%start_steps, c_trj%start_steps, [c_trj%n_start_steps])
    else
      allocate(f_trj%start_steps(0))
    end if

    f_trj%trj_format = c_trj%trj_format
    f_trj%trj_type   = c_trj%trj_type
    f_trj%trj_natom  = c_trj%trj_natom
  end subroutine c2f_trj_info

  subroutine c2f_sel_info(c_inp, inp)
    use string_mod
    use select_mod
    use conv_f_c_util
    type(c_s_sel_info), intent(in) :: c_inp
    type(s_sel_info), intent(out) :: inp
    character(kind=c_char), pointer :: c_groups_array(:,:) => null()
    character(kind=c_char), pointer :: c_mole_names_array(:,:) => null()
    integer :: i

    if (c_associated(c_inp%groups)) then
      call c_f_pointer(c_inp%groups, c_groups_array, [c_inp%n_groups, MaxLineLong])
      allocate(character(len=MaxLineLong) :: inp%groups(c_inp%n_groups))

      do i = 1, c_inp%n_groups
        call cptr_to_fstring(c_loc(c_groups_array(i, :)), inp%groups(i))
      end do
    else
      allocate(character(len=MaxLineLong) :: inp%groups(0))
    end if

    if (c_associated(c_inp%mole_names)) then
      call c_f_pointer(c_inp%mole_names, c_mole_names_array, [c_inp%n_mole_names, MaxLine])
      allocate(character(len=MaxLine) :: inp%mole_names(c_inp%n_mole_names))
      do i = 1, c_inp%n_mole_names
        call cptr_to_fstring(c_loc(c_mole_names_array(i, :)), inp%mole_names(i))
      end do
    else
      allocate(character(len=MaxLine) :: inp%mole_names(0))
    end if
  end subroutine c2f_sel_info

end module ctrl_c_mod
