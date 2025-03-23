module define_molecule
  use, intrinsic :: iso_c_binding
  use input_mod
  use molecules_mod
  use fileio_pdb_mod
  use fileio_psf_mod
  use constants_mod      ! Likely contains MaxFilename
  use molecules_str_mod  ! Contains s_molecule definition
  use s_molecule_c_mod
  implicit none

  private
  public :: define_molecule_from_file

  ! Define MaxFilename if it"s not available from constants_mod
  integer, parameter :: MaxFilename = 256  ! Adjust this value as needed

contains
  subroutine define_molecule_from_file( &
          pdb_filename, &
          top_filename, &
          gpr_filename, &
          psf_filename, &
          ref_filename, &
          fit_filename, &
          out_mol) &
      bind(C, name="define_molecule_from_file")
    use fileio_top_mod
    use fileio_par_mod
    use fileio_gpr_mod
    use fileio_pdb_mod
    use fileio_psf_mod
    use fileio_crd_mod
    use fileio_prmtop_mod
    use fileio_ambcrd_mod
    use fileio_grotop_mod
    use fileio_grocrd_mod
    use fileio_mode_mod
    use conv_f_c_util
    implicit none
    ! Input parameters
    character(kind=c_char), intent(in) :: pdb_filename(*)
    character(kind=c_char), intent(in) :: top_filename(*)
    character(kind=c_char), intent(in) :: gpr_filename(*)
    character(kind=c_char), intent(in) :: psf_filename(*)
    character(kind=c_char), intent(in) :: ref_filename(*)
    character(kind=c_char), intent(in) :: fit_filename(*)
    ! Output parameters
    type(s_molecule_c), intent(out) :: out_mol
    ! Local variables
    type(s_inp_info) :: inp_info
    type(s_pdb) :: pdb
    type(s_top) :: top
    type(s_gpr) :: gpr
    type(s_psf) :: psf
    type(s_pdb) :: ref
    type(s_pdb) :: fit
    type(s_molecule) :: molecule
    character(MaxFilename) :: filename

    if (pdb_filename(1) /= c_null_char) then
      call c2f_string(pdb_filename, filename)
      inp_info%pdbfile = trim(filename)
      call input_files(inp_info, pdb=pdb)
    end if
    if (top_filename(1) /= c_null_char) then
      call c2f_string(top_filename, filename)
      inp_info%topfile = trim(filename)
      call input_files(inp_info, top=top)
    end if
    if (gpr_filename(1) /= c_null_char) then
      call c2f_string(gpr_filename, filename)
      inp_info%gprfile = trim(filename)
      call input_files(inp_info, gpr=gpr)
    end if
    if (psf_filename(1) /= c_null_char) then
      call c2f_string(psf_filename, filename)
      inp_info%psffile = trim(filename)
      call input_files(inp_info, psf=psf)
    end if
    if (ref_filename(1) /= c_null_char) then
      call c2f_string(ref_filename, filename)
      inp_info%reffile = trim(filename)
      call input_files(inp_info, ref=ref)
    end if
    if (fit_filename(1) /= c_null_char) then
      call c2f_string(fit_filename, filename)
      inp_info%fitfile = trim(filename)
      call input_files(inp_info, fit=fit)
    end if
    call define_molecules(molecule, &
        pdb=pdb, top=top, gpr=gpr, psf=psf, ref=ref, fit=fit)
    call f2c_s_molecule(molecule, out_mol)
  end subroutine define_molecule_from_file

end module define_molecule
