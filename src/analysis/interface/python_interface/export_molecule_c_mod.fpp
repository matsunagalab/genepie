module export_molecule_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod

  use molecules_str_mod
  implicit none
contains

  subroutine export_pdb_to_string_c( &
          molecule, out_pdb_ptr) &
      bind(C, name="export_pdb_to_string_c")
    use conv_f_c_util
    implicit none
    type(s_molecule_c), intent(in) :: molecule
    type(c_ptr), intent(out) :: out_pdb_ptr

    type(s_molecule) :: f_molecule
    character(len=:), allocatable :: out_pdb_f
    character(kind=c_char), pointer :: out_pdb_c(:)

    call c2f_s_molecule(molecule, f_molecule)
    call export_pdb_to_string(f_molecule, out_pdb_f)

    if (allocated(out_pdb_f)) then
      call f2c_string(out_pdb_f, out_pdb_c)
      out_pdb_ptr = c_loc(out_pdb_c(1))
    else
      out_pdb_ptr = c_null_ptr
    end if
  end subroutine export_pdb_to_string_c

  subroutine export_pdb_to_string( &
          molecule, out_pdb)
    use fileio_pdb_mod
    use internal_file_type_mod
    use molecules_mod
    implicit none
    type(s_molecule), intent(inout) :: molecule
    character(len=:), allocatable, intent(out) :: out_pdb
    type(s_pdb)              :: pdb_out

    call export_molecules(molecule, pdb=pdb_out)
    call write_pdb_to_string(out_pdb, pdb_out)
    call dealloc_pdb_all(pdb_out)
  end subroutine export_pdb_to_string
end module export_molecule_c_mod
