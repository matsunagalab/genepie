!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Module   selection_c_mod
!! @brief   C interface for GENESIS atom selection
!! @authors YS
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module selection_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use conv_f_c_util
  use select_atoms_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use constants_mod
  use error_mod, only: s_error, error_init, error_set, &
                       ERROR_SELECTION, error_to_c

  implicit none

  public :: selection_c
  public :: deallocate_selection_c

  ! Module-level pointer for selection results
  integer, pointer, save :: idx_ptr(:) => null()

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    selection_c
  !> @brief        C interface for atom selection using GENESIS selection syntax
  !! @param[in]    molecule_c      : molecule structure (C)
  !! @param[in]    selection_str   : selection expression (C string)
  !! @param[in]    selection_len   : length of selection string
  !! @param[out]   indices         : pointer to selected atom indices
  !! @param[out]   n_indices       : number of selected atoms
  !! @param[out]   status          : error status (0 = success)
  !! @param[out]   msg             : error message
  !! @param[in]    msglen          : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine selection_c(molecule_c, selection_str, selection_len, &
                         indices, n_indices, status, msg, msglen) &
      bind(C, name="selection_c")
    implicit none

    ! Arguments
    type(s_molecule_c), intent(in) :: molecule_c
    character(kind=c_char), intent(in) :: selection_str(*)
    integer(c_int), value :: selection_len
    type(c_ptr), intent(out) :: indices
    integer(c_int), intent(out) :: n_indices
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_molecule) :: molecule
    type(s_selatoms) :: selatoms
    character(len=1024) :: sel_str_f
    type(s_error) :: err
    integer :: i

    ! Initialize
    call error_init(err)
    status = 0
    indices = c_null_ptr
    n_indices = 0

    ! Convert C string to Fortran string
    call c2f_string(selection_str, sel_str_f)

    ! Convert s_molecule_c to s_molecule (creates a copy with char arrays)
    call c2f_s_molecule(molecule_c, molecule)

    ! Perform atom selection
    call select_atom(molecule, trim(sel_str_f), selatoms)

    ! Check if selection returned any atoms
    if (.not. allocated(selatoms%idx) .or. size(selatoms%idx) == 0) then
      call error_set(err, ERROR_SELECTION, &
                     "Selection returned no atoms: " // trim(sel_str_f))
      call error_to_c(err, status, msg, msglen)
      call dealloc_molecules_all(molecule)
      return
    end if

    ! Deallocate previous results if any
    if (associated(idx_ptr)) then
      deallocate(idx_ptr)
      nullify(idx_ptr)
    end if

    ! Copy results to module-level pointer
    n_indices = size(selatoms%idx)
    allocate(idx_ptr(n_indices))
    idx_ptr(:) = selatoms%idx(:)
    indices = c_loc(idx_ptr(1))

    ! Cleanup
    call dealloc_selatoms(selatoms)
    call dealloc_molecules_all(molecule)

  end subroutine selection_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_selection_c
  !> @brief        Deallocate selection results
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine deallocate_selection_c() bind(C, name="deallocate_selection_c")
    implicit none

    if (associated(idx_ptr)) then
      deallocate(idx_ptr)
      nullify(idx_ptr)
    end if

  end subroutine deallocate_selection_c

end module selection_c_mod
