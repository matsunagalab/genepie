!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Module   remd_convert_c_mod
!! @brief   C interface for sorting REMD trajectories into memory
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module remd_convert_c_mod
  use, intrinsic :: iso_c_binding
  use rc_convert_mod, only: remd_count_frames, remd_fill_frames
  use error_mod
  implicit none
  private

  public :: remd_convert_counts_c
  public :: remd_convert_counts_body
  public :: remd_convert_fill_c
  public :: remd_convert_fill_body

  type :: t_remd_counts_ctx
    type(c_ptr)    :: rem_ptr = c_null_ptr
    integer(c_int) :: nreplica = 0
    integer(c_int) :: name_len = 0
    integer(c_int) :: convert_type = 0
    type(c_ptr)    :: ids_ptr = c_null_ptr
    integer(c_int) :: n_ids = 0
    integer(c_int) :: nsteps = 0
    integer(c_int) :: exchange_period = 0
    integer(c_int) :: crdout_period = 0
    integer(c_int) :: trjout_period = 0
    type(c_ptr)    :: counts_ptr = c_null_ptr
    integer(c_int) :: n_out = 0
    type(s_error)  :: err
  end type t_remd_counts_ctx

  type :: t_remd_fill_ctx
    type(c_ptr)    :: dcd_ptr = c_null_ptr
    type(c_ptr)    :: rem_ptr = c_null_ptr
    integer(c_int) :: nreplica = 0
    integer(c_int) :: name_len = 0
    integer(c_int) :: convert_type = 0
    type(c_ptr)    :: ids_ptr = c_null_ptr
    integer(c_int) :: n_ids = 0
    integer(c_int) :: nsteps = 0
    integer(c_int) :: exchange_period = 0
    integer(c_int) :: crdout_period = 0
    integer(c_int) :: trjout_period = 0
    type(c_ptr)    :: selected_ptr = c_null_ptr
    integer(c_int) :: n_selected = 0
    integer(c_int) :: natom = 0
    type(c_ptr)    :: counts_ptr = c_null_ptr
    integer(c_int) :: n_out = 0
    type(c_ptr)    :: coords_ptr = c_null_ptr
    type(c_ptr)    :: box_ptr = c_null_ptr
    type(s_error)  :: err
  end type t_remd_fill_ctx

contains

  subroutine remd_convert_counts_c( &
          rem_names, nreplica, name_len, &
          convert_type, convert_ids, n_ids, &
          nsteps, exchange_period, crdout_period, trjout_period, &
          frame_counts, n_out, &
          status, msg, msglen) &
          bind(C, name="remd_convert_counts_c")
    implicit none
    character(kind=c_char), intent(in), target :: rem_names(*)
    integer(c_int), value :: nreplica, name_len, convert_type
    integer(c_int), value :: n_ids, nsteps
    integer(c_int), value :: exchange_period, crdout_period, trjout_period
    integer(c_int), value :: n_out, msglen
    type(c_ptr), value :: convert_ids
    type(c_ptr), value :: frame_counts
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    type(t_remd_counts_ctx), target :: c

    c%rem_ptr = c_loc(rem_names)
    c%nreplica = nreplica
    c%name_len = name_len
    c%convert_type = convert_type
    c%ids_ptr = convert_ids
    c%n_ids = n_ids
    c%nsteps = nsteps
    c%exchange_period = exchange_period
    c%crdout_period = crdout_period
    c%trjout_period = trjout_period
    c%counts_ptr = frame_counts
    c%n_out = n_out
    call run_guarded(remd_convert_counts_body, c, c%err, status, msg, msglen)
  end subroutine remd_convert_counts_c

  subroutine remd_convert_counts_body(ctx) bind(C, name="remd_convert_counts_body")
    implicit none
    type(c_ptr), value :: ctx
    type(t_remd_counts_ctx), pointer :: c
    character(kind=c_char), pointer :: rem(:)
    integer(c_int), pointer :: ids(:), counts(:)
    integer :: dummy_id(1)

    call c_f_pointer(ctx, c)
    call c_f_pointer(c%rem_ptr, rem, [c%nreplica * c%name_len])
    call c_f_pointer(c%counts_ptr, counts, [c%n_out])
    dummy_id(1) = 1
    if (c%n_ids > 0) then
      call c_f_pointer(c%ids_ptr, ids, [c%n_ids])
      call remd_count_frames(rem, c%nreplica, c%name_len, &
                             c%convert_type, ids, c%n_ids, &
                             c%nsteps, c%exchange_period, &
                             c%crdout_period, c%trjout_period, &
                             counts, c%n_out)
    else
      call remd_count_frames(rem, c%nreplica, c%name_len, &
                             c%convert_type, dummy_id, 0, &
                             c%nsteps, c%exchange_period, &
                             c%crdout_period, c%trjout_period, &
                             counts, c%n_out)
    end if
  end subroutine remd_convert_counts_body

  subroutine remd_convert_fill_c( &
          dcd_names, rem_names, nreplica, name_len, &
          convert_type, convert_ids, n_ids, &
          nsteps, exchange_period, crdout_period, trjout_period, &
          selected, n_selected, natom, &
          frame_counts, n_out, coords_ptrs, box_ptrs, &
          status, msg, msglen) &
          bind(C, name="remd_convert_fill_c")
    implicit none
    character(kind=c_char), intent(in), target :: dcd_names(*)
    character(kind=c_char), intent(in), target :: rem_names(*)
    integer(c_int), value :: nreplica, name_len, convert_type
    integer(c_int), value :: n_ids, nsteps
    integer(c_int), value :: exchange_period, crdout_period, trjout_period
    integer(c_int), value :: n_selected, natom, n_out, msglen
    type(c_ptr), value :: convert_ids, selected, frame_counts
    type(c_ptr), value :: coords_ptrs, box_ptrs
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    type(t_remd_fill_ctx), target :: c

    c%dcd_ptr = c_loc(dcd_names)
    c%rem_ptr = c_loc(rem_names)
    c%nreplica = nreplica
    c%name_len = name_len
    c%convert_type = convert_type
    c%ids_ptr = convert_ids
    c%n_ids = n_ids
    c%nsteps = nsteps
    c%exchange_period = exchange_period
    c%crdout_period = crdout_period
    c%trjout_period = trjout_period
    c%selected_ptr = selected
    c%n_selected = n_selected
    c%natom = natom
    c%counts_ptr = frame_counts
    c%n_out = n_out
    c%coords_ptr = coords_ptrs
    c%box_ptr = box_ptrs
    call run_guarded(remd_convert_fill_body, c, c%err, status, msg, msglen)
  end subroutine remd_convert_fill_c

  subroutine remd_convert_fill_body(ctx) bind(C, name="remd_convert_fill_body")
    implicit none
    type(c_ptr), value :: ctx
    type(t_remd_fill_ctx), pointer :: c
    character(kind=c_char), pointer :: dcd(:), rem(:)
    integer(c_int), pointer :: ids(:), selected(:), counts(:)
    type(c_ptr), pointer :: coords(:), boxes(:)
    integer :: dummy_id(1)

    call c_f_pointer(ctx, c)
    call c_f_pointer(c%dcd_ptr, dcd, [c%nreplica * c%name_len])
    call c_f_pointer(c%rem_ptr, rem, [c%nreplica * c%name_len])
    call c_f_pointer(c%selected_ptr, selected, [c%n_selected])
    call c_f_pointer(c%counts_ptr, counts, [c%n_out])
    call c_f_pointer(c%coords_ptr, coords, [c%n_out])
    call c_f_pointer(c%box_ptr, boxes, [c%n_out])
    dummy_id(1) = 1
    if (c%n_ids > 0) then
      call c_f_pointer(c%ids_ptr, ids, [c%n_ids])
      call remd_fill_frames(dcd, rem, c%nreplica, c%name_len, &
                            c%convert_type, ids, c%n_ids, &
                            c%nsteps, c%exchange_period, &
                            c%crdout_period, c%trjout_period, &
                            selected, c%n_selected, c%natom, &
                            counts, c%n_out, coords, boxes)
    else
      call remd_fill_frames(dcd, rem, c%nreplica, c%name_len, &
                            c%convert_type, dummy_id, 0, &
                            c%nsteps, c%exchange_period, &
                            c%crdout_period, c%trjout_period, &
                            selected, c%n_selected, c%natom, &
                            counts, c%n_out, coords, boxes)
    end if
  end subroutine remd_convert_fill_body

end module remd_convert_c_mod
