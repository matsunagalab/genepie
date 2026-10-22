!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   source_args_mod
!> @brief   Trajectory-source arguments shared by the bind(C) wrappers
!! @authors Claude Code
!
!  A Python STrajectories is either an in-memory (coords, pbc_boxes) pair or
!  the metadata of a lazily read DCD file. Wrappers that accept both receive
!  the same group of C arguments; this module packs them into t_source_args
!  and turns them into an s_trj_source with the validation the lazy RMSD/Rg
!  wrappers perform.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module source_args_mod
  use, intrinsic :: iso_c_binding
  use s_trajectories_c_mod
  use trj_source_mod
  use error_mod
  use string_mod
  use constants_mod
  implicit none
  private

  !> Trajectory-source arguments of a bind(C) entry point (plain copies).
  type, public :: t_source_args
    type(c_ptr) :: trajes_ptr = c_null_ptr       ! s_trajectories_c (memory)
    integer     :: is_lazy = 0                   ! 1: read the DCD lazily
    character(MaxFilename) :: filename_f = ''
    integer     :: trj_type = 0
    integer     :: dcd_natom_expected = 0
    type(c_ptr) :: source_selection_ptr = c_null_ptr
    integer     :: n_source_selection = 0
    integer     :: ana_period = 1
  end type t_source_args

  public :: pack_source_args
  public :: init_source_from_args

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pack_source_args
  !> @brief        Copy the trajectory-source C arguments into t_source_args
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pack_source_args(args, trajes_ptr, is_lazy, dcd_filename, &
                              filename_len, trj_type, dcd_natom_expected, &
                              source_selection_ptr, n_source_selection, &
                              ana_period)
    implicit none

    type(t_source_args),    intent(inout) :: args
    type(c_ptr),            intent(in)    :: trajes_ptr
    integer(c_int),         intent(in)    :: is_lazy
    character(kind=c_char), intent(in)    :: dcd_filename(*)
    integer(c_int),         intent(in)    :: filename_len
    integer(c_int),         intent(in)    :: trj_type
    integer(c_int),         intent(in)    :: dcd_natom_expected
    type(c_ptr),            intent(in)    :: source_selection_ptr
    integer(c_int),         intent(in)    :: n_source_selection
    integer(c_int),         intent(in)    :: ana_period

    args%trajes_ptr           = trajes_ptr
    args%is_lazy              = is_lazy
    args%filename_f           = ''
    if (is_lazy /= 0 .and. filename_len > 0) &
      call c_filename_to_fortran(dcd_filename, filename_len, args%filename_f)
    args%trj_type             = trj_type
    args%dcd_natom_expected   = dcd_natom_expected
    args%source_selection_ptr = source_selection_ptr
    args%n_source_selection   = n_source_selection
    args%ana_period           = ana_period

  end subroutine pack_source_args

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_source_from_args
  !> @brief        Initialise an s_trj_source (memory or lazy DCD) from the
  !!               packed arguments; failures are recorded in err only
  !! @authors      Claude Code
  !! @param[in]    args   : packed trajectory-source arguments
  !! @param[in]    natom  : atoms of the molecule the frames must match
  !! @param[in]    who    : wrapper name used in error messages
  !! @param[inout] source : trajectory source (released by the caller)
  !! @param[inout] err    : error state
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_source_from_args(args, natom, who, source, err)
    implicit none

    type(t_source_args), intent(in)    :: args
    integer,             intent(in)    :: natom
    character(*),        intent(in)    :: who
    type(s_trj_source),  intent(inout) :: source
    type(s_error),       intent(inout) :: err

    type(s_trajectories_c), pointer :: trajes
    integer, pointer :: source_selection_f(:)
    integer :: init_status

    if (args%ana_period <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     who // ': ana_period must be positive')
      return
    end if

    if (args%is_lazy /= 0) then

      if (args%n_source_selection /= natom .or. &
          .not. c_associated(args%source_selection_ptr)) then
        call error_set(err, ERROR_INVALID_PARAM, &
                       who // ': invalid source selection')
        return
      end if
      call c_f_pointer(args%source_selection_ptr, source_selection_f, &
                       [args%n_source_selection])
      if (any(source_selection_f < 1) .or. &
          any(source_selection_f > args%dcd_natom_expected)) then
        call error_set(err, ERROR_INVALID_PARAM, &
                       who // ': source selection out of range')
        return
      end if

      call init_source_lazy_dcd(source, trim(args%filename_f), &
                                args%trj_type, args%ana_period, &
                                source_selection_f, &
                                args%n_source_selection, init_status)
      if (init_status /= 0) then
        call error_set(err, init_status, &
                       who // ': unable to initialize DCD source')
        return
      end if

      if (source%dcd_natom /= args%dcd_natom_expected) then
        call error_set(err, ERROR_ATOM_COUNT, &
                       who // ': atom count mismatch')
        return
      end if

    else

      if (.not. c_associated(args%trajes_ptr)) then
        call error_set(err, ERROR_INVALID_PARAM, &
                       who // ': trajectory pointer is null')
        return
      end if
      call c_f_pointer(args%trajes_ptr, trajes)
      if (trajes%natom /= natom) then
        call error_set(err, ERROR_ATOM_COUNT, &
                       who // ': trajectory and molecule atom counts differ')
        return
      end if
      if (trajes%nframe <= 0 .or. .not. c_associated(trajes%coords)) then
        call error_set(err, ERROR_NO_DATA, &
                       who // ': trajectory has no frames')
        return
      end if

      call init_source_memory(source, trajes%coords, trajes%pbc_boxes, &
                              trajes%natom, trajes%nframe, args%ana_period)

    end if

  end subroutine init_source_from_args

end module source_args_mod
