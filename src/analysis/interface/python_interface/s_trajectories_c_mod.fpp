module s_trajectories_c_mod
  use, intrinsic :: iso_c_binding
  use trajectory_str_mod
  implicit none
  private

  ! multi frame s_trajectory
  type, public, bind(C) :: s_trajectories_c
    type(C_ptr) :: coords      ! doubles[nframe][natom][3]
    type(C_ptr) :: pbc_boxes   ! double[nframe][3][3]
    integer(C_int) :: nframe
    integer(C_int) :: natom
  end type s_trajectories_c

  public :: get_frame

contains

  ! Accessor for single frame
  subroutine get_frame(trajs_fort, frame_idx, traj_fort)
    type(s_trajectories_c), intent(in) :: trajs_fort
    integer, intent(in) :: frame_idx
    type(s_trajectory), intent(inout) :: traj_fort
    real(C_double), pointer :: coords(:,:,:), boxes(:,:,:)

    if (allocated(traj_fort%coord)) then
      if (.not. all(shape(traj_fort%coord) == [3, trajs_fort%natom])) then
        deallocate(traj_fort%coord)
        allocate(traj_fort%coord(3,trajs_fort%natom))
      end if
    else
      allocate(traj_fort%coord(3,trajs_fort%natom))
    end if
    call C_F_POINTER(trajs_fort%coords, coords, [3, trajs_fort%natom, trajs_fort%nframe])
    call C_F_POINTER(trajs_fort%pbc_boxes, boxes, [3, 3, trajs_fort%nframe])
    traj_fort%pbc_box = boxes(:, :, frame_idx)
    traj_fort%coord = coords(:, :, frame_idx)
  end subroutine
end module s_trajectories_c_mod
