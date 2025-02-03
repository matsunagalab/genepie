module s_trajectories_c_mod
  use, intrinsic :: iso_c_binding
  use trajectory_str_mod
  implicit none
  private

  ! multi frame s_trajectory
  type, public, bind(C) :: s_trajectories_c
    type(c_ptr) :: coords      ! doubles[nframe][natom][3]
    type(c_ptr) :: pbc_boxes   ! double[nframe][3][3]
    integer(c_int) :: nframe
    integer(c_int) :: natom
  end type s_trajectories_c

  public :: get_frame
  public :: set_frame
  public :: init_empty_s_trajectories_c
  public :: deallocate_s_trajectories_c
  public :: deallocate_s_trajectories_c_array

contains

  subroutine init_empty_s_trajectories_c(trajs_fort, natom, nframe) &
      bind(C, name="init_empty_s_trajectories_c")
    implicit none
    type(s_trajectories_c), intent(inout) :: trajs_fort
    integer(c_int), intent(in) :: natom
    integer(c_int), intent(in) :: nframe
    real(c_double), pointer :: coords(:,:,:), pbc_boxes(:,:,:)
    integer(c_int) :: i

    allocate(coords(3, natom, nframe))
    allocate(pbc_boxes(3, 3, nframe))

    pbc_boxes(:,:,:) = 0.0_C_double
    do i = 1, nframe
      pbc_boxes(1,1,i) = 1.0_C_double
      pbc_boxes(2,2,i) = 1.0_C_double
      pbc_boxes(3,3,i) = 1.0_C_double
    end do
    trajs_fort%coords = C_LOC(coords)
    trajs_fort%pbc_boxes = C_LOC(pbc_boxes)
    trajs_fort%nframe = nframe
    trajs_fort%natom = natom
  end subroutine init_empty_s_trajectories_c

  subroutine deallocate_s_trajectories_c(trajs_fort) &
      bind(C, name="deallocate_s_trajectories_c")
    implicit none
    type(s_trajectories_c), intent(inout) :: trajs_fort
    real(c_double), pointer :: coords(:,:,:), boxes(:,:,:)
    call C_F_POINTER(trajs_fort%coords, coords, [3, trajs_fort%natom, trajs_fort%nframe])
    deallocate(coords)
    call C_F_POINTER(trajs_fort%pbc_boxes, boxes, [3, 3, trajs_fort%nframe])
    deallocate(boxes)
  end subroutine deallocate_s_trajectories_c

  subroutine deallocate_s_trajectories_c_array(trajs_array, len) &
      bind(C, name="deallocate_s_trajectories_c_array")
    implicit none
    type(c_ptr), intent(inout) :: trajs_array
    integer(c_int), intent(in) :: len
    type(s_trajectories_c), pointer :: buf(:)
    integer :: i

    call C_F_POINTER(trajs_array, buf, [len])
    do i = 1, len
      call deallocate_s_trajectories_c(buf(i))
    end do
    deallocate(buf)
  end subroutine deallocate_s_trajectories_c_array

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

  subroutine set_frame(out_traj, src_traj, frame_idx)
    type(s_trajectories_c), intent(inout) :: out_traj
    type(s_trajectory), intent(in) :: src_traj
    integer, intent(in) :: frame_idx

    real(C_double), pointer :: coords(:,:,:), boxes(:,:,:)

    call C_F_POINTER(out_traj%coords, coords, [3, out_traj%natom, out_traj%nframe])
    call C_F_POINTER(out_traj%pbc_boxes, boxes, [3, 3, out_traj%nframe])
    coords(:, :, frame_idx) = src_traj%coord
    boxes(:, :, frame_idx) = src_traj%pbc_box
  end subroutine
end module s_trajectories_c_mod
