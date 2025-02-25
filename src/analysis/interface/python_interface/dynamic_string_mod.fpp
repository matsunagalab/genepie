!> @file dynamic_string_mod.fpp
!> @brief This module provides utilities for dynamically appending strings.
!> It supports efficient memory management and newline handling.

module dynamic_string_mod
    implicit none

contains

    !> @brief Appends a new string to an existing dynamic string without adding a newline.
    !>
    !> This subroutine dynamically resizes the target string if necessary to accommodate
    !> the new content. The buffer size is doubled when needed to ensure efficiency.
    !>
    !> @param[inout] string The target dynamic string to which the new part will be appended.
    !>                      It must be allocatable and of type character(len=:).
    !> @param[in] new_part The new string to append to the target string.
    subroutine append_to_dynamic_string(string, new_part)
        character(len=:), allocatable, intent(inout) :: string
        character(len=*), intent(in) :: new_part      !< The new string to append
        integer :: required_length
        integer :: current_length
        character(len=:), allocatable :: temp_string  ! Temporary variable for reallocation

        ! Calculate the required length
        if (allocated(string)) then
            current_length = len(string)
            required_length = len_trim(string) + len_trim(new_part)
        else
            current_length = 0
            required_length = len_trim(new_part)
        end if

        ! Resize the buffer if necessary
        if (required_length > current_length) then
            do while (required_length > current_length)
                current_length = max(current_length * 2, 10)  ! Double the size (minimum size is 10)
            end do
            ! Use a temporary variable to preserve the existing content
            if (allocated(string)) then
                allocate(character(len=current_length) :: temp_string)
                temp_string = string
                deallocate(string)
                allocate(character(len=current_length) :: string)
                string = temp_string
                deallocate(temp_string)
            else
                allocate(character(len=current_length) :: string)
            end if
        end if

        ! Append the new part to the string
        if (allocated(string)) then
            string = string(:len_trim(string)) // trim(new_part)
        else
            string = trim(new_part)
        end if
    end subroutine append_to_dynamic_string

    !> @brief Appends a new string to an existing dynamic string with a newline at the end.
    !>
    !> This subroutine appends the new string followed by a newline character (CHAR(10)).
    !> Internally, it calls `append_to_dynamic_string` to handle the actual appending.
    !>
    !> @param[inout] string The target dynamic string to which the new part will be appended.
    !>                      It must be allocatable and of type character(len=:).
    !> @param[in] new_part The new string to append to the target string.
    subroutine append_to_dynamic_string_ln(string, new_part)
        character(len=:), allocatable, intent(inout) :: string
        character(len=*), intent(in) :: new_part      !< The new string to append
        character(len=len_trim(new_part) + 1) :: part_with_newline

        ! Add a newline character to the new part
        part_with_newline = trim(new_part) // CHAR(10)

        ! Call the base subroutine to append the modified string
        call append_to_dynamic_string(string, part_with_newline)
    end subroutine append_to_dynamic_string_ln

end module dynamic_string_mod
