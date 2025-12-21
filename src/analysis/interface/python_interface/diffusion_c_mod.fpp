!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  ma_main
!! @brief   MSD analysis
!! @authors Donatas Surblys (DS), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module diffusion_c_mod

  use, intrinsic :: iso_c_binding
  use constants_mod
  use diffusion_impl_mod
  use da_setup_mod
  use da_control_mod
  use da_option_str_mod
  use input_str_mod
  use output_str_mod
  use error_mod
  use string_mod
  use messages_mod

  implicit none

contains

  subroutine diffusion_analysis_c(c_msd, msd_dim1, msd_dim2, &
          ctrl_text, ctrl_len, &
          out_data, status, msg, msglen) &
        bind(C, name="diffusion_analysis_c")
    use conv_f_c_util
    implicit none
    type(c_ptr), intent(in) :: c_msd
    integer(c_int), intent(in) :: msd_dim1
    integer(c_int), intent(in) :: msd_dim2
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out) :: out_data
    integer(c_int),          intent(out) :: status
    character(kind=c_char),  intent(out) :: msg(*)
    integer(c_int),          value       :: msglen

    real(wp), pointer :: f_msd(:,:)
    real(wp), pointer :: f_out_data(:,:)

    type(s_error) :: err

    call error_init(err)
    call C_F_POINTER(c_msd, f_msd, [msd_dim2, msd_dim1])
    allocate(f_out_data( &
        (size(f_msd, dim=1) - 1) * 2 + 1, &
        size(f_msd, dim=2)))
    call diffusion_analysis_main(f_msd, ctrl_text, ctrl_len, f_out_data, err)
    if (error_has(err)) then
      call error_to_c(err, status, msg, msglen)
      return
    end if
    status = 0
    if (msglen > 0) msg(1) = c_null_char

    out_data = f2c_double_array(f_out_data)
  end subroutine diffusion_analysis_c

  subroutine diffusion_analysis_main(msd_data, ctrl_text, ctrl_len, out_data, err)
    use, intrinsic :: iso_c_binding
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer, intent(in) :: ctrl_len
    real(wp), intent(in) :: msd_data(:,:)
    real(wp), intent(out) :: out_data(:,:)
    type(s_error), intent(inout) :: err
    ! local variables
    type(s_ctrl_data)      :: ctrl_data
    type(s_input)          :: input
    type(s_output)         :: output
    type(s_option)         :: option


    ! [Step1] Read control file
    !
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '

    call control_from_string(ctrl_text, ctrl_len, ctrl_data)


    ! [Step2] Set relevant variables and structures
    !
    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '

    call setup(ctrl_data, input, output, option)


    ! [Step3] Analysis of mean square dispacement
    !
    write(MsgOut,'(A)') '[STEP3] Analysis of mean square dispacement'
    call analyze(msd_data, input, option, out_data, err)


    ! [Step4] Deallocate memory
    !
    write(MsgOut,'(A)') '[STEP5] Deallocate memory'
    write(MsgOut,'(A)') ' '

    call dealloc_option(option)

  end subroutine diffusion_analysis_main

end module diffusion_c_mod
