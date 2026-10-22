!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pa_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pa_analyze_mod

  use pa_option_str_mod
  use fitting_mod
  use trj_source_mod
  use result_sink_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use input_str_mod
  use molecules_mod
  use molecules_str_mod
  use fileio_pdb_mod
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: analyze
  public  :: analyze_prjcrd_unified
  public  :: prjcrd_vcv_size
  private :: assign_mass
  private :: sum_project

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      prjcrd_vcv_size
  !> @brief        dimension of the PC space
  !!               (3 * analysis atoms for GLOBAL, 3 * fitting atoms for LOCAL)
  !! @authors      Claude Code
  !! @param[in]    fitting : fitting information
  !! @param[in]    option  : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function prjcrd_vcv_size(fitting, option) result(nvcv)

    ! formal arguments
    type(s_fitting),         intent(in)    :: fitting
    type(s_option),          intent(in)    :: option
    integer                                :: nvcv


    if (option%vcv_matrix == VcvMatrixGlobal) then
      nvcv = size(option%analysis_atom%idx) * 3
    else
      nvcv = size(fitting%fitting_atom%idx) * 3
    end if

    return

  end function prjcrd_vcv_size

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories (CLI driver)
  !! @authors      NT, Claude Code
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory file list information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] option     : option information
  !! @param[inout] input      : input information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, fitting, option, &
                     input, output)

    ! formal arguments
    type(s_molecule),         intent(inout) :: molecule
    type(s_trj_list), target, intent(inout) :: trj_list
    type(s_fitting),          intent(inout) :: fitting
    type(s_option),           intent(inout) :: option
    type(s_input),            intent(inout) :: input
    type(s_output),           intent(inout) :: output

    ! local variables
    type(s_pdb)              :: pdb_in
    type(s_trj_source)       :: source
    type(s_result_sink)      :: prj_sink
    real(wp)                 :: ttl_mass, rmsf
    integer                  :: natom, nsela, nself, nstru
    integer                  :: iatom, isel, nvcv
    integer                  :: val_in, vec_in
    integer                  :: i, j, alloc_stat

    real(wp),   allocatable  :: ave_coord (:,:)
    real(wp),   allocatable  :: eig_value (:)
    real(wp),   allocatable  :: eig_vector(:,:)


    if (option%check_only) &
      return

    natom = molecule%num_atoms
    nsela = size(option%analysis_atom%idx)
    nself = size(fitting%fitting_atom%idx)
    nvcv  = prjcrd_vcv_size(fitting, option)


    ! allocate memory
    !
    allocate(ave_coord (3,natom),    &
             eig_value (nvcv),       &
             eig_vector(nvcv, nvcv), stat=alloc_stat)

    if (alloc_stat /= 0) &
      call error_msg_alloc


    ! averaged coordinates: analysis atoms from pdb_avefile, fitting atoms
    ! from pdb_aftfile (both are exports of the same averaged structure)
    !
    ave_coord(1:3,1:natom) = 0.0_wp

    call input_pdb(input%pdb_avefile, pdb_in)

    if (nsela /= size(pdb_in%atom_name)) &
      call error_msg('Analyze> different count between pdb_avefile atom '// &
                     'and analysis_atom selection.')

    do isel = 1, nsela
      iatom = option%analysis_atom%idx(isel)
      ave_coord(1:3,iatom) = pdb_in%atom_coord(1:3,isel)
    end do

    call dealloc_pdb_all(pdb_in)

    call input_pdb(input%pdb_aftfile, pdb_in)

    if (nself /= size(pdb_in%atom_name)) &
      call error_msg('Analyze> different count between pdb_aftfile atom '// &
                     'and fitting_atom selection.')

    do isel = 1, nself
      iatom = fitting%fitting_atom%idx(isel)
      ave_coord(1:3,iatom) = pdb_in%atom_coord(1:3,isel)
    end do

    call dealloc_pdb_all(pdb_in)

    ! setup eigen-vector and eigen-val
    call open_file(val_in, input%valfile, IOFileInput)
    call open_file(vec_in, input%vecfile, IOFileInput)

    do i = 1, nvcv
      read(val_in,*,err=900,end=900) eig_value(i)
      do j = 1, nvcv
        read(vec_in,*,err=900,end=900) eig_vector(j,i)
      end do
    end do

    call close_file(val_in)
    call close_file(vec_in)


    ! projection (writes one row per structure to prjfile)
    !
    call init_sink_file_rows(prj_sink, output%prjfile, option%num_pca, &
                             '(i9,100(f10.5,1x))')

    call init_source_file(source, trj_list, natom)
    call analyze_prjcrd_unified(molecule, source, fitting, option, &
                                ave_coord, eig_vector, prj_sink, nstru)
    call finalize_source(source)
    call finalize_sink(prj_sink)


    ! rmsf of each mode (mass assigned by the core when mass_weight is on)
    !
    if (fitting%mass_weight) then
      ttl_mass = 0.0_wp
      do  isel = 1, nsela
        iatom = option%analysis_atom%idx(isel)
        ttl_mass = ttl_mass + molecule%mass(iatom)
      end do
    else
      ttl_mass = 1.0_wp
    end if

    write(MsgOut,'(A)') 'Analyze> rmsf (A)'

    do i = 1, option%num_pca
      rmsf = sqrt(eig_value(i) / ttl_mass)
      write(MsgOut,*) i, rmsf
    end do
    write(MsgOut,*)    ' '


    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [prjfile] ' // trim(output%prjfile)
    write(MsgOut,'(A)') '    Column 1: Snapshot index'
    write(MsgOut,'(A)') '    Column 2: Projected coordinates on the PC1 axis'
    write(MsgOut,'(A)') '       :                    :'
    write(MsgOut,'(A)') '    Column N: Projected coordinates on the PC(N-1) axis'
    write(MsgOut,'(A)') ''

    ! deallocate memory
    !
    deallocate(eig_vector, eig_value, ave_coord)

    return

900 call error_msg('Analyze> mismatch VEC/ VAL file from atom selection count')

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_prjcrd_unified
  !> @brief        projection of the fitted, mass-weighted fluctuations about
  !!               an averaged structure onto principal modes (shared by the
  !!               CLI and the Python interface)
  !! @authors      NT, Claude Code
  !! @param[inout] molecule   : molecule information (mass may be assigned)
  !! @param[inout] source     : trajectory source (file, memory or lazy DCD)
  !! @param[inout] fitting    : fitting information
  !! @param[in]    option     : option information (num_pca modes are used)
  !! @param[in]    ave_coord  : averaged coordinates (3, natom), not mass
  !!                            weighted; only the analysis and fitting
  !!                            atoms are read
  !! @param[in]    eig_vector : eigenvectors (nvcv, >= num_pca); column k is
  !!                            the k-th principal mode,
  !!                            nvcv = prjcrd_vcv_size(fitting, option)
  !! @param[inout] prj_sink   : receives num_pca projections per structure
  !! @param[out]   nstru_out  : number of structures analysed
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_prjcrd_unified(molecule, source, fitting, option, &
                                    ave_coord, eig_vector, prj_sink, nstru_out)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_trj_source),      intent(inout) :: source
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(in)    :: option
    real(wp),                intent(in)    :: ave_coord(:,:)
    real(wp),                intent(in)    :: eig_vector(:,:)
    type(s_result_sink),     intent(inout) :: prj_sink
    integer,                 intent(out)   :: nstru_out

    ! local variables
    type(s_trajectory)       :: trajectory
    real(wp)                 :: ttl_mass
    integer                  :: natom, nvcv, nsela, nstru
    integer                  :: iatom, isel, alloc_stat, frame_status

    real(wp),   allocatable  :: sqrt_mass (:)
    real(wp),   allocatable  :: ave_w     (:,:)
    real(wp),   allocatable  :: trj_coord (:,:)
    real(wp),   allocatable  :: vector(:)
    real(wp),   allocatable  :: project(:)


    natom = molecule%num_atoms
    nsela = size(option%analysis_atom%idx)
    nvcv  = prjcrd_vcv_size(fitting, option)

    if (option%num_pca < 1) &
      call error_msg('Analyze_Prjcrd_Unified> num_pca must be positive')
    if (option%num_pca > nvcv) &
      call error_msg('Analyze_Prjcrd_Unified> num_pca exceeds the number of modes')
    if (size(ave_coord,1) /= 3 .or. size(ave_coord,2) < natom) &
      call error_msg('Analyze_Prjcrd_Unified> ave_coord must be (3, natom)')
    if (size(eig_vector,1) < nvcv .or. size(eig_vector,2) < option%num_pca) &
      call error_msg('Analyze_Prjcrd_Unified> eig_vector is too small')

    allocate(sqrt_mass (natom),   &
             ave_w     (3,natom), &
             trj_coord (3,natom), &
             vector    (nvcv),    &
             project   (option%num_pca), stat=alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    ! setup mass
    !
    if (fitting%mass_weight) then
      call assign_mass(molecule)
      do iatom = 1, natom
        sqrt_mass(iatom) = sqrt(molecule%mass(iatom))
      end do
    else
      do iatom = 1, natom
        sqrt_mass(iatom) = 1.0_wp
      end do
    end if

    ! mass-weighted averaged coordinates
    !
    do iatom = 1, natom
      ave_w(1:3,iatom) = ave_coord(1:3,iatom) * sqrt_mass(iatom)
    end do

    ! setup total mass
    !
    if (fitting%mass_weight) then
      ttl_mass = 0.0_wp
      do  isel = 1, nsela
        iatom = option%analysis_atom%idx(isel)
        ttl_mass = ttl_mass + molecule%mass(iatom)
      end do
    else
      ttl_mass = 1.0_wp
    end if

    project(1:option%num_pca) = 0.0_wp

    ! analysis loop
    !
    nstru = 0
    call reset_source(source)

    do while (has_more_frames(source))

      call get_next_frame(source, trajectory, frame_status)
      if (frame_status /= 0) exit

      nstru = nstru + 1
      write(MsgOut,*) '      number of structures = ', nstru

      ! geometrical fitting
      !
      call run_fitting(fitting,             &
                       molecule%atom_coord, &
                       trajectory%coord,    &
                       trajectory%coord)

      ! mass-weighted fitting
      !
      do iatom = 1, natom
        trj_coord(1:3,iatom) = trajectory%coord(1:3,iatom)*sqrt_mass(iatom)
      end do

      call run_fitting(fitting,   &
                       ave_w,     &
                       trj_coord, &
                       trj_coord)

      if (option%vcv_matrix == VcvMatrixGlobal) then
        call sum_project(option%analysis_atom%idx, &
                         option%num_pca, ttl_mass, &
                         trj_coord, ave_w, eig_vector, vector, project)
      else
        call sum_project(fitting%fitting_atom%idx, &
                         option%num_pca, ttl_mass, &
                         trj_coord, ave_w, eig_vector, vector, project)
      end if

      call write_result_row(prj_sink, project)
      write(MsgOut, '(i9,100(f10.5,1x))') nstru, project(1:option%num_pca)

    end do

    nstru_out = nstru

    deallocate(project, vector, trj_coord, ave_w, sqrt_mass)

    return

  end subroutine analyze_prjcrd_unified

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    assign_mass
  !> @brief        assign mass
  !! @authors      NT
  !! @param[inout] molecule   : molecule information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine assign_mass(molecule)

    ! parameters
    real(wp),                parameter     :: MassH   =  1.008000_wp
    real(wp),                parameter     :: MassC   = 12.011000_wp
    real(wp),                parameter     :: MassN   = 14.007000_wp
    real(wp),                parameter     :: MassO   = 15.999000_wp
    real(wp),                parameter     :: MassS   = 32.060000_wp
    real(wp),                parameter     :: MassP   = 30.974000_wp
    real(wp),                parameter     :: MassMG  = 24.305000_wp
    real(wp),                parameter     :: MassZN  = 65.370000_wp

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule

    ! local variables
    integer                  :: i


    write(MsgOut,'(a)'),      'WARNING: atom mass is not assigned.'
    write(MsgOut,'(a)'),      '   uses default mass.'
    write(MsgOut,'(a,f9.6)')  '      1) H   : ', MassH
    write(MsgOut,'(a,f9.6)')  '      2) C   : ', MassC
    write(MsgOut,'(a,f9.6)')  '      3) N   : ', MassN
    write(MsgOut,'(a,f9.6)')  '      4) O   : ', MassO
    write(MsgOut,'(a,f9.6)')  '      5) S   : ', MassS
    write(MsgOut,'(a,f9.6)')  '      6) P   : ', MassP
    write(MsgOut,'(a,f9.6)')  '      7) MG  : ', MassMG
    write(MsgOut,'(a,f9.6)')  '      8) ZN  : ', MassZN


    do i = 1, molecule%num_atoms

      if (molecule%atom_name(i)(1:1) == 'H') then
        molecule%mass(i) = MassH
      else if (molecule%atom_name(i)(1:1) == 'C') then
        molecule%mass(i) = MassC
      else if (molecule%atom_name(i)(1:1) == 'N') then
        molecule%mass(i) = MassN
      else if (molecule%atom_name(i)(1:1) == 'O') then
        molecule%mass(i) = MassO
      else if (molecule%atom_name(i)(1:1) == 'S') then
        molecule%mass(i) = MassS
      else if (molecule%atom_name(i)(1:1) == 'P') then
        molecule%mass(i) = MassP
      else if (molecule%atom_name(i)(1:2) == 'MG') then
        molecule%mass(i) = MassMG
      else if (molecule%atom_name(i)(1:2) == 'ZN') then
        molecule%mass(i) = MassZN
      else
        write(MsgOut,'(a,a)') 'Assign_Mass> Unknown atom :', &
             molecule%atom_name(i)
      end if

    end do

    write(MsgOut,'(a)') ''

    return

  end subroutine assign_mass

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    sum_project
  !> @brief        sum project
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine sum_project(idx, num_pca, ttl_mass, trj_coord, ave_coord, &
                         eig_vector, vector, project)

    ! formal arguments
    integer,                 intent(in)    :: idx(:)
    integer,                 intent(in)    :: num_pca
    real(wp),                intent(in)    :: ttl_mass
    real(wp),                intent(in)    :: trj_coord(:,:)
    real(wp),                intent(in)    :: ave_coord(:,:)
    real(wp),                intent(in)    :: eig_vector(:,:)
    real(wp),                intent(inout) :: vector(:)
    real(wp),                intent(inout) :: project(:)

    ! local variables
    integer                  :: i, j, k, iatom, nsel, nvcv


    nsel = size(idx)
    nvcv  = nsel * 3

    do i = 1, nsel
      iatom = idx(i)
      vector(3*(i-1)+1) = trj_coord(1,iatom) - ave_coord(1,iatom)
      vector(3*(i-1)+2) = trj_coord(2,iatom) - ave_coord(2,iatom)
      vector(3*(i-1)+3) = trj_coord(3,iatom) - ave_coord(3,iatom)
    end do

    do i = 1, num_pca
      project(i) = 0.0_wp
      ! j = nvcv - i + 1
      do k = 1, nvcv
        project(i) = project(i) + vector(k) * eig_vector(k,i)
      end do
      project(i) = project(i) / sqrt(ttl_mass)
    end do

    return

  end subroutine sum_project

end module pa_analyze_mod
