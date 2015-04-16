module mesh_test

  ! Tests for mesh module

  use kinds_module
  use mpi_module
  use fruit
  use fson
  use mesh_module

  implicit none
  private

#include <petsc-finclude/petsc.h90>

public :: test_init

PetscReal, parameter :: tol = 1.e-6_dp

contains

!------------------------------------------------------------------------

  subroutine test_init

    ! Mesh init test

    use eos_module, only: max_primary_variable_name_length
    use cell_module
    use face_module

    character(max_mesh_filename_length) :: filename = "data/mesh_test_init.json"
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    character(max_primary_variable_name_length), allocatable :: primary_variable_names(:)
    Vec :: x
    type(cell_type) :: cell
    type(face_type) :: face
    PetscInt :: global_solution_dof, global_cell_dof, global_face_dof, num_primary
    PetscInt :: dim, facedof
    DM :: dm_face
    PetscSection :: section
    DMLabel :: ghost_label
    PetscReal, pointer :: fg(:)
    PetscInt :: i, f, offset, fstart, fend, ghost_face, ghost_cell(2)
    PetscInt, pointer :: cells(:)
    PetscErrorCode :: ierr
    character(len = 24) :: msg
    PetscInt, parameter :: expected_dim = 3, num_cells = 3, num_faces = 16
    PetscReal, parameter :: area = 200._dp
    
    primary_variable_names = [character(max_primary_variable_name_length) :: &
         "Pressure", "Temperature"]

    if (mpi%rank == mpi%output_rank) then
       json => fson_parse(filename)
    end if

    call mesh%init(json, primary_variable_names)

    call DMGetDimension(mesh%dm, dim, ierr)
    if (mpi%rank == mpi%output_rank) then
       call assert_equals(expected_dim, dim, "mesh dimension")
    end if

    call DMCreateGlobalVector(mesh%dm, x, ierr)
    call VecGetSize(x, global_solution_dof, ierr)
    if (mpi%rank == mpi%output_rank) then
       num_primary = size(primary_variable_names)
       call assert_equals(num_cells * num_primary, global_solution_dof, "global solution dof")
    end if
    call VecDestroy(x, ierr)

    facedof = face%dof()
    call VecGetSize(mesh%face_geom, global_face_dof, ierr)
    if (mpi%rank == mpi%output_rank) then
       call assert_equals(num_faces * facedof, global_face_dof, "global face geometry dof")
    end if

    ! Check face geometry:
    call VecGetDM(mesh%face_geom, dm_face, ierr)
    call VecGetArrayF90(mesh%face_geom, fg, ierr)
    call DMGetDefaultSection(dm_face, section, ierr)
    call DMPlexGetHeightStratum(mesh%dm, 1, fstart, fend, ierr)
    call DMPlexGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    do f = fstart, fend-1
       call DMLabelGetValue(ghost_label, f, ghost_face, ierr); CHKERRQ(ierr)
       if (ghost_face < 0) then
          call DMPlexGetSupport(mesh%dm, f, cells, ierr)
          do i = 1, 2
             call DMLabelGetValue(ghost_label, cells(i), ghost_cell(i), ierr)
          end do
          if (all(ghost_cell < 0)) then
             call section_offset(section, f, offset, ierr); CHKERRQ(ierr)
             call face%assign(fg, offset)
             write(*,'(/a, i1, a, i2, a, i2, a, i2, 1x, i2, a, f6.2, a, f6.2)') 'rank ', mpi%rank, &
                  ' f ', f, ' offset ', offset, ' fstart, fend ', fstart, fend, &
                  ' area ', face%area, ' z ', face%centroid(3)
             write(msg, '(a, i2)') 'face area ', f
             call assert_equals(area, face%area, tol, msg)
          end if
       end if
    end do
    call VecRestoreArrayF90(mesh%face_geom, fg, ierr)

    call mesh%destroy()
    deallocate(primary_variable_names)

  end subroutine test_init

!------------------------------------------------------------------------

end module mesh_test
