module mesh_test

  ! Tests for mesh module

  use kinds_module
  use mpi_module
  use fruit
  use fson
  use mesh_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

public :: test_mesh_init
public :: test_mesh_order_vector

PetscReal, parameter :: tol = 1.e-6_dp

contains

!------------------------------------------------------------------------

  subroutine test_mesh_init

    ! Mesh init test

    use eos_module, only: max_primary_variable_name_length
    use dm_utils_module, only: section_offset
    use boundary_module, only: open_boundary_label_name
    use fson_mpi_module
    use cell_module
    use face_module

    character(max_mesh_filename_length) :: filename = "data/mesh/test_init.json"
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    character(max_primary_variable_name_length), allocatable :: primary_variable_names(:)
    Vec :: x
    type(face_type) :: face
    PetscInt :: global_solution_dof, num_primary
    PetscInt :: dim, facedof
    DM :: dm_face
    PetscSection :: section
    DMLabel :: ghost_label
    PetscReal, pointer :: fg(:)
    PetscInt :: f, offset, fstart, fend, ghost_face
    PetscReal :: dist(2)
    PetscErrorCode :: ierr
    character(len = 24) :: msg
    PetscInt, parameter :: expected_dim = 3, num_cells = 3, num_faces = 16
    PetscReal, parameter :: face_area = 200._dp
    PetscReal, parameter :: face_distance(2, 19:20) = &
         reshape([5._dp, 10._dp, 10._dp, 15._dp], [2,2])
    PetscReal, parameter :: face_centroid(3, 19:20) = &
         reshape([5._dp, 10._dp, 50._dp, 5._dp, 10._dp, 30._dp], [3,2])
    
    primary_variable_names = [character(max_primary_variable_name_length) :: &
         "Pressure", "Temperature"]

    json => fson_parse_mpi(filename)
    call mesh%init(json)
    call fson_destroy_mpi(json)

    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(primary_variable_names)

    call DMGetDimension(mesh%dm, dim, ierr); CHKERRQ(ierr)
    if (mpi%rank == mpi%output_rank) then
       call assert_equals(expected_dim, dim, "mesh dimension")
    end if

    call DMGetGlobalVector(mesh%dm, x, ierr); CHKERRQ(ierr)
    call VecGetSize(x, global_solution_dof, ierr); CHKERRQ(ierr)
    if (mpi%rank == mpi%output_rank) then
       num_primary = size(primary_variable_names)
       call assert_equals(num_cells * num_primary, global_solution_dof, "global solution dof")
    end if
    call DMRestoreGlobalVector(mesh%dm, x, ierr); CHKERRQ(ierr)

    call face%init()
    facedof = face%dof()

    ! Check face geometry:
    call VecGetDM(mesh%face_geom, dm_face, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(mesh%face_geom, fg, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(dm_face, section, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(mesh%dm, 1, fstart, fend, ierr); CHKERRQ(ierr)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    do f = fstart, fend-1
       call DMLabelGetValue(ghost_label, f, ghost_face, ierr); CHKERRQ(ierr)
       if (ghost_face < 0) then
          call section_offset(section, f, offset, ierr); CHKERRQ(ierr)
          call face%assign(fg, offset)
          write(msg, '(a, i2)') 'face area ', f
          call assert_equals(face_area, face%area, tol, msg)
          dist = face%distance
          if (face%normal(3) > 0._dp) then
             dist = [dist(2), dist(1)]
          end if
          write(msg, '(a, i2)') 'face distance ', f
          call assert_equals(0._dp, norm2(dist - face_distance(:,f)), tol, msg)
          write(msg, '(a, i2)') 'face centroid ', f
          call assert_equals(0._dp, norm2(face%centroid - face_centroid(:,f)), tol, msg)
       end if
    end do
    call face%destroy()
    call VecRestoreArrayF90(mesh%face_geom, fg, ierr)

    call mesh%destroy()
    deallocate(primary_variable_names)

  end subroutine test_mesh_init

!------------------------------------------------------------------------

  subroutine test_mesh_order_vector

    use eos_module, only: max_primary_variable_name_length
    use boundary_module, only: open_boundary_label_name
    use fson_mpi_module
    use dm_utils_module, only: global_vec_range_start, &
         global_section_offset

    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    character(max_primary_variable_name_length), allocatable :: &
         primary_variable_names(:)
    Vec :: v
    PetscViewer :: viewer
    PetscErrorCode :: ierr
    PetscInt :: c, ghost, order, offset, range_start, diff
    DMLabel :: ghost_label, order_label
    PetscReal, pointer :: v_array(:)
    PetscSection :: section
    character(len = 32) :: msg

    primary_variable_names = [character(max_primary_variable_name_length) :: &
         "Pressure"]

    json => fson_parse_mpi(str = '{"mesh": "data/mesh/col100.exo"}')
    call mesh%init(json)
    call fson_destroy_mpi(json)

    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(primary_variable_names)

    call DMGetGlobalVector(mesh%dm, v, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(v, "vec", ierr); CHKERRQ(ierr)

    call PetscViewerHDF5Open(mpi%comm, "data/mesh/vec100.h5", FILE_MODE_READ, &
         viewer, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5PushGroup(viewer, "/", ierr); CHKERRQ(ierr)
    call VecLoad(v, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

    call mesh%order_vector(v)

    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call DMGetLabel(mesh%dm, cell_order_label_name, order_label, ierr)
    call VecGetArrayF90(v, v_array, ierr); CHKERRQ(ierr)
    call DMGetDefaultGlobalSection(mesh%dm, section, ierr); CHKERRQ(ierr)
    call global_vec_range_start(v, range_start)
    diff = 0
    do c = mesh%start_cell, mesh%end_interior_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          call DMLabelGetValue(order_label, c, order, ierr); CHKERRQ(ierr)
          call global_section_offset(section, c, range_start, offset, ierr)
          diff = diff + abs(order - nint(v_array(offset)))
       end if
    end do
    call VecGetArrayF90(v, v_array, ierr); CHKERRQ(ierr)
    call DMRestoreGlobalVector(mesh%dm, v, ierr); CHKERRQ(ierr)

    write(msg, '(a, i3)') 'order difference, rank ', mpi%rank
    call assert_equals(0, diff, msg)

    call mesh%destroy()
    deallocate(primary_variable_names)

  end subroutine test_mesh_order_vector

!------------------------------------------------------------------------

end module mesh_test
