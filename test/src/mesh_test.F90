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

PetscReal, parameter :: tol = 1.e-6_dp

contains

!------------------------------------------------------------------------

  subroutine test_mesh_init

    ! Mesh init test

    use dm_utils_module, only: section_offset
    use fson_mpi_module
    use cell_module
    use face_module

    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    character(len = 11), allocatable :: primary(:)
    Vec :: x
    type(face_type) :: face
    PetscInt :: global_solution_dof, num_primary
    PetscInt :: dim, facedof
    DM :: dm_face
    PetscSection :: section
    DMLabel :: ghost_label, cell_order_label
    PetscReal, pointer :: fg(:)
    PetscInt :: f, offset, fstart, fend, ghost_face, i, order(2), gf
    PetscInt, pointer :: cells(:)
    PetscReal :: dist(2)
    PetscErrorCode :: ierr
    character(len = 24) :: msg
    PetscInt, parameter :: expected_dim = 3, num_cells = 3, num_faces = 16
    PetscReal, parameter :: face_area = 200._dp
    PetscReal, parameter :: face_distance(2, 2) = &
         reshape([5._dp, 10._dp, 10._dp, 15._dp], [2,2])
    PetscReal, parameter :: face_centroid(3, 2) = &
         reshape([5._dp, 10._dp, 50._dp, 5._dp, 10._dp, 30._dp], [3,2])
    
    primary = ["Pressure   ", "Temperature"]

    json => fson_parse_mpi(str = '{"mesh": "data/mesh/block3.exo"}')
    call mesh%init(json)
    call fson_destroy_mpi(json)

    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(primary)

    call DMGetDimension(mesh%dm, dim, ierr); CHKERRQ(ierr)
    if (mpi%rank == mpi%output_rank) then
       call assert_equals(expected_dim, dim, "mesh dimension")
    end if

    call DMGetGlobalVector(mesh%dm, x, ierr); CHKERRQ(ierr)
    call VecGetSize(x, global_solution_dof, ierr); CHKERRQ(ierr)
    if (mpi%rank == mpi%output_rank) then
       num_primary = size(primary)
       call assert_equals(num_cells * num_primary, global_solution_dof, &
            "global solution dof")
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
    call DMGetLabel(mesh%dm, cell_order_label_name, cell_order_label, ierr)
    CHKERRQ(ierr)
    do f = fstart, fend - 1
       call DMLabelGetValue(ghost_label, f, ghost_face, ierr); CHKERRQ(ierr)
       if (ghost_face < 0) then
          call section_offset(section, f, offset, ierr); CHKERRQ(ierr)
          call face%assign_geometry(fg, offset)
          write(msg, '(a, i2)') 'face area ', f
          call assert_equals(face_area, face%area, tol, msg)
          dist = face%distance
          call DMPlexGetSupport(mesh%dm, f, cells, ierr); CHKERRQ(ierr)
          do i = 1, 2
             call DMLabelGetValue(cell_order_label, cells(i), order(i), ierr)
             CHKERRQ(ierr)
          end do
          if (all(order == [0,1])) then
             gf = 1
          else if (all(order == [1,2])) then
             gf = 2
          else
             gf = 0
          end if
          if (gf > 0) then
             write(msg, '(a, i2)') 'face distance ', f
             call assert_equals(0._dp, norm2(dist - face_distance(:, gf)), tol, msg)
             write(msg, '(a, i2)') 'face centroid ', f
             call assert_equals(0._dp, norm2(face%centroid - face_centroid(:, gf)), &
                  tol, msg)
          end if
       end if
    end do
    call face%destroy()
    call VecRestoreArrayF90(mesh%face_geom, fg, ierr); CHKERRQ(ierr)

    call mesh%destroy()
    deallocate(primary)

  end subroutine test_mesh_init

!------------------------------------------------------------------------

end module mesh_test
