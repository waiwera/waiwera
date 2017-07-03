module mesh_test

  ! Tests for mesh module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fruit
  use fson
  use mesh_module

  implicit none
  private

  public :: test_mesh_init, test_2d_cartesian_geometry, &
       test_2d_radial_geometry, test_mesh_face_permeability_direction

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
    Vec :: x
    type(face_type) :: face
    PetscInt :: global_solution_dof
    PetscInt, parameter :: dof = 2
    PetscInt :: dim
    DM :: dm_face
    PetscSection :: section
    DMLabel :: ghost_label, cell_order_label
    PetscReal, pointer, contiguous :: fg(:)
    PetscInt :: f, offset, fstart, fend, ghost_face, i, order(2), gf
    PetscInt, pointer :: cells(:)
    PetscReal :: dist(2)
    PetscErrorCode :: ierr
    PetscMPIInt :: rank
    character(len = 24) :: msg
    PetscInt, parameter :: expected_dim = 3, num_cells = 3, num_faces = 16
    PetscReal, parameter :: face_area = 200._dp
    PetscReal, parameter :: face_distance(2, 2) = &
         reshape([5._dp, 10._dp, 10._dp, 15._dp], [2,2])
    PetscReal, parameter :: face_centroid(3, 2) = &
         reshape([5._dp, 10._dp, 50._dp, 5._dp, 10._dp, 30._dp], [3,2])
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    
    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{"mesh": "data/mesh/block3.exo"}')
    call mesh%init(json)
    call fson_destroy_mpi(json)

    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(dof, gravity)

    call DMGetDimension(mesh%dm, dim, ierr); CHKERRQ(ierr)
    if (rank == 0) then
       call assert_equals(expected_dim, dim, "mesh dimension")
    end if

    call DMGetGlobalVector(mesh%dm, x, ierr); CHKERRQ(ierr)
    call VecGetSize(x, global_solution_dof, ierr); CHKERRQ(ierr)
    if (rank == 0) then
       call assert_equals(num_cells * dof, global_solution_dof, &
            "global solution dof")
    end if
    call DMRestoreGlobalVector(mesh%dm, x, ierr); CHKERRQ(ierr)

    call face%init()

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

  end subroutine test_mesh_init

!------------------------------------------------------------------------

  subroutine test_2d_cartesian_geometry
    ! Test 2D Cartesian volumes and areas

    use fson_mpi_module
    use cell_module
    use face_module
    use dm_utils_module, only: section_offset, local_vec_section

    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    PetscInt, parameter :: dof = 2
    PetscReal, pointer, contiguous :: cell_geom_array(:), face_geom_array(:)
    PetscSection :: cell_geom_section, face_geom_section
    PetscInt :: c, offset, f
    type(cell_type) :: cell
    type(face_type) :: face
    PetscErrorCode :: ierr
    PetscMPIInt :: rank
    character(50) :: msg
    PetscBool :: volumes_OK, areas_OK
    PetscReal, parameter :: expected_cell_vol = 62500._dp
    PetscReal, parameter :: expected_face_area = 2500._dp
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "data/mesh/2D.msh",' // &
         '"thickness": 100.}}')
    call mesh%init(json)
    call fson_destroy_mpi(json)
    call mesh%configure(dof, gravity)

    call local_vec_section(mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    volumes_OK = PETSC_TRUE

    do c = mesh%start_cell, mesh%end_cell - 1
       if (mesh%ghost_cell(c) < 0) then
          call section_offset(cell_geom_section, c, offset, &
               ierr); CHKERRQ(ierr)
          call cell%assign_geometry(cell_geom_array, offset)
          volumes_OK = (volumes_OK .and. &
               abs(cell%volume - expected_cell_vol) <= tol)
       end if
    end do

    write(msg, '(a, i3)') 'cell volumes on proc', rank
    call assert_true(volumes_OK, msg)

    call VecRestoreArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    call local_vec_section(mesh%face_geom, face_geom_section)
    call VecGetArrayReadF90(mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)
    call face%init()
    areas_OK = PETSC_TRUE

    do f = mesh%start_face, mesh%end_face - 1
       if (mesh%ghost_face(f) < 0) then
          call section_offset(face_geom_section, f, offset, &
               ierr); CHKERRQ(ierr)
          call face%assign_geometry(face_geom_array, offset)
          areas_OK = (areas_OK .and. &
               abs(face%area - expected_face_area) <= tol)
       end if
    end do

    write(msg, '(a, i3)') 'face areas on proc', rank
    call assert_true(areas_OK, msg)

    call VecRestoreArrayReadF90(mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)
    call face%destroy()

    call mesh%destroy()

  end subroutine test_2d_cartesian_geometry

!------------------------------------------------------------------------

  subroutine test_2d_radial_geometry
    ! Test 2D radial volumes and areas

    use fson_mpi_module
    use cell_module
    use face_module
    use dm_utils_module, only: section_offset, local_vec_section
    use utils_module, only: pi

    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    PetscInt, parameter :: dof = 2
    PetscReal, pointer, contiguous :: cell_geom_array(:), face_geom_array(:)
    PetscSection :: cell_geom_section, face_geom_section
    PetscInt :: c, offset, f
    type(cell_type) :: cell
    type(face_type) :: face
    PetscErrorCode :: ierr
    PetscMPIInt :: rank
    character(50) :: msg
    PetscBool :: volumes_OK, areas_OK
    PetscReal :: r, r1, r2, l, expected_cell_vol, expected_face_area
    PetscBool :: vertical
    PetscReal, parameter :: dr = 300._dp / 12, dy = 200._dp / 8
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "data/mesh/2D.msh",' // &
         '"radial": true}}')
    call mesh%init(json)
    call fson_destroy_mpi(json)
    call mesh%configure(dof, gravity)

    call local_vec_section(mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    volumes_OK = PETSC_TRUE

    do c = mesh%start_cell, mesh%end_cell - 1
       if (mesh%ghost_cell(c) < 0) then
          call section_offset(cell_geom_section, c, offset, &
               ierr); CHKERRQ(ierr)
          call cell%assign_geometry(cell_geom_array, offset)
          r = cell%centroid(1)
          r1 = r - 0.5_dp * dr
          r2 = r + 0.5_dp * dr
          expected_cell_vol = pi * (r2 * r2 - r1 * r1) * dy
          volumes_OK = (volumes_OK .and. &
               abs(cell%volume - expected_cell_vol) <= tol)
       end if
    end do

    write(msg, '(a, i3)') 'cell volumes on proc', rank
    call assert_true(volumes_OK, msg)

    call VecRestoreArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    call local_vec_section(mesh%face_geom, face_geom_section)
    call VecGetArrayReadF90(mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)
    call face%init()
    areas_OK = PETSC_TRUE

    do f = mesh%start_face, mesh%end_face - 1
       if (mesh%ghost_face(f) < 0) then
          call section_offset(face_geom_section, f, offset, &
               ierr); CHKERRQ(ierr)
          call face%assign_geometry(face_geom_array, offset)
          r = face%centroid(1)
          vertical = (abs(face%normal(2)) > tol)
          if (vertical) then
             l = dr
          else
             l = dy
          end if
          expected_face_area = 2._dp * pi * r * l
          areas_OK = (areas_OK .and. &
               abs(face%area - expected_face_area) <= tol)
       end if
    end do

    write(msg, '(a, i3)') 'face areas on proc', rank
    call assert_true(areas_OK, msg)

    call VecRestoreArrayReadF90(mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)
    call face%destroy()

    call mesh%destroy()

  end subroutine test_2d_radial_geometry

!------------------------------------------------------------------------

  subroutine test_mesh_face_permeability_direction
    ! Test face permeability direction

    use fson_mpi_module
    use dm_utils_module, only: local_vec_section, section_offset
    use face_module

    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    PetscInt, parameter :: dof = 2
    PetscInt :: f, offset
    PetscErrorCode :: ierr
    PetscSection :: face_geom_section
    type(face_type) :: face
    PetscReal :: dist
    PetscReal, pointer, contiguous :: face_geom_array(:)
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscReal, parameter :: position(3) = [1750._dp, 2000._dp, 400._dp]
    PetscReal, parameter :: tol = 1.e-6
    PetscInt, parameter :: expected_direction = 1

    json => fson_parse_mpi(str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo", ' // &
         '"faces": [' // &
         '{"cells": [16, 23], "permeability direction": 1}' // &
         ']}}')
    call mesh%init(json)

    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(dof, gravity)
    call mesh%override_face_properties(json)
    call fson_destroy_mpi(json)

    call local_vec_section(mesh%face_geom, face_geom_section)
    call VecGetArrayReadF90(mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)
    call face%init()

    do f = mesh%start_face, mesh%end_face - 1
       if (mesh%ghost_face(f) < 0) then
          call section_offset(face_geom_section, f, offset, ierr); CHKERRQ(ierr)
          call face%assign_geometry(face_geom_array, offset)
          dist = norm2(face%centroid - position)
          if (dist <= tol) then
             call assert_equals(expected_direction, &
                  nint(face%permeability_direction), &
                  'face [16, 23] direction')
          end if
       end if
    end do

    call face%destroy()
    call VecRestoreArrayReadF90(mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)
    call mesh%destroy()

  end subroutine test_mesh_face_permeability_direction

!------------------------------------------------------------------------

end module mesh_test
