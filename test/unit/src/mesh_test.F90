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
       test_2d_radial_geometry, test_mesh_face_permeability_direction, &
       test_setup_minc_dm, test_rock_assignment, test_minc_rock, &
       test_dm_cell_index

  PetscReal, parameter :: tol = 1.e-6_dp

contains

!------------------------------------------------------------------------

  subroutine test_mesh_init

    ! Mesh init test

    use cell_order_module, only: cell_order_label_name
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
    PetscErrorCode :: ierr, err
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
    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(dof, gravity, json, err = err)
    call fson_destroy_mpi(json)

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
    PetscErrorCode :: ierr, err
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
    call mesh%configure(dof, gravity, json, err = err)
    call fson_destroy_mpi(json)

    call local_vec_section(mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    volumes_OK = PETSC_TRUE

    do c = mesh%strata(0)%start, mesh%strata(0)%end - 1
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

    do f = mesh%strata(1)%start, mesh%strata(1)%end - 1
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
    PetscErrorCode :: ierr, err
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
    call mesh%configure(dof, gravity, json, err = err)
    call fson_destroy_mpi(json)

    call local_vec_section(mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    volumes_OK = PETSC_TRUE

    do c = mesh%strata(0)%start, mesh%strata(0)%end - 1
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

    do f = mesh%strata(1)%start, mesh%strata(1)%end - 1
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
    PetscErrorCode :: ierr, err
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
    call mesh%configure(dof, gravity, json, err = err)
    call mesh%override_face_properties(json)
    call fson_destroy_mpi(json)

    call local_vec_section(mesh%face_geom, face_geom_section)
    call VecGetArrayReadF90(mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)
    call face%init()

    do f = mesh%strata(1)%start, mesh%strata(1)%end - 1
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

  subroutine test_setup_minc_dm
    ! Test setup_minc_dm

    use fson_mpi_module

    type(fson_value), pointer :: json
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"zones": ["all"], "fracture": {"volume": 0.1}}}}')
    call minc_test('all', json, 1, 2 * 49, 1, [49, 49])
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo",' // &
         '  "zones": {"left": {"x": [0, 1500]}},' // &
         '  "minc": {"zones": ["left"], "matrix": {"volume": 0.9}}}}')
    call minc_test('partial', json, 1, 63, 1, [49, 14])
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo",' // &
         '  "zones": {"left": {"x": [0, 1500]}, "right": {"-": "left"}},' // &
         '  "minc": [{"zones": ["left"], "fracture": {"volume": 0.1}}, ' // &
         '           {"zones": ["right"], "fracture": {"volume": 0.1}, ' // &
         '            "matrix": {"volume": [0.3, 0.6]}}]}}')
    call minc_test('two-zone', json, 2, 133, 2, [49, 49, 35])
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo",' // &
         '  "zones": {"left": {"x": [0, 1500]}, ' // &
         '            "right corner": {"x": [2500, 4500], "y": [3000, 4500]}},' // &
         '  "minc": [{"zones": ["right corner"], "fracture": {"volume": 0.1}}, ' // &
         '   {"zones": ["left"], "matrix": {"volume": [0.3, 0.6]}}]}}')
    call minc_test('two-zone partial', json, 2, 83, 2, [49, 20, 14])
    call fson_destroy_mpi(json)

  contains

!........................................................................

    subroutine minc_test(name, json, expected_num_zones, expected_num_cells, &
         expected_max_level, expected_num_minc_level_cells)

      use minc_module, only: minc_level_label_name, minc_zone_label_name
      use cell_module
      use face_module
      use dm_utils_module, only: section_offset, local_vec_section
      use cell_order_module, only: cell_order_label_name

      character(*), intent(in) :: name
      type(fson_value), pointer, intent(in out) :: json
      PetscInt, intent(in) :: expected_num_zones, expected_num_cells
      PetscInt, intent(in) :: expected_max_level
      PetscInt, intent(in) :: expected_num_minc_level_cells(0: expected_max_level)
      ! Locals:
      type(mesh_type) :: mesh
      PetscInt, parameter :: dof = 2
      PetscInt :: num_cells, num_minc_zones, m, num_local, num, max_num_levels
      PetscErrorCode :: err
      PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
      character(48) :: str
      PetscInt :: iminc, num_minc_zone_cells, i, c
      IS :: minc_IS
      PetscInt, pointer :: minc_cells(:)
      Vec :: orig_cell_geom, orig_face_geom
      type(cell_type) :: orig_cell, cell
      type(face_type) :: face
      DMLabel :: ghost_label, minc_level_label, cell_order_label
      PetscInt, parameter :: nc = 1, np = 1 ! dummy values for cell init
      PetscSection :: orig_cell_section, cell_section, face_section
      PetscReal, pointer, contiguous :: orig_cell_geom_array(:), cell_geom_array(:)
      PetscReal, pointer, contiguous :: face_geom_array(:)
      PetscInt :: orig_cell_offset, cell_offset, ghost, cell_order, minc_cell_order
      PetscInt :: face_offset, cell_p, face_p, h
      PetscReal :: expected_vol, expected_area
      PetscInt :: ic(expected_max_level)
      PetscReal, parameter :: tol = 1.e-6_dp

      call mesh%init(json)
      call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
      call mesh%configure(dof, gravity, json, err = err)
      call assert_equals(0, err, name // ": minc config error")
      call fson_destroy_mpi(json)

      if (rank == 0) then
         call assert_true(mesh%has_minc, name // ": mesh has minc")
      end if
      num_minc_zones = size(mesh%minc)
      max_num_levels = maxval(mesh%minc%num_levels)
      if (rank == 0) then
         call assert_equals(expected_num_zones, &
              num_minc_zones, name // ": num minc zones")
         call assert_equals(expected_max_level, &
              max_num_levels, name // ": num minc levels")
      end if

      num_cells = total_interior_cell_count(mesh)
      if (rank == 0) then
         call assert_equals(expected_num_cells, &
              num_cells, name  // ": num cells")
      end if

      do m = 0, expected_max_level
         call DMGetStratumSize(mesh%minc_dm, minc_level_label_name, m, &
               num_local, ierr); CHKERRQ(ierr)
         call MPI_reduce(num_local, num, 1, MPI_INTEGER, MPI_SUM, &
              0, PETSC_COMM_WORLD, ierr)
         if (rank == 0) then
            write(str, '(a, i1)') ": num minc points, level ", m
            call assert_equals(expected_num_minc_level_cells(m), &
                 num, name  // str)
         end if
      end do

      ! Test geometry:
      call DMPlexComputeGeometryFVM(mesh%original_dm, orig_cell_geom, &
           orig_face_geom, ierr); CHKERRQ(ierr)
      call orig_cell%init(nc, np)
      call cell%init(nc, np)
      call face%init(nc, np)
      call local_vec_section(orig_cell_geom, orig_cell_section)
      call VecGetArrayReadF90(orig_cell_geom, orig_cell_geom_array, ierr)
      CHKERRQ(ierr)
      call local_vec_section(mesh%cell_geom, cell_section)
      call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
      CHKERRQ(ierr)
      call local_vec_section(mesh%face_geom, face_section)
      call VecGetArrayReadF90(mesh%face_geom, face_geom_array, ierr)
      CHKERRQ(ierr)
      call DMGetLabel(mesh%original_dm, "ghost", ghost_label, ierr)
      CHKERRQ(ierr)
      call DMGetLabel(mesh%minc_dm, minc_level_label_name, minc_level_label, ierr)
      CHKERRQ(ierr)
      call DMGetLabel(mesh%minc_dm, cell_order_label_name, &
         cell_order_label, ierr)
    CHKERRQ(ierr)
      ic = 0
      h = 0
      do iminc = 1, size(mesh%minc)
         associate(minc => mesh%minc(iminc))
           call DMGetStratumSize(mesh%original_dm, minc_zone_label_name, iminc, &
                num_minc_zone_cells, ierr); CHKERRQ(ierr)
           if (num_minc_zone_cells > 0) then
              call DMGetStratumIS(mesh%original_dm, minc_zone_label_name, &
                   iminc, minc_IS, ierr); CHKERRQ(ierr)
              call ISGetIndicesF90(minc_IS, minc_cells, ierr); CHKERRQ(ierr)
              do i = 1, num_minc_zone_cells
                 c = minc_cells(i)
                 call DMLabelGetValue(ghost_label, c, ghost, ierr)
                 if (ghost < 0) then
                    call DMLabelGetValue(cell_order_label, c, cell_order, ierr)
                    CHKERRQ(ierr)
                    call section_offset(orig_cell_section, c, orig_cell_offset, ierr)
                    CHKERRQ(ierr)
                    call orig_cell%assign_geometry(orig_cell_geom_array, orig_cell_offset)
                    cell_p = mesh%original_strata(h)%minc_point(c, 0)
                    call section_offset(cell_section, cell_p, cell_offset, ierr)
                    CHKERRQ(ierr)
                    call cell%assign_geometry(cell_geom_array, cell_offset)
                    expected_vol = orig_cell%volume * minc%volume(1)
                    write(str, '(a, a, i3, a)') name, ": fracture volume(", cell_order, ")"
                    call assert_equals(expected_vol, cell%volume, tol, str)
                    do m = 1, minc%num_levels
                       cell_p = mesh%original_strata(h)%minc_point(ic(m), m)
                       call section_offset(cell_section, cell_p, &
                            cell_offset, ierr); CHKERRQ(ierr)
                       call cell%assign_geometry(cell_geom_array, cell_offset)
                       expected_vol = orig_cell%volume * minc%volume(m + 1)
                       write(str, '(a, a, i3, a, i1, a)') name, ": minc volume(", &
                            cell_order, ", ", m, ")"
                       call assert_equals(expected_vol, cell%volume, tol, str)
                       face_p = mesh%original_strata(h + 1)%minc_point(ic(m), m)
                       call section_offset(face_section, face_p, &
                            face_offset, ierr); CHKERRQ(ierr)
                       call face%assign_geometry(face_geom_array, face_offset)
                       expected_area = orig_cell%volume * minc%connection_area(m)
                       write(str, '(a, a, i3, a, i1, a, i1, a)') name, &
                            ": minc area(", cell_order, ", ", m-1, ":", m, ")"
                       call assert_equals(expected_area, face%area, tol, str)
                       write(str, '(a, a, i3, a, i1, a, i1, a)') name, ": minc distance 1(", &
                            cell_order, ", ", m-1, ":", m, ")"
                       call assert_equals(minc%connection_distance(m), face%distance(1), &
                            tol, str)
                       write(str, '(a, a, i3, a, i1, a, i1, a)') name, ": minc distance 2(", &
                            cell_order, ", ", m-1, ":", m, ")"
                       call assert_equals(minc%connection_distance(m + 1), face%distance(2), &
                            tol, str)
                       call DMLabelGetValue(cell_order_label, cell_p, minc_cell_order, ierr)
                       CHKERRQ(ierr)
                       write(str, '(a, a, i3, a, i1, a)') name, ": minc cell order(", &
                            cell_order, ", ", m, ")"
                       call assert_equals(cell_order, minc_cell_order, str)
                       ic(m) = ic(m) + 1
                    end do
                 end if
              end do
           end if
         end associate
      end do
      call VecRestoreArrayReadF90(orig_cell_geom, orig_cell_geom_array, ierr)
      call VecRestoreArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
      call VecRestoreArrayReadF90(mesh%face_geom, face_geom_array, ierr)
      call VecDestroy(orig_cell_geom, ierr); CHKERRQ(ierr)
      call VecDestroy(orig_face_geom, ierr); CHKERRQ(ierr)
      call cell%destroy()
      call orig_cell%destroy()
      call face%destroy()
      call mesh%destroy()

    end subroutine minc_test

!........................................................................

    PetscInt function total_interior_cell_count(mesh) result(n)
      type(mesh_type), intent(in) :: mesh
      PetscInt :: c, n_local
      n_local = 0
      do c = mesh%strata(0)%start, mesh%strata(0)%end - 1
         if (mesh%ghost_cell(c) < 0) then
            n_local = n_local + 1
         end if
      end do
      call MPI_reduce(n_local, n, 1, MPI_INTEGER, MPI_SUM, &
         0, PETSC_COMM_WORLD, ierr)
    end function total_interior_cell_count

  end subroutine test_setup_minc_dm

!------------------------------------------------------------------------

  subroutine test_rock_assignment
    ! rock assignment

    use fson_mpi_module
    use rock_module
    use dm_utils_module, only: global_section_offset, global_vec_section
    use dictionary_module

    type(fson_value), pointer :: json
    character(:), allocatable :: json_str
    PetscErrorCode :: ierr
    PetscMPIInt :: rank
    PetscInt, parameter :: num_rocktypes = 2

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)

    json_str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo"}, ' // &
         '"rock": {"types": [' // &
         '  {"name": "rock1", "porosity": 0.1, ' // &
         '   "cells": [0,1,2,3,4,5,6,7,8,9,' // &
         '    10,11,12,13,14,15,16,17,18,19,20]},' // &
         '  {"name": "rock2", "porosity": 0.2, ' // &
         '    "cells": [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,' // &
         '    36,37,38,39,40,41,42,43,44,45,46,47,48]},' // &
         '  ]}}'

    call rock_test_case(json_str, [21, 28], "cells")

    json_str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo", ' // &
         '"zones": {' // &
         '"left_zone": {"x": [0, 3000]},' // &
         '"right_zone": {"-": "left_zone"}' // &
         '}}, ' // &
         '"rock": {"types": [' // &
         '  {"name": "rock1", "porosity": 0.1, ' // &
         '   "zones": ["left_zone"]}, ' // &
         '  {"name": "rock2", "porosity": 0.2, ' // &
         '    "zones": ["right_zone"]}' // &
         '  ]}}'

    call rock_test_case(json_str, [35, 14], "zones")

    json_str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo", ' // &
         '"zones": {' // &
         '"zone4": {"-": "zone3"},' // &
         '"zone3": {"+": "zone1", "-": "zone2"},' // &
         '"zone1": {"x": [0, 3000]},' // &
         '"zone2": {"x": [1500, 2500], "y": [1500, 2500]}' // &
         '}}, ' // &
         '"rock": {"types": [' // &
         '  {"name": "rock1", "porosity": 0.1, ' // &
         '   "zones": ["zone3"]}, ' // &
         '  {"name": "rock2", "porosity": 0.2, ' // &
         '    "zones": ["zone4"]}' // &
         '  ]}}'

    call rock_test_case(json_str, [31, 18], "zones")
  contains

    subroutine rock_test_case(json_str, expected_count, title)

      character(*), intent(in) :: json_str
      PetscInt, intent(in) :: expected_count(num_rocktypes)
      character(*), intent(in) :: title
      ! Locals:
      type(mesh_type) :: mesh
      PetscInt, parameter :: dof = 2
      Vec :: rock_vector
      type(dictionary_type) :: rock_dict
      type(rock_type) :: rock
      PetscReal, contiguous, pointer :: rock_array(:)
      PetscSection :: section
      DMLabel :: ghost_label
      PetscInt :: rock_range_start, start_cell, end_cell
      PetscInt :: c, ghost, offset, irock
      PetscErrorCode :: ierr, err
      PetscInt :: rock_count_local(num_rocktypes), rock_count(num_rocktypes)
      PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
      PetscReal, parameter :: tol = 1.e-6

      json => fson_parse_mpi(str = json_str)
      call mesh%init(json)

      call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
      call mesh%configure(dof, gravity, json, err = err)
      call rock_dict%init(owner = PETSC_TRUE)
      call setup_rock_vector(json, mesh%dm, rock_vector, rock_dict, &
           rock_range_start, mesh%ghost_cell, err = err)
      call assert_equals(0, err, "setup rock vector error")
      call fson_destroy_mpi(json)

      call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
      call global_vec_section(rock_vector, section)
      call rock%init()
      call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
      CHKERRQ(ierr)
      call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)

      rock_count_local = 0

      do c = start_cell, end_cell - 1
         call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
         if (ghost < 0) then
            call global_section_offset(section, c, rock_range_start, &
                 offset, ierr); CHKERRQ(ierr)
            call rock%assign(rock_array, offset)
            if (rock%porosity < 0.15_dp) then
               irock = 1
            else
               irock = 2
            end if
            rock_count_local(irock) = rock_count_local(irock) + 1
         end if
      end do

      call MPI_reduce(rock_count_local, rock_count, num_rocktypes, &
           MPI_INTEGER, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)

      if (rank == 0) then
         call assert_equals(expected_count, rock_count, num_rocktypes, &
              "rock counts: " // trim(title))
      end if

      call rock%destroy()
      call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)

      call VecDestroy(rock_vector, ierr); CHKERRQ(ierr)
      call mesh%destroy()
      call rock_dict%destroy()

    end subroutine rock_test_case

  end subroutine test_rock_assignment

!------------------------------------------------------------------------

  subroutine test_minc_rock
    ! MINC rock assignment

    use fson_mpi_module
    use dictionary_module
    use rock_module
    use minc_module, only: minc_level_label_name
    use dm_utils_module, only: global_section_offset, global_vec_section

    PetscMPIInt :: rank
    PetscErrorCode :: ierr
    character(:), allocatable :: json_str

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json_str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"zones": "all", ' // &
         '           "fracture": {"volume": 0.1, "planes": 3, ' // &
         '                        "spacing": 100, "rock": {"type": "fracture"}}, ' // &
         '           "matrix": {"volume": 0.9, "rock": {"type": "matrix"}}}},' // &
         ' "rock": {"types": [{"name": "original", "porosity": 0.1, "zones": "all"}, ' // &
         '                    {"name": "fracture", "porosity": 0.6}, ' // &
         '                    {"name": "matrix", "porosity": 0.02}]}}'

    call minc_rock_test_case(json_str, "case 1", 1, 0.6_dp, 0.02_dp)

    json_str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"zones": "all", ' // &
         '           "fracture": {"volume": 0.1, "planes": 3, ' // &
         '                        "spacing": 100, "rock": {"type": "fracture"}}, ' // &
         '           "matrix": {"volume": [0.3, 0.6], "rock": {"type": "matrix"}}}},' // &
         ' "rock": {"types": [{"name": "original", "porosity": 0.1, "zones": "all"}, ' // &
         '                    {"name": "fracture", "porosity": 0.6}, ' // &
         '                    {"name": "matrix"}]}}'

    call minc_rock_test_case(json_str, "case 2", 2, 0.6_dp, 2._dp / 45._dp)

  contains

    subroutine minc_rock_test_case(json_str, title, num_levels, &
         expected_fracture_porosity, expected_matrix_porosity)

      character(*), intent(in) :: json_str
      character(*), intent(in) :: title
      PetscInt, intent(in) :: num_levels
      PetscReal, intent(in) :: expected_fracture_porosity, expected_matrix_porosity
      ! Locals:
      type(fson_value), pointer :: json
      type(mesh_type) :: mesh
      PetscInt, parameter :: dof = 2
      Vec :: rock_vector
      type(dictionary_type) :: rock_dict
      PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
      PetscErrorCode :: err
      PetscInt :: rock_range_start, c, i, m, num_cells, offset
      PetscReal, contiguous, pointer :: rock_array(:)
      type(rock_type) :: rock
      PetscSection :: section
      IS :: minc_IS
      PetscInt, pointer :: minc_points(:)
      PetscReal :: expected_porosity(0 : num_levels)
      character(8) :: levelstr
      PetscReal, parameter :: tol = 1.e-6

      expected_porosity = expected_matrix_porosity
      expected_porosity(0) = expected_fracture_porosity

      json => fson_parse_mpi(str = json_str)
      call mesh%init(json)

      call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
      call mesh%configure(dof, gravity, json, err = err)
      call assert_equals(0, err, title // " mesh configure error")
      call rock_dict%init(owner = PETSC_TRUE)
      call setup_rock_vector(json, mesh%dm, rock_vector, rock_dict, &
           rock_range_start, mesh%ghost_cell, err = err)
      call assert_equals(0, err, title // " setup rock vector error")
      call mesh%setup_minc_rock_properties(json, rock_vector, &
           rock_dict, rock_range_start, err = err)
      call assert_equals(0, err, title // " setup MINC rock properties error")
      call fson_destroy_mpi(json)

      call VecGetArrayReadF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
      call global_vec_section(rock_vector, section)
      call rock%init()

      do m = 0, num_levels
         write(levelstr, '(a, i1)') ' level ', m
         call DMGetStratumSize(mesh%minc_dm, minc_level_label_name, m, &
              num_cells, ierr); CHKERRQ(ierr)
         if (num_cells > 0) then
            call DMGetStratumIS(mesh%minc_dm, minc_level_label_name, &
                 m, minc_IS, ierr); CHKERRQ(ierr)
            call ISGetIndicesF90(minc_IS, minc_points, ierr); CHKERRQ(ierr)
            do i = 1, size(minc_points)
               c = minc_points(i)
               if (mesh%ghost_cell(c) < 0) then
                  call global_section_offset(section, c, rock_range_start, &
                       offset, ierr); CHKERRQ(ierr)
                  call rock%assign(rock_array, offset)
                  call assert_equals(expected_porosity(m), rock%porosity, tol, &
                       title // levelstr)
               end if
            end do
            call ISRestoreIndicesF90(minc_IS, minc_points, ierr); CHKERRQ(ierr)
         end if
      end do

      call VecRestoreArrayReadF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
      call VecDestroy(rock_vector, ierr); CHKERRQ(ierr)
      call mesh%destroy()
      call rock_dict%destroy()
      call rock%destroy()

    end subroutine minc_rock_test_case

  end subroutine test_minc_rock

!------------------------------------------------------------------------

  subroutine test_dm_cell_index
    ! DM cell index

    use fson_mpi_module

    PetscMPIInt :: rank
    PetscInt :: i
    PetscErrorCode :: ierr
    character(:), allocatable :: json_str

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json_str = &
         '{"mesh": {"filename": "data/mesh/3x3grid.exo"}}'
    call cell_index_test_case(json_str, 'no boundary', [(i, i = 0, 8)])

    json_str = &
         '{"mesh": {"filename": "data/mesh/3x3grid.exo"}, ' // &
         ' "boundaries": [{"faces": {' // &
         ' "cells": [0, 3], "normal": [-1, 0, 0]' // &
         '}}]}'
    call cell_index_test_case(json_str, 'boundary', [(i, i = 0, 8)])

    json_str = &
         '{"mesh": {"filename": "data/mesh/3x3grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"zones": ["all"], "fracture": {"volume": 0.1}}}}'
    call cell_index_test_case(json_str, 'minc no boundary', &
         [[(i, i = 0, 8)], [(i, i = 0, 8)]])

    json_str = &
         '{"mesh": {"filename": "data/mesh/3x3grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"zones": ["all"], "fracture": {"volume": 0.1}}}, ' // &
         '  "boundaries": [{"faces": {' // &
         '    "cells": [0, 3], "normal": [-1, 0, 0]}}]' // &
         '}'
    call cell_index_test_case(json_str, 'minc boundary', &
         [[(i, i = 0, 8)], [(i, i = 0, 8)]])

  contains

    subroutine cell_index_test_case(json_str, title, expected_order)

      use dm_utils_module, only: global_vec_range_start, global_vec_section, &
           global_section_offset
      use cell_order_module, only: cell_order_label_name
      use IFC67_module
      use eos_w_module

      character(*), intent(in) :: json_str
      character(*), intent(in) :: title
      PetscInt, intent(in) :: expected_order(:)
      ! Locals:
      type(fson_value), pointer :: json
      type(mesh_type) :: mesh
      PetscInt, parameter :: dof = 1
      PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
      PetscInt :: c, start_cell, end_cell, ghost, order, i
      PetscInt :: range_start, offset, count
      Vec :: v, v0
      VecScatter :: scatter
      IS :: is0
      PetscReal, pointer :: varray(:)
      PetscInt, pointer :: isarray(:)
      PetscInt, allocatable :: idx(:)
      PetscSection :: section
      DMLabel :: ghost_label, order_label
      type(IFC67_type) :: thermo
      type(eos_w_type) :: eos
      PetscErrorCode :: err

      json => fson_parse_mpi(str = json_str)
      call mesh%init(json)
      call thermo%init()
      call eos%init(json, thermo)
      call mesh%setup_boundaries(json, eos)
      call mesh%configure(dof, gravity, json, err = err)
      call fson_destroy_mpi(json)
      
      call DMGetGlobalVector(mesh%dm, v, ierr); CHKERRQ(ierr)

      call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
      CHKERRQ(ierr)
      call global_vec_range_start(v, range_start)
      call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
      call DMGetLabel(mesh%dm, cell_order_label_name, order_label, ierr); CHKERRQ(ierr)
      call VecGetArrayF90(v, varray, ierr); CHKERRQ(ierr)
      call global_vec_section(v, section)

      count = 0
      do c = start_cell, end_cell - 1
         call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
         if (ghost < 0) then
            call global_section_offset(section, c, range_start, &
                 offset, ierr); CHKERRQ(ierr)
            call DMLabelGetValue(order_label, c, order, ierr); CHKERRQ(ierr)
            varray(offset) = dble(order)
            count = count + 1
         end if
      end do
      call VecRestoreArrayF90(v, varray, ierr); CHKERRQ(ierr)

      call VecScatterCreateToZero(v, scatter, v0, ierr); CHKERRQ(ierr)
      call VecScatterBegin(scatter, v, v0, INSERT_VALUES, SCATTER_FORWARD, ierr)
      call VecScatterEnd(scatter, v, v0, INSERT_VALUES, SCATTER_FORWARD, ierr)

      call ISAllGather(mesh%cell_index, is0, ierr); CHKERRQ(ierr)
      call ISGetIndicesF90(is0, isarray, ierr); CHKERRQ(ierr)
      call VecGetArrayF90(v0, varray, ierr); CHKERRQ(ierr)

      if (rank == 0) then
         associate(n => size(isarray))
           allocate(idx(n))
           idx = 0
           do i = 1, n
              idx(i) = int(varray(isarray(i) + 1))
           end do
           call assert_equals(expected_order, idx, n, title // ': cell order')
           deallocate(idx)
         end associate
      end if

      call VecRestoreArrayF90(v0, varray, ierr); CHKERRQ(ierr)
      call ISRestoreIndicesF90(is0, isarray, ierr); CHKERRQ(ierr)
      call ISDestroy(is0, ierr); CHKERRQ(ierr)

      call VecScatterDestroy(scatter, ierr); CHKERRQ(ierr)
      call VecDestroy(v0, ierr); CHKERRQ(ierr)
      call DMRestoreGlobalVector(mesh%dm, v, ierr); CHKERRQ(ierr)
      call mesh%destroy()
      call eos%destroy()
      call thermo%destroy()

    end subroutine cell_index_test_case

  end subroutine test_dm_cell_index

!------------------------------------------------------------------------

end module mesh_test
