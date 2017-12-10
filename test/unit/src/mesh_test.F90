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
       test_setup_minc_dm, test_minc_rock, &
       test_rock_assignment, test_cell_order

  PetscReal, parameter :: tol = 1.e-6_dp

contains

!------------------------------------------------------------------------

  subroutine test_mesh_init

    ! Mesh init test

    use dm_utils_module, only: section_offset, local_to_natural_cell_index
    use fson_mpi_module
    use cell_module
    use face_module
    use IAPWS_module
    use eos_we_module

    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(mesh_type) :: mesh
    Vec :: x
    type(face_type) :: face
    PetscInt :: global_solution_dof

    PetscInt :: dim
    DM :: dm_face
    PetscSection :: section
    DMLabel :: ghost_label
    PetscReal, pointer, contiguous :: fg(:)
    PetscInt :: f, offset, fstart, fend, ghost_face, order(2), gf
    PetscInt, pointer :: cells(:)
    PetscReal :: dist(2)
    PetscErrorCode :: ierr, err
    PetscMPIInt :: rank
    character(len = 24) :: msg
    PetscViewer :: viewer
    ISLocalToGlobalMapping :: l2g
    PetscInt, parameter :: expected_dim = 3, num_cells = 3, num_faces = 16
    PetscReal, parameter :: face_area = 200._dp
    PetscReal, parameter :: face_distance(2, 2) = &
         reshape([5._dp, 10._dp, 10._dp, 15._dp], [2,2])
    PetscReal, parameter :: face_centroid(3, 2) = &
         reshape([5._dp, 10._dp, 50._dp, 5._dp, 10._dp, 30._dp], [3,2])
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    
    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call thermo%init()
    call eos%init(json, thermo)
    viewer = PETSC_NULL_VIEWER

    json => fson_parse_mpi(str = '{"mesh": "data/mesh/block3.exo"}')
    call mesh%init(json)
    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(eos, gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)

    call DMGetDimension(mesh%dm, dim, ierr); CHKERRQ(ierr)
    if (rank == 0) then
       call assert_equals(expected_dim, dim, "mesh dimension")
    end if

    call DMGetGlobalVector(mesh%dm, x, ierr); CHKERRQ(ierr)
    call VecGetSize(x, global_solution_dof, ierr); CHKERRQ(ierr)
    if (rank == 0) then
       call assert_equals(num_cells * eos%num_primary_variables, &
            global_solution_dof, "global solution dof")
    end if
    call DMRestoreGlobalVector(mesh%dm, x, ierr); CHKERRQ(ierr)

    call face%init()
    call DMGetLocalToGlobalMapping(mesh%dm, l2g, ierr); CHKERRQ(ierr)

    ! Check face geometry:
    call VecGetDM(mesh%face_geom, dm_face, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(mesh%face_geom, fg, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(dm_face, section, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(mesh%dm, 1, fstart, fend, ierr); CHKERRQ(ierr)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
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
          order = local_to_natural_cell_index(mesh%cell_order, l2g, cells)
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
    call eos%destroy()
    call thermo%destroy()

  end subroutine test_mesh_init

!------------------------------------------------------------------------

  subroutine test_2d_cartesian_geometry
    ! Test 2D Cartesian volumes and areas

    use fson_mpi_module
    use cell_module
    use face_module
    use dm_utils_module, only: section_offset, local_vec_section
    use IAPWS_module
    use eos_we_module

    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(mesh_type) :: mesh
    PetscViewer :: viewer
    PetscReal, pointer, contiguous :: cell_geom_array(:), face_geom_array(:)
    PetscSection :: cell_geom_section, face_geom_section
    PetscInt :: c, offset, f, start_cell, end_cell, start_face, end_face
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
    call thermo%init()
    call eos%init(json, thermo)
    viewer = PETSC_NULL_VIEWER

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "data/mesh/2D.msh",' // &
         '"thickness": 100.}}')
    call mesh%init(json)
    call mesh%configure(eos, gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)

    call local_vec_section(mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    volumes_OK = PETSC_TRUE

    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)

    do c = start_cell, end_cell - 1
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

    call DMPlexGetHeightStratum(mesh%dm, 1, start_face, end_face, ierr)
    CHKERRQ(ierr)
    do f = start_face, end_face - 1
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
    call eos%destroy()
    call thermo%destroy()

  end subroutine test_2d_cartesian_geometry

!------------------------------------------------------------------------

  subroutine test_2d_radial_geometry
    ! Test 2D radial volumes and areas

    use fson_mpi_module
    use cell_module
    use face_module
    use dm_utils_module, only: section_offset, local_vec_section
    use utils_module, only: pi
    use IAPWS_module
    use eos_we_module

    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    PetscViewer :: viewer
    PetscReal, pointer, contiguous :: cell_geom_array(:), face_geom_array(:)
    PetscSection :: cell_geom_section, face_geom_section
    PetscInt :: c, offset, f, start_cell, end_cell, start_face, end_face
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
    call thermo%init()
    call eos%init(json, thermo)
    viewer = PETSC_NULL_VIEWER

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "data/mesh/2D.msh",' // &
         '"radial": true}}')
    call mesh%init(json)
    call mesh%configure(eos, gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)

    call local_vec_section(mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    volumes_OK = PETSC_TRUE
    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)

    do c = start_cell, end_cell - 1
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
    call DMPlexGetHeightStratum(mesh%dm, 1, start_face, end_face, ierr)
    CHKERRQ(ierr)

    do f = start_face, end_face - 1
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
    call eos%destroy()
    call thermo%destroy()

  end subroutine test_2d_radial_geometry

!------------------------------------------------------------------------

  subroutine test_mesh_face_permeability_direction
    ! Test face permeability direction

    use fson_mpi_module
    use dm_utils_module, only: local_vec_section, section_offset
    use face_module
    use IAPWS_module
    use eos_we_module

    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(mesh_type) :: mesh
    PetscViewer :: viewer
    PetscInt :: f, offset, start_face, end_face
    PetscErrorCode :: ierr, err
    PetscSection :: face_geom_section
    type(face_type) :: face
    PetscReal :: dist
    PetscReal, pointer, contiguous :: face_geom_array(:)
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscReal, parameter :: position(3) = [1750._dp, 2000._dp, 400._dp]
    PetscReal, parameter :: tol = 1.e-6
    PetscInt, parameter :: expected_direction = 1

    call thermo%init()
    call eos%init(json, thermo)
    viewer = PETSC_NULL_VIEWER

    json => fson_parse_mpi(str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo", ' // &
         '"faces": [' // &
         '{"cells": [16, 23], "permeability direction": 1}' // &
         ']}}')
    call mesh%init(json)

    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(eos, gravity, json, viewer = viewer, err = err)
    call mesh%override_face_properties(json)
    call fson_destroy_mpi(json)

    call local_vec_section(mesh%face_geom, face_geom_section)
    call VecGetArrayReadF90(mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)
    call face%init()
    call DMPlexGetHeightStratum(mesh%dm, 1, start_face, end_face, ierr)
    CHKERRQ(ierr)

    do f = start_face, end_face - 1
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
    call eos%destroy()
    call thermo%destroy()

  end subroutine test_mesh_face_permeability_direction

!------------------------------------------------------------------------

  subroutine test_setup_minc_dm
    ! Test setup_minc_dm

    use fson_mpi_module

    character(:), allocatable :: json_str
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json_str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"zones": ["all"], "fracture": {"volume": 0.1}}}}'
    call minc_test('all', json_str, 1, 2 * 49, 1, [49, 49])

    json_str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo",' // &
         '  "zones": {"left": {"x": [0, 1500]}},' // &
         '  "minc": {"zones": ["left"], "matrix": {"volume": 0.9}}}}'
    call minc_test('partial', json_str, 1, 63, 1, [49, 14])

    json_str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo",' // &
         '  "zones": {"left": {"x": [0, 1500]}, "right": {"-": "left"}},' // &
         '  "minc": [{"zones": ["left"], "fracture": {"volume": 0.1}}, ' // &
         '           {"zones": ["right"], "fracture": {"volume": 0.1}, ' // &
         '            "matrix": {"volume": [0.3, 0.6]}}]}}'
    call minc_test('two-zone', json_str, 2, 133, 2, [49, 49, 35])

    json_str = &
         '{"mesh": {"filename": "data/mesh/7x7grid.exo",' // &
         '  "zones": {"left": {"x": [0, 1500]}, ' // &
         '            "right corner": {"x": [2500, 4500], "y": [3000, 4500]}},' // &
         '  "minc": [{"zones": ["right corner"], "fracture": {"volume": 0.1}}, ' // &
         '   {"zones": ["left"], "matrix": {"volume": [0.3, 0.6]}}]}}'
    call minc_test('two-zone partial', json_str, 2, 83, 2, [49, 20, 14])

  contains

!........................................................................

    subroutine minc_test(name, json_str, expected_num_zones, expected_num_cells, &
         expected_max_level, expected_num_minc_level_cells)

      use minc_module, only: minc_level_label_name, minc_zone_label_name
      use cell_module
      use face_module
      use dm_utils_module, only: section_offset, local_vec_section
      use cell_order_module, only: cell_order_label_name
      use IAPWS_module
      use eos_we_module

      character(*), intent(in) :: name
      character(*), intent(in) :: json_str
      PetscInt, intent(in) :: expected_num_zones, expected_num_cells
      PetscInt, intent(in) :: expected_max_level
      PetscInt, intent(in) :: expected_num_minc_level_cells(0: expected_max_level)
      ! Locals:
      type(fson_value), pointer :: json, orig_json
      type(mesh_type) :: mesh, orig_mesh
      character(:), allocatable :: orig_json_str
      type(IAPWS_type) :: thermo
      type(eos_we_type) :: eos
      PetscViewer :: viewer
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
      PetscInt :: orig_cell_offset, cell_offset, ghost, cell_order
      PetscInt :: face_offset, cell_p, face_p, h
      PetscReal :: expected_vol, expected_area
      PetscInt :: ic(expected_max_level)
      PetscReal, parameter :: tol = 1.e-6_dp

      viewer = PETSC_NULL_VIEWER
      call thermo%init()
      json => fson_parse_mpi(str = json_str)
      call eos%init(json, thermo)
      call mesh%init(json)
      call mesh%configure(eos, gravity, json, viewer = viewer, err = err)
      call assert_equals(0, err, name // ": minc config error")
      call fson_destroy_mpi(json)

      orig_json_str = get_orig_json_str(json_str)
      orig_json => fson_parse_mpi(str = orig_json_str)
      call orig_mesh%init(orig_json)
      call orig_mesh%configure(eos, gravity, orig_json, viewer = viewer, err = err)
      call fson_destroy_mpi(orig_json)

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
         call DMGetStratumSize(mesh%dm, minc_level_label_name, m, &
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
      call DMPlexComputeGeometryFVM(orig_mesh%dm, orig_cell_geom, &
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
      call DMGetLabel(orig_mesh%dm, "ghost", ghost_label, ierr)
      CHKERRQ(ierr)
      call DMGetLabel(mesh%dm, minc_level_label_name, minc_level_label, ierr)
      CHKERRQ(ierr)
      call DMGetLabel(mesh%dm, cell_order_label_name, &
         cell_order_label, ierr)
    CHKERRQ(ierr)
      ic = 0
      h = 0
      do iminc = 1, size(mesh%minc)
         associate(minc => mesh%minc(iminc))
           call DMGetStratumSize(orig_mesh%dm, minc_zone_label_name, iminc, &
                num_minc_zone_cells, ierr); CHKERRQ(ierr)
           if (num_minc_zone_cells > 0) then
              call DMGetStratumIS(orig_mesh%dm, minc_zone_label_name, &
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
                    cell_p = mesh%strata(h)%minc_point(c, 0)
                    call section_offset(cell_section, cell_p, cell_offset, ierr)
                    CHKERRQ(ierr)
                    call cell%assign_geometry(cell_geom_array, cell_offset)
                    expected_vol = orig_cell%volume * minc%volume(1)
                    write(str, '(a, a, i3, a)') name, ": fracture volume(", cell_order, ")"
                    call assert_equals(expected_vol, cell%volume, tol, str)
                    do m = 1, minc%num_levels
                       cell_p = mesh%strata(h)%minc_point(ic(m), m)
                       call section_offset(cell_section, cell_p, &
                            cell_offset, ierr); CHKERRQ(ierr)
                       call cell%assign_geometry(cell_geom_array, cell_offset)
                       expected_vol = orig_cell%volume * minc%volume(m + 1)
                       write(str, '(a, a, i3, a, i1, a)') name, ": minc volume(", &
                            cell_order, ", ", m, ")"
                       call assert_equals(expected_vol, cell%volume, tol, str)
                       face_p = mesh%strata(h + 1)%minc_point(ic(m), m)
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
      call orig_mesh%destroy()
      call eos%destroy()
      call thermo%destroy()

    end subroutine minc_test

!........................................................................

    function get_orig_json_str(json_str) result(orig_str)
      ! Gets JSON string with MINC specification removed.

      character(*), intent(in) :: json_str
      character(len = len(json_str)) :: orig_str
      ! Locals:
      PetscInt :: i

      orig_str = json_str
      i = index(json_str, '"minc"')
      if (i > 0) then
         associate(n => len(json_str))
           orig_str(i:n) = '"minc": {}}}'
         end associate
      end if

    end function get_orig_json_str

!........................................................................

    PetscInt function total_interior_cell_count(mesh) result(n)
      ! Count interior cells.

      type(mesh_type), intent(in) :: mesh
      ! Locals:
      PetscInt :: c, n_local, start_cell, end_cell
      PetscErrorCode :: ierr

      call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
      CHKERRQ(ierr)

      n_local = 0
      do c = start_cell, end_cell - 1
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
    use IAPWS_module
    use eos_we_module

    type(fson_value), pointer :: json
    character(:), allocatable :: json_str
    PetscErrorCode :: ierr
    PetscMPIInt :: rank
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    PetscViewer :: viewer
    PetscInt, parameter :: num_rocktypes = 2

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call thermo%init()
    call eos%init(json, thermo)
    viewer = PETSC_NULL_VIEWER

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

    call eos%destroy()
    call thermo%destroy()

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
      call rock_dict%init(owner = PETSC_TRUE)
      call mesh%configure(eos, gravity, json, viewer = viewer, err = err)
      call setup_rock_vector(json, mesh%dm, mesh%cell_order, rock_vector, &
           rock_dict, rock_range_start, mesh%ghost_cell, err = err)
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
    use IAPWS_module
    use eos_we_module

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
      type(IAPWS_type) :: thermo
      type(eos_we_type) :: eos
      PetscViewer :: viewer
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

      viewer = PETSC_NULL_VIEWER
      expected_porosity = expected_matrix_porosity
      expected_porosity(0) = expected_fracture_porosity

      call thermo%init()
      json => fson_parse_mpi(str = json_str)
      call eos%init(json, thermo)
      call mesh%init(json)

      call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
      call mesh%configure(eos, gravity, json, viewer = viewer, err = err)
      call assert_equals(0, err, title // " mesh configure error")
      call rock_dict%init(owner = PETSC_TRUE)
      call setup_rock_vector(json, mesh%dm, mesh%cell_order, rock_vector, rock_dict, &
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
         call DMGetStratumSize(mesh%dm, minc_level_label_name, m, &
              num_cells, ierr); CHKERRQ(ierr)
         if (num_cells > 0) then
            call DMGetStratumIS(mesh%dm, minc_level_label_name, &
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

  subroutine test_cell_order
    ! cell_order AO

    use fson_mpi_module
    use IAPWS_module
    use eos_we_module
    use dm_utils_module, only: local_to_natural_cell_index, &
         global_vec_section, global_section_offset, &
         global_vec_range_start, local_to_natural_cell_index

    PetscMPIInt :: rank
    character(:), allocatable :: json_str
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json_str = '{"mesh": {' // &
         '"filename": "data/mesh/7x7grid.exo"}}'
    call cell_order_test_case(json_str, ' no bdy')

    json_str = '{"mesh": {' // &
         '"filename": "data/mesh/7x7grid.exo"}, ' // &
         '"boundaries": [{"faces": {"cells": [0, 1, 2, 3, 4, 5], ' // &
         '  "normal": [0, -1, 0]}}]' // &
         '}'
    call cell_order_test_case(json_str, 'bdy')

  contains

    subroutine cell_order_test_case(json_str, title)
      
      character(*), intent(in) :: json_str
      character(*), intent(in) :: title
      ! Locals:
      type(mesh_type) :: mesh
      type(IAPWS_type) :: thermo
      type(eos_we_type) :: eos
      type(fson_value), pointer :: json
      PetscViewer :: viewer
      PetscErrorCode :: err
      PetscInt :: c, start_cell, end_cell, offset, range_start
      PetscInt :: n, count, bs
      PetscInt, allocatable :: label_order(:), order(:)
      DMLabel :: label
      ISLocalToGlobalMapping :: l2g
      Vec :: v, v0
      PetscSection :: section
      PetscReal, pointer :: v_array(:)
      VecScatter :: scatter
      IS :: index0
      PetscInt, pointer :: ind(:)
      PetscInt, allocatable :: val(:)
      PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
      character(20), parameter :: label_name = "cell order"
      
      json => fson_parse_mpi(str = json_str)
      call thermo%init()
      call eos%init(json, thermo)
      viewer = PETSC_NULL_VIEWER
      call mesh%init(json)

      call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
      CHKERRQ(ierr)
      ! Create order label on serial DM:
      call DMCreateLabel(mesh%dm, label_name, ierr); CHKERRQ(ierr)
      do c = start_cell, end_cell - 1
         call DMSetLabelValue(mesh%dm, label_name, c, c, ierr); CHKERRQ(ierr)
      end do

      call mesh%configure(eos, gravity, json, viewer = viewer, err = err)
      call fson_destroy_mpi(json)

      ! Test distributed label values against mesh cell%order AO:
      call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
      CHKERRQ(ierr)
      allocate(label_order(start_cell: end_cell - 1), &
           order(start_cell: end_cell - 1))
      call DMGetLabel(mesh%dm, label_name, label, ierr); CHKERRQ(ierr)
      call DMGetLocalToGlobalMapping(mesh%dm, l2g, ierr); CHKERRQ(ierr)
      do c = start_cell, end_cell - 1
         call DMLabelGetValue(label, c, label_order(c), ierr); CHKERRQ(ierr)
      end do
      order = local_to_natural_cell_index(mesh%cell_order, l2g, &
           [(c, c = start_cell, end_cell - 1)])
      call assert_equals(label_order, order, end_cell - start_cell, &
           "cell order " // trim(title))
      deallocate(label_order, order)

      ! Test cell index IS:
      call DMGetGlobalVector(mesh%dm, v, ierr); CHKERRQ(ierr)
      call VecSet(v, 0._dp, ierr); CHKERRQ(ierr)
      call global_vec_section(v, section)
      call VecGetArrayF90(v, v_array, ierr); CHKERRQ(ierr)
      call global_vec_range_start(v, range_start)
      do c = start_cell, end_cell - 1
         if (mesh%ghost_cell(c) < 0) then
             call global_section_offset(section, c, &
                  range_start, offset, ierr); CHKERRQ(ierr)
             v_array(offset) = dble(local_to_natural_cell_index(&
                  mesh%cell_order, l2g, c))
         end if
      end do
      call VecRestoreArrayF90(v, v_array, ierr); CHKERRQ(ierr)
      call VecScatterCreateToZero(v, scatter, v0, ierr); CHKERRQ(ierr)
      call VecScatterBegin(scatter, v, v0, INSERT_VALUES, &
           SCATTER_FORWARD, ierr); CHKERRQ(ierr)
      call VecScatterEnd(scatter, v, v0, INSERT_VALUES, &
           SCATTER_FORWARD, ierr); CHKERRQ(ierr)
      call VecScatterDestroy(scatter, ierr); CHKERRQ(ierr)
      call ISAllGather(mesh%cell_index, index0, ierr); CHKERRQ(ierr)
      call ISGetIndicesF90(index0, ind, ierr); CHKERRQ(ierr)
      call VecGetArrayF90(v0, v_array, ierr); CHKERRQ(ierr)
      call VecGetBlockSize(v, bs, ierr); CHKERRQ(ierr)
      if (rank == 0) then
         count = size(ind)
         allocate(val(0: count - 1))
         do n = 0, count - 1
            val(n) = int(v_array(ind(n + 1) * bs + 1))
         end do
         call assert_equals([(n, n = 0, count - 1)], val, count, &
              "cell index " // trim(title))
         deallocate(val)
      end if
      call ISRestoreIndicesF90(index0, ind, ierr); CHKERRQ(ierr)
      call VecRestoreArrayF90(v0, v_array, ierr); CHKERRQ(ierr)
      call VecDestroy(v0, ierr); CHKERRQ(ierr)
      call DMRestoreGlobalVector(mesh%dm, v, ierr); CHKERRQ(ierr)

      call mesh%destroy()
      call eos%destroy()
      call thermo%destroy()

    end subroutine cell_order_test_case

  end subroutine test_cell_order

!------------------------------------------------------------------------

end module mesh_test
