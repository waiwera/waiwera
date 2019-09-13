module mesh_test

  ! Tests for mesh module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use fson
  use mesh_module

  implicit none
  private

  character(len = 512) :: data_path

  public :: setup, teardown
  public :: test_mesh_init, test_2d_cartesian_geometry, &
       test_2d_radial_geometry, test_mesh_face_permeability_direction, &
       test_setup_minc_dm, test_minc_rock, &
       test_rock_assignment, test_cell_natural_global, test_minc_cell_natural_global, &
       test_global_to_fracture_natural, test_redistribute

contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: ios

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()

    call get_environment_variable('WAIWERA_TEST_DATA_PATH', &
         data_path, status = ios)
    if (ios /= 0) data_path = ''

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown()

    PetscErrorCode :: ierr

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

  subroutine mesh_geometry_sanity_check(mesh, test, title)
    !! Returns true if mesh cell and face geometry vectors pass a
    !! basic sanity check.

    use cell_module, only: cell_type
    use face_module, only: face_type
    use dm_utils_module, only: local_vec_section, section_offset, &
         dm_get_num_partition_ghost_points, dm_get_end_interior_cell


    class(mesh_type), intent(in) :: mesh
    class(unit_test_type), intent(in out) :: test
    character(*), intent(in) :: title
    ! Locals:
    PetscInt :: start_cell, end_cell, c, offset
    PetscInt :: num_ghost_cells, end_non_ghost_cell, end_interior_cell
    PetscInt :: start_face, end_face, f, num_cells
    PetscSection :: cell_geom_section, face_geom_section
    PetscReal, pointer, contiguous :: cell_geom_array(:), face_geom_array(:)
    type(cell_type) :: cell
    type(face_type) :: face
    PetscInt :: dirn
    character(80) :: msg
    PetscErrorCode :: ierr
    PetscReal, parameter :: tol = 1.e-6_dp

    call local_vec_section(mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    if (end_cell > start_cell) then
       end_interior_cell = dm_get_end_interior_cell(mesh%dm, end_cell)
       num_ghost_cells = dm_get_num_partition_ghost_points(mesh%dm, 0)
       end_non_ghost_cell = end_interior_cell - num_ghost_cells
       call test%assert(all(mesh%ghost_cell(start_cell: end_non_ghost_cell - 1) < 0), &
            trim(title) // ": ghost cell array for non-ghost cells")
       call test%assert(all(mesh%ghost_cell(end_non_ghost_cell: end_interior_cell - 1) > 0), &
            trim(title) // ": ghost cell array for ghost cells")
       call cell%init(1, 1)
       offset = section_offset(cell_geom_section, end_cell - 1)
       call test%assert(offset + cell%dof - 1, size(cell_geom_array), &
            trim(title) // ": cell geom array size")
       do c = start_cell, end_interior_cell - 1
          if (mesh%ghost_cell(c) < 0) then
             offset = section_offset(cell_geom_section, c)
             call cell%assign_geometry(cell_geom_array, offset)
             write(msg, '(a, i4, a, e10.4)') " : cell ", c, " volume = ", cell%volume
             call test%assert(cell%volume > tol, trim(title) // trim(msg))
          end if
       end do
       call cell%destroy()
    end if

    call VecRestoreArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    call local_vec_section(mesh%face_geom, face_geom_section)
    call VecGetArrayReadF90(mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)

    call DMPlexGetHeightStratum(mesh%dm, 1, start_face, end_face, ierr)
    CHKERRQ(ierr)
    call face%init(1, 1)
    do f = start_face, end_face - 1
       if (mesh%ghost_face(f) < 0) then
          call DMPlexGetSupportSize(mesh%dm, f, num_cells, ierr); CHKERRQ(ierr)
          write(msg, '(a, i0, a)') " : face ", f, " num cells"
          call test%assert(2, num_cells, trim(title) // trim(msg))
          offset = section_offset(face_geom_section, f)
          call face%assign_geometry(face_geom_array, offset)
          write(msg, '(a, i0, a, e10.4)') " : face ", f, " area = ", face%area
          call test%assert(face%area > tol, trim(title) // trim(msg))
          write(msg, '(a, i0, a, e10.4)') " : face ", f, " distance12 = ", face%distance12
          call test%assert(abs(face%distance12) > tol, trim(title) // trim(msg))
          write(msg, '(a, i0, a)') " : face ", f, " distance12 = sum(distance)"
          call test%assert(face%distance12, sum(face%distance), trim(title) // trim(msg))
          dirn = nint(face%permeability_direction)
          write(msg, '(a, i0, a, i2)') " : face ", f, " perm dirn =", dirn
          call test%assert(((1 <= dirn) .and. (dirn <= 3)), trim(title) // trim(msg))
       end if
    end do

    call face%destroy()
    call VecRestoreArrayReadF90(mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)

  end subroutine mesh_geometry_sanity_check

!------------------------------------------------------------------------

  subroutine test_mesh_init(test)

    ! Mesh init test

    use dm_utils_module, only: section_offset, local_to_natural_cell_index
    use fson_mpi_module
    use cell_module
    use face_module
    use IAPWS_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(mesh_type) :: mesh
    Vec :: x
    type(face_type) :: face
    PetscInt :: global_solution_dof

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

    json => fson_parse_mpi(str = '{"mesh": "' // trim(adjustl(data_path)) // 'mesh/block3.exo"}')
    call eos%init(json, thermo)
    call mesh%init(eos, json)
    call mesh%configure(gravity, json, err = err)
    call fson_destroy_mpi(json)
    call mesh%construct_ghost_cells(gravity)
    call mesh%destroy_distribution_data()

    if (rank == 0) then
       call test%assert(expected_dim, mesh%dim, "mesh dimension")
    end if
    call mesh_geometry_sanity_check(mesh, test, "")

    call DMGetGlobalVector(mesh%dm, x, ierr); CHKERRQ(ierr)
    call VecGetSize(x, global_solution_dof, ierr); CHKERRQ(ierr)
    if (rank == 0) then
       call test%assert(num_cells * eos%num_primary_variables, &
            global_solution_dof, "global solution dof")
    end if
    call DMRestoreGlobalVector(mesh%dm, x, ierr); CHKERRQ(ierr)

    call face%init()
    call DMGetLocalToGlobalMapping(mesh%dm, l2g, ierr); CHKERRQ(ierr)

    ! Check face geometry:
    call VecGetDM(mesh%face_geom, dm_face, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(mesh%face_geom, fg, ierr); CHKERRQ(ierr)
    call DMGetSection(dm_face, section, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(mesh%dm, 1, fstart, fend, ierr); CHKERRQ(ierr)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    CHKERRQ(ierr)
    do f = fstart, fend - 1
       call DMLabelGetValue(ghost_label, f, ghost_face, ierr); CHKERRQ(ierr)
       if (ghost_face < 0) then
          offset = section_offset(section, f)
          call face%assign_geometry(fg, offset)
          write(msg, '(a, i2)') 'face area ', f
          call test%assert(face_area, face%area, msg)
          dist = face%distance
          call DMPlexGetSupport(mesh%dm, f, cells, ierr); CHKERRQ(ierr)
          order = local_to_natural_cell_index(mesh%cell_natural_global, l2g, cells)
          if (all(order == [0,1])) then
             gf = 1
          else if (all(order == [1,2])) then
             gf = 2
          else
             gf = 0
          end if
          if (gf > 0) then
             write(msg, '(a, i2)') 'face distance ', f
             call test%assert(face_distance(:, gf), dist, msg)
             write(msg, '(a, i2)') 'face centroid ', f
             call test%assert(face_centroid(:, gf), face%centroid, msg)
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

  subroutine test_2d_cartesian_geometry(test)
    ! Test 2D Cartesian volumes and areas

    use fson_mpi_module
    use cell_module
    use face_module
    use dm_utils_module, only: section_offset, local_vec_section
    use IAPWS_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(mesh_type) :: mesh
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
    PetscReal, parameter :: tol = 1.e-6_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call thermo%init()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/2D.msh",' // &
         '"thickness": 100.}}')
    call eos%init(json, thermo)
    call mesh%init(eos, json)
    call mesh%configure(gravity, json, err = err)
    call mesh%construct_ghost_cells(gravity)

    call fson_destroy_mpi(json)
    call mesh%destroy_distribution_data()

    call mesh_geometry_sanity_check(mesh, test, "")

    call local_vec_section(mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    volumes_OK = PETSC_TRUE

    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)

    do c = start_cell, end_cell - 1
       if (mesh%ghost_cell(c) < 0) then
          offset = section_offset(cell_geom_section, c)
          call cell%assign_geometry(cell_geom_array, offset)
          volumes_OK = (volumes_OK .and. &
               abs(cell%volume - expected_cell_vol) <= tol)
       end if
    end do

    write(msg, '(a, i3)') 'cell volumes on proc', rank
    call test%assert(volumes_OK, msg)

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
          offset = section_offset(face_geom_section, f)
          call face%assign_geometry(face_geom_array, offset)
          areas_OK = (areas_OK .and. &
               abs(face%area - expected_face_area) <= tol)
       end if
    end do

    write(msg, '(a, i3)') 'face areas on proc', rank
    call test%assert(areas_OK, msg)

    call VecRestoreArrayReadF90(mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)
    call face%destroy()

    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()

  end subroutine test_2d_cartesian_geometry

!------------------------------------------------------------------------

  subroutine test_2d_radial_geometry(test)
    ! Test 2D radial volumes and areas

    use fson_mpi_module
    use cell_module
    use face_module
    use dm_utils_module, only: section_offset, local_vec_section
    use utils_module, only: pi
    use IAPWS_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
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
    PetscReal, parameter :: tol = 1.e-6_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call thermo%init()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/2D.msh",' // &
         '"radial": true}}')
    call eos%init(json, thermo)
    call mesh%init(eos, json)
    call mesh%configure(gravity, json, err = err)
    call mesh%construct_ghost_cells(gravity)

    call fson_destroy_mpi(json)
    call mesh%destroy_distribution_data()

    call mesh_geometry_sanity_check(mesh, test, "")

    call local_vec_section(mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    volumes_OK = PETSC_TRUE
    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)

    do c = start_cell, end_cell - 1
       if (mesh%ghost_cell(c) < 0) then
          offset = section_offset(cell_geom_section, c)
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
    call test%assert(volumes_OK, msg)

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
          offset = section_offset(face_geom_section, f)
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
    call test%assert(areas_OK, msg)

    call VecRestoreArrayReadF90(mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)
    call face%destroy()

    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()

  end subroutine test_2d_radial_geometry

!------------------------------------------------------------------------

  subroutine test_mesh_face_permeability_direction(test)
    ! Test face permeability direction

    use fson_mpi_module
    use dm_utils_module, only: local_vec_section, section_offset
    use face_module
    use IAPWS_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(mesh_type) :: mesh
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

    json => fson_parse_mpi(str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo", ' // &
         '"faces": [' // &
         '{"cells": [16, 23], "permeability direction": 1}' // &
         ']}}')
    call eos%init(json, thermo)
    call mesh%init(eos, json)

    call mesh%configure(gravity, json, err = err)
    call mesh%construct_ghost_cells(gravity)

    call mesh%override_face_properties()
    call fson_destroy_mpi(json)
    call mesh%destroy_distribution_data()

    call mesh_geometry_sanity_check(mesh, test, "")

    call local_vec_section(mesh%face_geom, face_geom_section)
    call VecGetArrayReadF90(mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)
    call face%init()
    call DMPlexGetHeightStratum(mesh%dm, 1, start_face, end_face, ierr)
    CHKERRQ(ierr)

    do f = start_face, end_face - 1
       if (mesh%ghost_face(f) < 0) then
          offset = section_offset(face_geom_section, f)
          call face%assign_geometry(face_geom_array, offset)
          dist = norm2(face%centroid - position)
          if (dist <= tol) then
             call test%assert(expected_direction, &
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

  PetscInt function total_interior_cell_count(mesh) result(n)
    ! Count mesh interior cells.

    use dm_utils_module, only: dm_get_end_interior_cell

    type(mesh_type), intent(in) :: mesh
    ! Locals:
    PetscInt :: c, n_local, start_cell, end_cell, end_interior_cell
    PetscErrorCode :: ierr

    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(mesh%dm, end_cell)

    n_local = 0
    do c = start_cell, end_interior_cell - 1
       if (mesh%ghost_cell(c) < 0) then
          n_local = n_local + 1
       end if
    end do
    call MPI_reduce(n_local, n, 1, MPI_INTEGER, MPI_SUM, &
         0, PETSC_COMM_WORLD, ierr)

  end function total_interior_cell_count

!------------------------------------------------------------------------

  PetscInt function total_interior_cell_count_sf(mesh) result(n)
    ! Count mesh interior cells using point SF.

    use dm_utils_module, only: dm_get_end_interior_cell

    type(mesh_type), intent(in) :: mesh
    ! Locals:
    PetscMPIInt :: np
    PetscInt :: c, n_local, start_cell, end_cell, end_interior_cell
    PetscSF :: point_sf
    PetscInt :: num_roots, num_leaves
    PetscInt, pointer :: local(:)
    type(PetscSFNode), pointer :: remote(:)
    PetscErrorCode :: ierr

    call MPI_comm_size(PETSC_COMM_WORLD, np, ierr)
    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(mesh%dm, end_cell)

    if (np == 1) then
       n_local = end_interior_cell - start_cell
    else
       n_local = 0
       call DMGetPointSF(mesh%dm, point_sf, ierr); CHKERRQ(ierr)
       call PetscSFGetGraph(point_sf, num_roots, num_leaves, &
            local, remote, ierr); CHKERRQ(ierr)
       do c = start_cell, end_interior_cell - 1
          if (local(c + 1) >= end_cell) then
             n_local = n_local + 1
          end if
       end do
    end if
    call MPI_reduce(n_local, n, 1, MPI_INTEGER, MPI_SUM, &
         0, PETSC_COMM_WORLD, ierr)

  end function total_interior_cell_count_sf

!------------------------------------------------------------------------

  PetscInt function total_interior_cell_label_count(mesh, label_name, &
       label_value) result(n)
    ! Count mesh interior cells with specified label value.

    use dm_utils_module, only: dm_get_end_interior_cell

    type(mesh_type), intent(in) :: mesh
    character(*), intent(in) :: label_name
    PetscInt, intent(in) :: label_value
    ! Locals:
    PetscInt :: c, n_local, start_cell, end_cell, end_interior_cell, val
    DMLabel :: label
    PetscErrorCode :: ierr

    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(mesh%dm, end_cell)
    call DMGetLabel(mesh%dm, label_name, label, ierr); CHKERRQ(ierr)

    n_local = 0
    do c = start_cell, end_interior_cell - 1
       if (mesh%ghost_cell(c) < 0) then
          call DMLabelGetValue(label, c, val, ierr); CHKERRQ(ierr)
          if (val == label_value) then
             n_local = n_local + 1
          end if
       end if
    end do
    call MPI_reduce(n_local, n, 1, MPI_INTEGER, MPI_SUM, &
         0, PETSC_COMM_WORLD, ierr)

  end function total_interior_cell_label_count

!------------------------------------------------------------------------

  subroutine dm_dag_sanity_check(dm, test, title)
    !! Sanity check for DM graph.

    DM, intent(in) :: dm
    class(unit_test_type), intent(in out) :: test
    character(*), intent(in) :: title
    ! Locals:
    PetscInt :: p, h, depth!, i
    PetscInt, allocatable :: pstart(:), pend(:)
    PetscInt, pointer :: cone(:), support(:)
    PetscMPIInt :: rank
    character(80) :: msg
    PetscErrorCode :: ierr

    call mpi_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call DMPlexGetDepth(dm, depth, ierr); CHKERRQ(ierr)
    allocate(pstart(0:depth), pend(0:depth))
    do h = 0, depth
       call DMPlexGetHeightStratum(dm, h, pstart(h), pend(h), &
            ierr); CHKERRQ(ierr)
    end do

    do h = 0, depth
       do p = pstart(h), pend(h) - 1
          call DMPlexGetCone(dm, p, cone, ierr)
          call DMPlexGetSupport(dm, p, support, ierr)
          if (h < depth) then
             write(msg, '(a, i0, a, i0)') ': cone for h:', h, 'p:', p
             call test%assert(all((pstart(h + 1) <= cone) .and. &
                  (cone < pend(h + 1))), trim(title) // msg)
          end if
          if (h > 0) then
             write(msg, '(a, i0, a, i0)') ': support for h:', h, 'p:', p
             call test%assert(all((pstart(h - 1) <= support) .and. &
                  (support < pend(h - 1))), trim(title) // msg)
          end if
          call DMPlexRestoreCone(dm, p, cone, ierr)
          call DMPlexRestoreSupport(dm, p, support, ierr)
       end do
    end do

    deallocate(pstart, pend)

  end subroutine dm_dag_sanity_check

!------------------------------------------------------------------------

  subroutine minc_mesh_natural_check(mesh, test, title)
    !! Checks natural cell index data structures on MINC mesh.

    use dm_utils_module, only: local_to_natural_cell_index
    use minc_module, only: minc_level_label_name

    class(mesh_type), intent(in) :: mesh
    class(unit_test_type), intent(in out) :: test
    character(*), intent(in) :: title
    ! Locals:
    ISLocalToGlobalMapping :: l2g
    DMLabel :: ghost_label, minc_label
    PetscInt :: start_cell, end_cell, c
    PetscInt :: ghost, minc_level, natural, parent_natural
    character(80) :: msg
    PetscErrorCode :: ierr

    call DMGetLocalToGlobalMapping(mesh%dm, l2g, ierr); CHKERRQ(ierr)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call DMGetLabel(mesh%dm, minc_level_label_name, minc_label, ierr)
    CHKERRQ(ierr)

    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, &
         ierr); CHKERRQ(ierr)

    do c = start_cell, end_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       call DMLabelGetValue(minc_label, c, minc_level, ierr); CHKERRQ(ierr)
       if ((ghost < 0) .and. (minc_level == 0)) then
          natural = local_to_natural_cell_index(mesh%cell_natural_global, l2g, c)
          parent_natural = mesh%local_to_parent_natural(c)
          write(msg, '(a, i0)') ': natural ', c
          call test%assert(parent_natural, natural, trim(title) // msg)
       end if
    end do

  end subroutine minc_mesh_natural_check

!------------------------------------------------------------------------

  subroutine test_setup_minc_dm(test)
    ! Test setup_minc_dm

    use fson_mpi_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    character(:), allocatable :: json_str
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"rock": {"zones": ["all"]}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}}}}}'
    call minc_test('all', json_str, 1, 2 * 49, 1, [49, 49])

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"left": {"x": [0, 1500]}},' // &
         '  "minc": {"rock": {"zones": ["left"]}, ' // &
         '           "geometry": {"matrix": {"volume": 0.9}}}}}'
    call minc_test('partial', json_str, 1, 63, 1, [49, 14])

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"left": {"x": [0, 1500]}, "right": {"-": "left"}},' // &
         '  "minc": [{"rock": {"zones": ["left"]}, ' // &
         '                     "geometry": {"fracture": {"volume": 0.1}}}, ' // &
         '           {"rock": {"zones": ["right"]}, ' // &
         '            "geometry": {"fracture": {"volume": 0.1}, ' // &
         '                         "matrix": {"volume": [0.3, 0.6]}}}]}}'
    call minc_test('two-zone', json_str, 2, 133, 2, [49, 49, 35])

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"left": {"x": [0, 1500]}, ' // &
         '            "right corner": {"x": [2500, 4500], "y": [3000, 4500]}},' // &
         '  "minc": [{"rock": {"zones": ["right corner"]}, ' // &
         '            "geometry": {"fracture": {"volume": 0.1}}}, ' // &
         '           {"rock": {"zones": ["left"]}, ' // &
         '            "geometry": {"matrix": {"volume": [0.3, 0.6]}}}]}}'
    call minc_test('two-zone partial', json_str, 2, 83, 2, [49, 20, 14])

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"left": {"x": [0, 1500]}, ' // &
         '            "right": {"-": "left"}},' // &
         '  "minc": {"rock": [{"zones": ["left"]}, {"zones": ["right"]}], ' // &
         '            "geometry": {"fracture": {"volume": 0.1}}}, ' // &
         '           }}'
    call minc_test('two sub-zone', json_str, 1, 2 * 49, 1, [49, 49])

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"left": {"x": [0, 1500]}, "right": {"-": "left"}},' // &
         '  "minc": [{"rock": {"types": ["rock1"]}, ' // &
         '                     "geometry": {"fracture": {"volume": 0.1}}}, ' // &
         '           {"rock": {"types": ["rock2"]}, ' // &
         '            "geometry": {"fracture": {"volume": 0.1}, ' // &
         '                         "matrix": {"volume": [0.3, 0.6]}}}]}, ' // &
         ' "rock": {"types": [{"name": "rock1", "zones": "left"}, ' // &
         '                    {"name": "rock2", "zones": ["right"]}]}' // &
         '}'
    call minc_test('rocktype', json_str, 2, 133, 2, [49, 49, 35])

  contains

!........................................................................

    subroutine minc_test(name, json_str, expected_num_zones, expected_num_cells, &
         expected_max_level, expected_num_minc_level_cells)

      use minc_module, only: minc_level_label_name, minc_zone_label_name
      use cell_module
      use face_module
      use dm_utils_module, only: section_offset, local_vec_section, &
           local_to_natural_cell_index
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
      PetscInt :: num_cells, num_cells_sf, num_minc_zones, m, num, max_num_levels
      PetscErrorCode :: err
      PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
      character(48) :: str
      PetscInt :: iminc, num_minc_zone_cells, i, c
      IS :: minc_IS
      PetscInt, pointer :: minc_cells(:)
      Vec :: orig_cell_geom, orig_face_geom
      type(cell_type) :: orig_cell, cell
      type(face_type) :: face
      ISLocalToGlobalMapping :: l2g
      DMLabel :: ghost_label
      PetscInt, parameter :: nc = 1, np = 1 ! dummy values for cell init
      PetscSection :: orig_cell_section, cell_section, face_section
      PetscReal, pointer, contiguous :: orig_cell_geom_array(:), cell_geom_array(:)
      PetscReal, pointer, contiguous :: face_geom_array(:)
      PetscInt :: orig_cell_offset, cell_offset, ghost, order
      PetscInt :: face_offset, cell_p, face_p, h
      PetscReal :: expected_vol, expected_area
      PetscInt :: ic(expected_max_level)

      call thermo%init()
      json => fson_parse_mpi(str = json_str)
      call eos%init(json, thermo)
      call mesh%init(eos, json)
      call mesh%configure(gravity, json, err = err)
      call mesh%construct_ghost_cells(gravity)

      call test%assert(0, err, name // ": minc config error")
      call fson_destroy_mpi(json)
      call mesh%destroy_distribution_data()

      call mesh_geometry_sanity_check(mesh, test, name)
      call dm_dag_sanity_check(mesh%dm, test, name)

      orig_json_str = get_orig_json_str(json_str)
      orig_json => fson_parse_mpi(str = orig_json_str)
      call orig_mesh%init(eos, orig_json)
      call orig_mesh%configure(gravity, orig_json, err = err)
      call orig_mesh%construct_ghost_cells(gravity)

      call fson_destroy_mpi(orig_json)
      call orig_mesh%destroy_distribution_data()

      call DMGetLocalToGlobalMapping(orig_mesh%dm, l2g, ierr); CHKERRQ(ierr)

      if (rank == 0) then
         call test%assert(mesh%has_minc, name // ": mesh has minc")
      end if
      num_minc_zones = size(mesh%minc)
      max_num_levels = maxval(mesh%minc%num_levels)
      if (rank == 0) then
         call test%assert(expected_num_zones, &
              num_minc_zones, name // ": num minc zones")
         call test%assert(expected_max_level, &
              max_num_levels, name // ": num minc levels")
      end if

      num_cells = total_interior_cell_count(mesh)
      num_cells_sf = total_interior_cell_count_sf(mesh)
      if (rank == 0) then
         call test%assert(expected_num_cells, &
              num_cells, name  // ": num cells")
         call test%assert(expected_num_cells, &
              num_cells_sf, name  // ": num cells SF")
      end if

      do m = 0, expected_max_level
         num = total_interior_cell_label_count(mesh, &
              minc_level_label_name, m)
         if (rank == 0) then
            write(str, '(a, i1)') ": num minc points, level ", m
            call test%assert(expected_num_minc_level_cells(m), &
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
                    order = local_to_natural_cell_index(orig_mesh%cell_natural_global, l2g, c)
                    orig_cell_offset = section_offset(orig_cell_section, c)
                    call orig_cell%assign_geometry(orig_cell_geom_array, orig_cell_offset)
                    cell_p = mesh%strata(h)%minc_point(c, 0)
                    cell_offset = section_offset(cell_section, cell_p)
                    call cell%assign_geometry(cell_geom_array, cell_offset)
                    expected_vol = orig_cell%volume * minc%volume(1)
                    write(str, '(a, a, i3, a)') name, ": fracture volume(", order, ")"
                    call test%assert(expected_vol, cell%volume, str)
                    do m = 1, minc%num_levels
                       cell_p = mesh%strata(h)%minc_point(ic(m), m)
                       cell_offset = section_offset(cell_section, cell_p)
                       call cell%assign_geometry(cell_geom_array, cell_offset)
                       expected_vol = orig_cell%volume * minc%volume(m + 1)
                       write(str, '(a, a, i3, a, i1, a)') name, ": minc volume(", &
                            order, ", ", m, ")"
                       call test%assert(expected_vol, cell%volume, str)
                       face_p = mesh%strata(h + 1)%minc_point(ic(m), m)
                       face_offset = section_offset(face_section, face_p)
                       call face%assign_geometry(face_geom_array, face_offset)
                       expected_area = orig_cell%volume * minc%connection_area(m)
                       write(str, '(a, a, i3, a, i1, a, i1, a)') name, &
                            ": minc area(", order, ", ", m-1, ":", m, ")"
                       call test%assert(expected_area, face%area, str)
                       write(str, '(a, a, i3, a, i1, a, i1, a)') name, ": minc distances(", &
                            order, ", ", m-1, ":", m, ")"
                       call test%assert(minc%connection_distance(m: m + 1), face%distance, &
                            str)
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

  end subroutine test_setup_minc_dm

!------------------------------------------------------------------------

  subroutine test_rock_assignment(test)
    ! rock assignment

    use fson_mpi_module
    use rock_module
    use dm_utils_module, only: global_section_offset, global_vec_section
    use dictionary_module
    use IAPWS_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    character(:), allocatable :: json_str
    PetscErrorCode :: ierr
    PetscMPIInt :: rank
    PetscInt, parameter :: num_rocktypes = 2

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo"}, ' // &
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
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo", ' // &
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
         '{"mesh": {"filename": "' //trim(adjustl(data_path)) // 'mesh/7x7grid.exo", ' // &
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
      type(IAPWS_type) :: thermo
      type(eos_we_type) :: eos
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

      call thermo%init()
      json => fson_parse_mpi(str = json_str)
      call eos%init(json, thermo)
      call mesh%init(eos, json)

      call rock_dict%init(owner = PETSC_TRUE)
      call mesh%configure(gravity, json, err = err)
      call mesh%construct_ghost_cells(gravity)

      call mesh_geometry_sanity_check(mesh, test, title)

      call setup_rock_vector(json, mesh%dm, rock_vector, &
           rock_dict, rock_range_start, mesh%ghost_cell, err = err)
      call test%assert(0, err, "setup rock vector error")
      call fson_destroy_mpi(json)
      call mesh%destroy_distribution_data()

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
            offset = global_section_offset(section, c, rock_range_start)
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
         call test%assert(expected_count, rock_count, &
              "rock counts: " // trim(title))
      end if

      call rock%destroy()
      call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)

      call VecDestroy(rock_vector, ierr); CHKERRQ(ierr)
      call mesh%destroy()
      call rock_dict%destroy()
      call eos%destroy()
      call thermo%destroy()

    end subroutine rock_test_case

  end subroutine test_rock_assignment

!------------------------------------------------------------------------

  subroutine test_minc_rock(test)
    ! MINC rock assignment

    use fson_mpi_module
    use dictionary_module
    use rock_module
    use minc_module, only: minc_level_label_name, minc_rocktype_zone_label_name
    use dm_utils_module, only: global_section_offset, global_vec_section, &
         natural_to_local_cell_index
    use IAPWS_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscErrorCode :: ierr
    character(:), allocatable :: json_str

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"rock": {"zones": "all", ' // &
         '                    "fracture": {"type": "fracture"}, ' // &
         '                    "matrix": {"type": "matrix"}}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1, "planes": 3, ' // &
         '                          "spacing": 100}, ' // &
         '                        "matrix": {"volume": 0.9}}}},' // &
         ' "rock": {"types": [{"name": "original", "porosity": 0.1, "zones": "all"}, ' // &
         '                    {"name": "fracture", "porosity": 0.6}, ' // &
         '                    {"name": "matrix", "porosity": 0.02}]}}'

    call minc_rock_test_case(json_str, "case 1", 1, [0.6_dp], [0.02_dp], [49])

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"rock": {"zones": "all", ' // &
         '                    "fracture": {"type": "fracture"}, ' // &
         '                    "matrix": {"type": "matrix"}}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1, "planes": 3, ' // &
         '                          "spacing": 100}, ' // &
         '                        "matrix": {"volume": [0.3, 0.6]}}}},' // &
         ' "rock": {"types": [{"name": "original", "porosity": 0.1, "zones": "all"}, ' // &
         '                    {"name": "fracture", "porosity": 0.6}, ' // &
         '                    {"name": "matrix"}]}}'

    call minc_rock_test_case(json_str, "case 2", 2, [0.6_dp], [2._dp / 45._dp], [49])

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}, "S": {"y": [0, 1500]}, "N": {"-": "S"}},' // &
         '  "minc": {"rock": [{"zones": "S", ' // &
         '                      "fracture": {"type": "fractureS"}, ' // &
         '                      "matrix": {"type": "matrixS"}}, ' // &
         '                    {"zones": "N", ' // &
         '                       "fracture": {"type": "fractureN"}, ' // &
         '                       "matrix": {"type": "matrixN"}}], ' // &
         '           "geometry": {"fracture": {"volume": 0.1, "planes": 3, ' // &
         '                          "spacing": 100}, ' // &
         '                        "matrix": {"volume": 0.9}}}},' // &
         ' "rock": {"types": [{"name": "original", "porosity": 0.1, "zones": "all"}, ' // &
         '                    {"name": "fractureS", "porosity": 0.6}, ' // &
         '                    {"name": "matrixS", "porosity": 0.02}, ' // &
         '                    {"name": "fractureN", "porosity": 0.7}, ' // &
         '                    {"name": "matrixN"}]}}'

    call minc_rock_test_case(json_str, "case 3", 1, [0.6_dp, 0.7_dp], &
         [0.02_dp, 1._dp / 30._dp], [14, 35])

  contains

    subroutine minc_rock_test_case(json_str, title, num_levels, &
         expected_fracture_porosity, expected_matrix_porosity, &
         expected_num_minc_rock_cells)

      character(*), intent(in) :: json_str
      character(*), intent(in) :: title
      PetscInt, intent(in) :: num_levels
      PetscReal, intent(in) :: expected_fracture_porosity(:), expected_matrix_porosity(:)
      PetscInt, intent(in) :: expected_num_minc_rock_cells(:)
      ! Locals:
      type(fson_value), pointer :: json
      type(mesh_type) :: mesh
      type(IAPWS_type) :: thermo
      type(eos_we_type) :: eos
      PetscInt :: num_local_minc_rock_cells, num_minc_rock_cells, num_minc_rocktypes
      Vec :: rock_vector
      type(dictionary_type) :: rock_dict
      PetscInt :: fracture_natural, fracture_local, r
      ISLocalToGlobalMapping :: l2g
      DMLabel :: minc_rocktype_label
      PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
      PetscErrorCode :: err
      PetscInt :: rock_range_start, c, i, m, num_cells, offset
      PetscReal, contiguous, pointer :: rock_array(:)
      type(rock_type) :: rock
      PetscSection :: section
      IS :: minc_IS
      PetscInt, pointer :: minc_points(:)
      PetscReal :: expected_porosity
      character(8) :: levelstr

      num_minc_rocktypes = size(expected_fracture_porosity)

      call thermo%init()
      json => fson_parse_mpi(str = json_str)
      call eos%init(json, thermo)
      call mesh%init(eos, json)
      call mesh%configure(gravity, json, err = err)
      call mesh%construct_ghost_cells(gravity)

      call test%assert(0, err, title // " mesh configure error")

      call rock_dict%init(owner = PETSC_TRUE)
      call setup_rock_vector(json, mesh%dm, rock_vector, rock_dict, &
           rock_range_start, mesh%ghost_cell, err = err)
      call test%assert(0, err, title // " setup rock vector error")
      call mesh%setup_minc_rock_properties(json, rock_vector, &
           rock_range_start, err = err)
      call test%assert(0, err, title // " setup MINC rock properties error")
      call fson_destroy_mpi(json)
      call mesh%destroy_distribution_data()

      call mesh_geometry_sanity_check(mesh, test, title)

      call VecGetArrayReadF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
      call global_vec_section(rock_vector, section)
      call rock%init()
      call DMGetLocalToGlobalMapping(mesh%dm, l2g, ierr); CHKERRQ(ierr)
      call DMGetLabel(mesh%dm, minc_rocktype_zone_label_name, minc_rocktype_label, ierr)

      do r = 1, num_minc_rocktypes
         call DMGetStratumSize(mesh%dm, minc_rocktype_zone_label_name, r, &
              num_local_minc_rock_cells, ierr); CHKERRQ(ierr)
         call mpi_reduce(num_local_minc_rock_cells, num_minc_rock_cells, 1, &
              MPI_INTEGER, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
         if (rank == 0) then
            call test%assert(expected_num_minc_rock_cells(r), num_minc_rock_cells, &
              title // " MINC rock cell count")
         end if
      end do

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
                  fracture_natural = mesh%local_to_parent_natural(c)
                  fracture_local = natural_to_local_cell_index(mesh%cell_natural_global, &
                       l2g, fracture_natural)
                  call DMLabelGetValue(minc_rocktype_label, fracture_local, &
                       r, ierr); CHKERRQ(ierr)
                  if (m == 0) then
                     expected_porosity = expected_fracture_porosity(r)
                  else
                     expected_porosity = expected_matrix_porosity(r)
                  end if
                  offset = global_section_offset(section, c, rock_range_start)
                  call rock%assign(rock_array, offset)
                  call test%assert(expected_porosity, rock%porosity, &
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

  subroutine test_cell_natural_global(test)
    ! cell_natural_global AO

    use fson_mpi_module
    use IAPWS_module
    use eos_we_module
    use dm_utils_module, only: local_to_natural_cell_index, &
         global_vec_section, global_section_offset, &
         global_vec_range_start, local_to_natural_cell_index

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    character(:), allocatable :: json_str
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json_str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo"}}'
    call cell_natural_global_test_case(json_str, ' no bdy')

    json_str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo"}, ' // &
         '"boundaries": [{"faces": {"cells": [0, 1, 2, 3, 4, 5], ' // &
         '  "normal": [0, -1, 0]}}]' // &
         '}'
    call cell_natural_global_test_case(json_str, 'bdy')

  contains

    subroutine cell_natural_global_test_case(json_str, title)
      
      character(*), intent(in) :: json_str
      character(*), intent(in) :: title
      ! Locals:
      type(mesh_type) :: mesh
      type(IAPWS_type) :: thermo
      type(eos_we_type) :: eos
      type(fson_value), pointer :: json
      PetscErrorCode :: err
      PetscInt :: c, start_cell, end_cell
      PetscInt, allocatable :: label_order(:), order(:)
      DMLabel :: label
      ISLocalToGlobalMapping :: l2g
      PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
      character(20), parameter :: label_name = "cell order"
      
      json => fson_parse_mpi(str = json_str)
      call thermo%init()
      call eos%init(json, thermo)
      call mesh%init(eos, json)

      call DMPlexGetHeightStratum(mesh%serial_dm, 0, start_cell, end_cell, ierr)
      CHKERRQ(ierr)
      ! Create order label on serial DM:
      call DMCreateLabel(mesh%serial_dm, label_name, ierr); CHKERRQ(ierr)
      do c = start_cell, end_cell - 1
         call DMSetLabelValue(mesh%serial_dm, label_name, c, c, ierr); CHKERRQ(ierr)
      end do

      call mesh%configure(gravity, json, err = err)
      call test%assert(0, err, "mesh config " // trim(title))

      call mesh%construct_ghost_cells(gravity)
      call mesh%destroy_distribution_data()
      call fson_destroy_mpi(json)

      call mesh_geometry_sanity_check(mesh, test, title)

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
      order = local_to_natural_cell_index(mesh%cell_natural_global, l2g, &
           [(c, c = start_cell, end_cell - 1)])
      call test%assert(label_order, order, "cell order " // trim(title))
      deallocate(label_order, order)

      call mesh%destroy()
      call eos%destroy()
      call thermo%destroy()

    end subroutine cell_natural_global_test_case

  end subroutine test_cell_natural_global

!------------------------------------------------------------------------

  subroutine test_minc_cell_natural_global(test)
    ! MINC cell order

    use fson_mpi_module
    use IAPWS_module
    use eos_we_module
    use dm_utils_module, only: dm_get_end_interior_cell, local_to_natural_cell_index
    use minc_module, only: minc_zone_label_name

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    character(:), allocatable :: json_str
    PetscInt, allocatable :: expected_minc_order(:,:)
    PetscInt :: c, i
    PetscInt, allocatable :: natural(:), minc_natural(:)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"rock": {"zones": ["all"]}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}}}}}'
    allocate(expected_minc_order(1, 0:48))
    expected_minc_order(1, :) = [(c + 49, c = 0, 48)]
    call minc_cell_natural_global_test_case(json_str, ' full no bdy', expected_minc_order)
    deallocate(expected_minc_order)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"rock": {"zones": ["all"]}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}}}}, ' // &
         '"boundaries": [{"faces": {"cells": [0, 1, 2, 3, 4, 5], ' // &
         '  "normal": [0, -1, 0]}}]' // &
         '}'
    allocate(expected_minc_order(1, 0:48))
    expected_minc_order(1, :) = [(c + 49, c = 0, 48)]
    call minc_cell_natural_global_test_case(json_str, ' full bdy', expected_minc_order)
    deallocate(expected_minc_order)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"sw": {"x": [0, 3000], "y": [0, 1500]}},' // &
         '  "minc": {"rock": {"zones": ["sw"]}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}}}}}'
    allocate(expected_minc_order(1, 0:48))
    expected_minc_order = 0
    expected_minc_order(1, 0:4) = [(c + 49, c = 0, 4)]
    expected_minc_order(1, 7:11) = [(c + 47, c = 7, 11)]
    call minc_cell_natural_global_test_case(json_str, ' partial no bdy', expected_minc_order)
    deallocate(expected_minc_order)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"sw": {"x": [0, 3000], "y": [0, 1500]}},' // &
         '  "minc": {"rock": {"zones": ["sw"]}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}}}},' // &
         '"boundaries": [{"faces": {"cells": [0, 1, 2, 3, 4, 5], ' // &
         '  "normal": [0, -1, 0]}}]' // &
         '}'
    allocate(expected_minc_order(1, 0:48))
    expected_minc_order = 0
    expected_minc_order(1, 0:4) = [(c + 49, c = 0, 4)]
    expected_minc_order(1, 7:11) = [(c + 47, c = 7, 11)]
    call minc_cell_natural_global_test_case(json_str, ' partial bdy', expected_minc_order)
    deallocate(expected_minc_order)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"sws": {"x": [0, 3000], "y": [0, 1000]}, ' // &
         '            "swn": {"x": [0, 3000], "y": [1000, 1500]}},' // &
         '  "minc": {"rock": [{"zones": ["swn"]}, {"zones": ["sws"]}], ' // &
         '           "geometry": {"fracture": {"volume": 0.1}}}},' // &
         '"boundaries": [{"faces": {"cells": [0, 1, 2, 3, 4, 5], ' // &
         '  "normal": [0, -1, 0]}}]' // &
         '}'
    allocate(expected_minc_order(1, 0:48))
    expected_minc_order = 0
    expected_minc_order(1, 0:4) = [(c + 49, c = 0, 4)]
    expected_minc_order(1, 7:11) = [(c + 47, c = 7, 11)]
    call minc_cell_natural_global_test_case(json_str, ' partial multi-rock bdy', expected_minc_order)
    deallocate(expected_minc_order)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"sw": {"x": [0, 3000], "y": [0, 1500]}, ' // &
         '           "ne": {"x": [3000, 4500], "y": [2000, 4500]}},' // &
         '  "minc": [{"rock": {"zones": ["sw"]}, ' // &
         '            "geometry": {"fracture": {"volume": 0.1}}}, ' // &
         '           {"rock": {"zones": ["ne"]}, ' // &
         '            "geometry": {"matrix": {"volume": [0.3, 0.6]}}}]}, ' // &
         '"boundaries": [{"faces": {"cells": [0, 1, 2, 3, 4, 5], ' // &
         '  "normal": [0, -1, 0]}}]' // &
         '}'
    allocate(expected_minc_order(2, 0:48))
    expected_minc_order = 0
    natural = [0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 26, 27, 33, 34, 40, 41, 47, 48]
    minc_natural = [(c, c = 49, 66)]
    do i = 1, size(natural)
       expected_minc_order(1, natural(i)) = minc_natural(i)
    end do
    natural = [26, 27, 33, 34, 40, 41, 47, 48]
    minc_natural = [(c, c = 67, 74)]
    do i = 1, size(natural)
       expected_minc_order(2, natural(i)) = minc_natural(i)
    end do
    call minc_cell_natural_global_test_case(json_str, ' multizone bdy', expected_minc_order)
    deallocate(expected_minc_order, natural, minc_natural)

  contains

    subroutine minc_cell_natural_global_test_case(json_str, title, expected_minc_order)

      character(*), intent(in) :: json_str
      character(*), intent(in) :: title
      PetscInt, intent(in) :: expected_minc_order(:, 0:)
      ! Locals:
      type(mesh_type) :: mesh
      type(IAPWS_type) :: thermo
      type(eos_we_type) :: eos
      type(fson_value), pointer :: json
      PetscInt :: c, p, m, natural, minc_natural, expected_natural
      PetscInt, allocatable :: ic(:)
      PetscInt :: iminc
      DMLabel :: minc_label
      ISLocalToGlobalMapping :: l2g
      character(32) :: natural_str
      PetscErrorCode :: err, ierr
      PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]

      json => fson_parse_mpi(str = json_str)
      call thermo%init()
      call eos%init(json, thermo)
      call mesh%init(eos, json)
      call mesh%configure(gravity, json, err = err)
      call mesh%construct_ghost_cells(gravity)

      call fson_destroy_mpi(json)
      call mesh%destroy_distribution_data()

      call mesh_geometry_sanity_check(mesh, test, title)

      call DMGetLabel(mesh%dm, minc_zone_label_name, minc_label, ierr)
      CHKERRQ(ierr)
      call DMGetLocalToGlobalMapping(mesh%dm, l2g, ierr); CHKERRQ(ierr)

      allocate(ic(maxval(mesh%minc%num_levels)))
      ic = 0
      do c = mesh%strata(0)%start, mesh%strata(0)%end_interior - 1
         if (mesh%ghost_cell(c) < 0) then
            call DMLabelGetValue(minc_label, c, iminc, ierr); CHKERRQ(ierr)
            if (iminc > 0) then
               associate(minc => mesh%minc(iminc))
                 do m = 1, minc%num_levels
                    p = mesh%strata(0)%minc_point(ic(m), m)
                    natural = local_to_natural_cell_index(mesh%cell_natural_global, l2g, c)
                    write(natural_str, '(i6)') natural
                    minc_natural = local_to_natural_cell_index(mesh%cell_natural_global, l2g, p)
                    expected_natural = expected_minc_order(m, natural)
                    call test%assert(expected_natural, minc_natural, &
                         'minc natural ' // title // ' ' // trim(natural_str))
                    ic(m) = ic(m) + 1
                 end do
               end associate
            end if
         end if
      end do

      call mesh%destroy()
      call eos%destroy()
      call thermo%destroy()
      deallocate(ic)

    end subroutine minc_cell_natural_global_test_case

  end subroutine test_minc_cell_natural_global

!------------------------------------------------------------------------

  subroutine test_global_to_fracture_natural(test)
    ! global_to_fracture_natural()

    use fson_mpi_module
    use IAPWS_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    character(:), allocatable :: json_str

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo"}}'
    call global_to_fracture_natural_test_case(test, json_str, 'single porosity', &
         [0, 10, 46], [0, 10, 46], [0, 0, 0])

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"rock": {"zones": ["all"]}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}}}}}'
    call global_to_fracture_natural_test_case(test, json_str, 'MINC all no bdy', &
         [0, 10, 46, 49, 59, 95], [0, 10, 46, 0, 10, 46], [0, 0, 0, 1, 1, 1])
    
    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"rock": {"zones": ["all"]}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}}}}, ' // &
         '"boundaries": [{"faces": {"cells": [0, 1, 2, 3, 4, 5], ' // &
         '  "normal": [0, -1, 0]}}]' // &
         '}'
    call global_to_fracture_natural_test_case(test, json_str, 'MINC all bdy', &
         [0, 10, 46, 49, 59, 95], [0, 10, 46, 0, 10, 46], [0, 0, 0, 1, 1, 1])

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"sw": {"x": [0, 3000], "y": [0, 1500]}},' // &
         '  "minc": {"rock": {"zones": ["sw"]}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}}}}}'
    call global_to_fracture_natural_test_case(test, json_str, 'MINC partial no bdy', &
         [0, 10, 46, 49, 58], [0, 10, 46, 0, 11], [0, 0, 0, 1, 1])

  contains

    subroutine global_to_fracture_natural_test_case(test, json_str, title, &
         natural_indices, expected_fracture_natural_indices, &
         expected_minc_levels)

      class(unit_test_type), intent(in out) :: test
      character(*), intent(in) :: json_str
      character(*), intent(in) :: title
      PetscInt, intent(in) :: natural_indices(:)
      PetscInt, intent(in) :: expected_fracture_natural_indices(:)
      PetscInt, intent(in) :: expected_minc_levels(:)
      ! Locals:
      type(mesh_type) :: mesh
      type(IAPWS_type) :: thermo
      type(eos_we_type) :: eos
      type(fson_value), pointer :: json
      PetscInt :: i, idx(1), natural, global, minc_level
      PetscErrorCode :: err, ierr
      character(2) :: natural_str
      PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]

      json => fson_parse_mpi(str = json_str)
      call thermo%init()
      call eos%init(json, thermo)
      call mesh%init(eos, json)
      call mesh%configure(gravity, json, err = err)
      call mesh%construct_ghost_cells(gravity)

      call fson_destroy_mpi(json)
      call mesh%destroy_distribution_data()

      call mesh_geometry_sanity_check(mesh, test, title)

      do i = 1, size(natural_indices)
         idx(1) = natural_indices(i)
         write(natural_str, '(i2)') idx(1)
         call AOApplicationToPetsc(mesh%cell_natural_global, 1, idx, ierr); CHKERRQ(ierr)
         global = idx(1)
         call mesh%global_to_parent_natural(global, natural, minc_level)
         call test%assert(expected_fracture_natural_indices(i), natural, &
              trim(title) // ' ' // trim(natural_str) // ' natural')
         call test%assert(expected_minc_levels(i), minc_level, trim(title) &
              // ' ' // trim(natural_str) // ' minc level')
      end do

      call mesh%destroy()
      call eos%destroy()
      call thermo%destroy()

    end subroutine global_to_fracture_natural_test_case

  end subroutine test_global_to_fracture_natural

!------------------------------------------------------------------------

  subroutine test_redistribute(test)
    ! MINC redistribute

    use fson_mpi_module
    use IAPWS_module
    use eos_we_module
    use minc_module, only: minc_level_label_name

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    character(:), allocatable :: json_str

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"rock": {"zones": ["all"]}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}}}}}'
    call redistribute_test(test, json_str, 'all', 1, 2 * 49, 1, [49, 49], 0)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"left": {"x": [0, 1500]}},' // &
         '  "minc": {"rock": {"zones": ["left"]}, ' // &
         '           "geometry": {"matrix": {"volume": 0.9}}}}}'
    call redistribute_test(test, json_str, 'partial', 1, 49 + 14, 1, [49, 14], 0)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"rock": {"zones": ["all"]}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}}}}, ' // &
         '"boundaries": [{"faces": {"cells": [0, 1, 2, 3, 4, 5, 6], ' // &
         '  "normal": [0, -1, 0]}}]' // &
         '}'
    call redistribute_test(test, json_str, 'all bdy', 1, 2 * 49, 1, [49, 49], 7)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo",' // &
         '  "zones": {"sw": {"x": [0, 3000], "y": [0, 1500]}, ' // &
         '           "ne": {"x": [3000, 4500], "y": [2000, 4500]}},' // &
         '  "minc": [{"rock": {"zones": ["sw"]}, ' // &
         '            "geometry": {"fracture": {"volume": 0.1}}}, ' // &
         '           {"rock": {"zones": ["ne"]}, ' // &
         '            "geometry": {"matrix": {"volume": [0.3, 0.6]}}}]}, ' // &
         '"boundaries": [{"faces": {"cells": [0, 1, 2, 3, 4, 5, 6], ' // &
         '  "normal": [0, -1, 0]}}]' // &
         '}'
    call redistribute_test(test, json_str, 'multizone bdy', 2, 49 + 10 + 16, &
         2, [49, 10 + 8, 8], 7)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col10.exo",' // &
         '  "zones": {"all": {"-": null}},' // &
         '  "minc": {"rock": {"zones": ["all"]}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}, ' // &
         '"matrix": {"volume": 0.9}}}}, ' // &
         '}'
    call redistribute_test(test, json_str, 'col no bdy', 1, 20, 1, [10, 10], 0)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col10.exo",' // &
         '  "zones": {"minc": {"z": [-600, -100]}},' // &
         '  "minc": {"rock": {"zones": ["minc"]}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}, ' // &
         '"matrix": {"volume": [0.3, 0.6]}}}}, ' // &
         '"boundaries": [{"faces": {"cells": [0], ' // &
         '  "normal": [0, 0, 1]}}]' // &
         '}'
    call redistribute_test(test, json_str, 'col bdy', 1, 20, 2, [10, 5, 5], 1)

  contains

    subroutine redistribute_test(test, json_str, title, expected_num_zones, &
         expected_num_cells, expected_max_level, expected_num_minc_level_cells, &
         expected_num_bdy_cells)

      class(unit_test_type), intent(in out) :: test
      character(*), intent(in) :: json_str
      character(*), intent(in) :: title
      PetscInt, intent(in) :: expected_num_zones, expected_num_cells
      PetscInt, intent(in) :: expected_max_level, expected_num_bdy_cells
      PetscInt, intent(in) :: expected_num_minc_level_cells(0: expected_max_level)
      ! Locals:
      type(mesh_type) :: mesh
      type(IAPWS_type) :: thermo
      type(eos_we_type) :: eos
      type(fson_value), pointer :: json
      PetscSF :: sf
      PetscMPIInt :: rank, np
      PetscInt :: num_minc_zones, max_num_levels, num_cells, num_cells_sf
      PetscInt :: m, num, num_local_bdy_cells, num_bdy_cells
      PetscErrorCode :: err, ierr
      character(48) :: str
      PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]

      call mpi_comm_rank(PETSC_COMM_WORLD, rank, ierr)
      call mpi_comm_size(PETSC_COMM_WORLD, np, ierr)

      json => fson_parse_mpi(str = json_str)
      call thermo%init()
      call eos%init(json, thermo)
      call mesh%init(eos, json)
      call mesh%configure(gravity, json, err = err)

      if (np > 1) then
         call mesh%redistribute(sf)
         call test%assert(sf .ne. PETSC_NULL_SF, title // ": null redistribution SF")
      end if
      call mesh%construct_ghost_cells(gravity)
      call fson_destroy_mpi(json)
      call mesh%destroy_distribution_data()

      call mesh_geometry_sanity_check(mesh, test, title)
      call dm_dag_sanity_check(mesh%dm, test, title)
      call minc_mesh_natural_check(mesh, test, title)

      if (rank == 0) then
         call test%assert(mesh%has_minc, title // ": mesh has minc")
      end if
      num_minc_zones = size(mesh%minc)
      max_num_levels = maxval(mesh%minc%num_levels)
      if (rank == 0) then
         call test%assert(expected_num_zones, &
              num_minc_zones, title // ": num minc zones")
         call test%assert(expected_max_level, &
              max_num_levels, title // ": num minc levels")
      end if

      num_cells = total_interior_cell_count(mesh)
      num_cells_sf = total_interior_cell_count_sf(mesh)
      if (rank == 0) then
         call test%assert(expected_num_cells, &
              num_cells, title  // ": num cells")
         call test%assert(expected_num_cells, &
              num_cells_sf, title  // ": num cells SF")
      end if

      do m = 0, expected_max_level
         num = total_interior_cell_label_count(mesh, &
              minc_level_label_name, m)
         if (rank == 0) then
            write(str, '(a, i1)') ": num minc points, level ", m
            call test%assert(expected_num_minc_level_cells(m), &
                 num, title  // str)
         end if
      end do

      call DMGetStratumSize(mesh%dm, boundary_ghost_label_name, 1, &
           num_local_bdy_cells, ierr); CHKERRQ(ierr)
      call mpi_reduce(num_local_bdy_cells, num_bdy_cells, 1, &
              MPI_INTEGER, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
      if (rank == 0) then
         call test%assert(expected_num_bdy_cells, num_bdy_cells, &
              title // " num boundary cells")
      end if
      if (np > 1) then
         call PetscSFDestroy(sf, ierr); CHKERRQ(ierr)
      end if
      call mesh%destroy()
      call eos%destroy()
      call thermo%destroy()

    end subroutine redistribute_test

  end subroutine test_redistribute

!------------------------------------------------------------------------

end module mesh_test
