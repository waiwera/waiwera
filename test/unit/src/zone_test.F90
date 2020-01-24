module zone_test

  ! Tests for zone module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use fson
  use zone_module
  use zone_label_module
  use fson_mpi_module
  use list_module

  implicit none
  private

  character(len = 512) :: data_path

  public :: setup, teardown
  public :: test_get_zone_type, test_cell_array, test_zone_depends, &
       test_cell_array_label, test_box_label, test_combine_label

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

  subroutine test_get_zone_type(test)

    ! get_zone_type

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    PetscInt :: zone_type
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '[1, 2, 3]')
    zone_type = get_zone_type_mpi(json)
    if (rank == 0) then
       call test%assert(ZONE_TYPE_CELL_ARRAY, zone_type, "array")
    end if
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"type": "array", "cells": [1, 2, 3]}')
    zone_type = get_zone_type_mpi(json)
    if (rank == 0) then
       call test%assert(ZONE_TYPE_CELL_ARRAY, zone_type, "type array")
    end if
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"cells": [1, 2, 3]}')
    zone_type = get_zone_type_mpi(json)
    if (rank == 0) then
       call test%assert(ZONE_TYPE_CELL_ARRAY, zone_type, "cells")
    end if
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"type": "foo"}')
    zone_type = get_zone_type_mpi(json)
    if (rank == 0) then
       call test%assert(-1, zone_type, "unknown type")
    end if
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"type": "box"}')
    zone_type = get_zone_type_mpi(json)
    if (rank == 0) then
       call test%assert(ZONE_TYPE_BOX, zone_type, "type box")
    end if
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"x": [0, 100]}')
    zone_type = get_zone_type_mpi(json)
    if (rank == 0) then
       call test%assert(ZONE_TYPE_BOX, zone_type, "x box")
    end if
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"r": [100, 200]}')
    zone_type = get_zone_type_mpi(json)
    if (rank == 0) then
       call test%assert(ZONE_TYPE_BOX, zone_type, "r box")
    end if
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"y": [100, 200], "z": [-50, 50]}')
    zone_type = get_zone_type_mpi(json)
    if (rank == 0) then
       call test%assert(ZONE_TYPE_BOX, zone_type, "y-z box")
    end if
    call fson_destroy_mpi(json)

  end subroutine test_get_zone_type

!------------------------------------------------------------------------

  subroutine test_cell_array(test)
    ! cell array

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(zone_cell_array_type) :: zone
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '[1, 2, 3]')
    if (rank == 0) then
       call zone%init_serial(1, 'zone1', json)
       call test%assert([1, 2, 3], zone%cells, 'array cells')
       call zone%destroy()
    end if
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"cells": [1, 2, 3]}')
    if (rank == 0) then
       call zone%init_serial(1, 'zone1', json)
       call test%assert([1, 2, 3], zone%cells, 'cells cells')
       call zone%destroy()
    end if
    call fson_destroy_mpi(json)

  end subroutine test_cell_array

!------------------------------------------------------------------------

  subroutine test_zone_depends(test)
    ! dependencies

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    class(zone_type), pointer :: zone
    character(max_zone_name_length), allocatable :: expected_depends(:)
    character(max_zone_name_length), allocatable :: depends(:)
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '[1, 2, 3]')
    expected_depends = [character(max_zone_name_length)::]
    allocate(zone_cell_array_type:: zone)
    call zone%init(1, 'zone', json)
    call zone%dependencies(depends)
    if (rank == 0) then
       call test%assert(expected_depends, depends, 'cells')
    end if
    call zone%destroy()
    deallocate(zone)
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"x": [0, 100]}')
    expected_depends = [character(max_zone_name_length)::]
    allocate(zone_box_type:: zone)
    call zone%init(1, 'zone', json)
    call zone%dependencies(depends)
    if (rank == 0) then
       call test%assert(expected_depends, depends, 'box')
    end if
    call zone%destroy()
    deallocate(zone)
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"+": ["zone1", "zone2"], ' // &
         '"*": "zone3", "-": "zone4"}')
    expected_depends = ["zone1", "zone2", "zone3", "zone4"]
    allocate(zone_combine_type:: zone)
    call zone%init(1, 'zone', json)
    call zone%dependencies(depends)
    if (rank == 0) then
       call test%assert(expected_depends, depends, 'combine')
    end if
    call zone%destroy()
    deallocate(zone)
    call fson_destroy_mpi(json)

  end subroutine test_zone_depends

!------------------------------------------------------------------------

  subroutine test_cell_array_label(test)
    ! cell array label

    use mesh_module
    use IAPWS_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    PetscMPIInt :: rank
    PetscErrorCode :: ierr, err

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call thermo%init()

    json => fson_parse_mpi(str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo", ' // &
         '"zones": {"zone1": [10, 15, 20, 27, 34, 44], ' // &
         '"zone2": [40, 30, 5]}}}')
    call eos%init(json, thermo)
    call mesh%init(eos, json)
    call mesh%configure(gravity, json, err = err)

    call test%assert(0, err, 'config error')
    call fson_destroy_mpi(json)
    call mesh%destroy_distribution_data()

    if (err == 0) then

       if (rank == 0) then
          call test%assert(2, mesh%zones%count, 'num zones')
       end if

       call zone_test(0, [10, 15, 20, 27, 34, 44])
       call zone_test(1, [40, 30, 5])

    end if

    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()

  contains

    subroutine zone_test(index, cells)

      PetscInt, intent(in) :: index
      PetscInt, intent(in) :: cells(:)
      ! Locals:
      type(list_node_type), pointer :: node
      PetscInt :: num_found, num_found_local, i
      PetscInt :: num_found_non_ghost, c, ghost
      DMLabel :: ghost_label
      PetscInt, pointer, contiguous :: points_array(:)
      IS :: points
      character(40) :: istr

      write(istr, '(a,i1,a)') '[', index, ']'

      associate(num_cells => size(cells))

        node => mesh%zones%get(index)
        select type(zone => node%data)
        type is (zone_cell_array_type)

           num_found_non_ghost = 0
           call DMGetStratumSize(mesh%dm, zone_label_name(zone%name), 1, &
                num_found_local, ierr); CHKERRQ(ierr)
           if (num_found_local > 0) then
              call DMGetStratumIS(mesh%dm, zone_label_name(zone%name), 1, &
                   points, ierr); CHKERRQ(ierr)
              call ISGetIndicesF90(points, points_array, ierr); CHKERRQ(ierr)
              call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
              do i = 1, num_found_local
                 c = points_array(i)
                 call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
                 if (ghost < 0) then
                    num_found_non_ghost = num_found_non_ghost + 1
                 end if
              end do
              call ISRestoreIndicesF90(points, points_array, ierr); CHKERRQ(ierr)
           end if
           call MPI_reduce(num_found_non_ghost, num_found, 1, MPI_INTEGER, &
                MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
           if (rank == 0) then
              call test%assert(num_cells, num_found, 'num found ' // istr)
           end if

        end select
      end associate

    end subroutine zone_test

  end subroutine test_cell_array_label

!------------------------------------------------------------------------

  subroutine test_box_label(test)
    ! box label

    use IAPWS_module
    use eos_we_module
    use mesh_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(mesh_type) :: mesh
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscMPIInt :: rank
    PetscErrorCode :: ierr, err

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call thermo%init()

    json => fson_parse_mpi(str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo", ' // &
         '"zones": {' // &
         '"xzone": {"x": [2000, 3000]}, ' // &
         '"all": {"type": "box"}, ' // &
         '"xyzone": {"x": [0, 2000], "y": [2500, 4500]}}}}')
    call eos%init(json, thermo)
    call mesh%init(eos, json)
    call mesh%configure(gravity, json, err = err)

    call test%assert(0, err, 'config error')
    call fson_destroy_mpi(json)
    call mesh%destroy_distribution_data()

    if (err == 0) then
       call mesh%zones%traverse(box_label_iterator)
    end if

    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()

  contains

    subroutine box_label_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: num_found, num_found_local, num_expected
      character(:), allocatable :: label_name

      select type (zone => node%data)
      class is (zone_type)
         label_name = zone_label_name(zone%name)
         call DMGetStratumSize(mesh%dm, label_name, 1, &
              num_found_local, ierr); CHKERRQ(ierr)
         call MPI_reduce(num_found_local, num_found, 1, MPI_INTEGER, &
              MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
         if (rank == 0) then
            select case (zone%name)
            case ('xzone')
               num_expected = 14
            case ('all')
               num_expected = 49
            case ('xyzone')
               num_expected = 9
            end select
            call test%assert(num_expected, num_found, &
                 'num found ' //  zone%name)
         end if
      end select
      stopped = PETSC_FALSE

    end subroutine box_label_iterator

  end subroutine test_box_label

!------------------------------------------------------------------------

  subroutine test_combine_label(test)
    ! combine label

    use mesh_module
    use IAPWS_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    PetscMPIInt :: rank
    PetscErrorCode :: ierr, err

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call thermo%init()

    json => fson_parse_mpi(str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/7x7grid.exo", ' // &
         '"zones": {' // &
         '"zone1": {"x": [2000, 3000]}, ' // &
         '"zone2": {"x": [3500, 4500]}, ' // &
         '"zone3": {"x": [2500, 4500], "y": [0, 1000]}, ' // &
         '"zone_plus": {"+": ["zone1", "zone2"]}, ' // &
         '"zone_minus": {"+": "zone_plus", "-": "zone3"}, ' // &
         '"zone_times": {"+": "zone_plus", "*": "zone3"}, ' // &
         '"zone_times2": {"*": ["zone_plus", "zone3"]}, ' // &
         '"all": {"-": null}}}}')
    call eos%init(json, thermo)
    call mesh%init(eos, json)
    call mesh%configure(gravity, json, err = err)

    call test%assert(0, err, 'config error')
    call fson_destroy_mpi(json)
    call mesh%destroy_distribution_data()

    if (err == 0) then
       call mesh%zones%traverse(combine_label_iterator)
    end if

    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()

  contains

    subroutine combine_label_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: num_found, num_found_local, num_expected
      character(:), allocatable :: label_name

      select type (zone => node%data)
      class is (zone_type)
         label_name = zone_label_name(zone%name)
         call DMGetStratumSize(mesh%dm, label_name, 1, &
              num_found_local, ierr); CHKERRQ(ierr)
         call MPI_reduce(num_found_local, num_found, 1, MPI_INTEGER, &
              MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
         if (rank == 0) then
            select case (zone%name)
            case ('zone1')
               num_expected = 14
            case ('zone2')
               num_expected = 7
            case ('zone3')
               num_expected = 3
            case ('zone_plus')
               num_expected = 21
            case ('zone_minus')
               num_expected = 19
            case ('zone_times', 'zone_times2')
               num_expected = 2
            case ('all')
               num_expected = 49
            end select
            call test%assert(num_expected, num_found, &
                 'num found: ' //  zone%name)
         end if
      end select
      stopped = PETSC_FALSE

    end subroutine combine_label_iterator

  end subroutine test_combine_label

!------------------------------------------------------------------------

end module zone_test
