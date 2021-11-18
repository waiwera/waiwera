module source_setup_test

  ! Tests for source setup module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use source_setup_module
  use eos_wge_module

  implicit none
  private

  character(len = 512) :: data_path

  public :: setup, teardown
  public :: test_setup_sources, test_source_index

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

  subroutine test_setup_sources(test)

    ! setup_sources() test

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use source_module
    use IAPWS_module
    use eos_module
    use dm_utils_module, only: global_vec_section, global_section_offset
    use utils_module, only: array_cumulative_sum, get_mpi_int_gather_array
    use tracer_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IAPWS_type) :: thermo
    type(eos_wge_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(source_type) :: source
    Vec :: fluid_vector, source_vector
    PetscReal, pointer, contiguous :: source_array(:)
    PetscSection :: source_section
    PetscInt :: fluid_range_start, source_range_start
    PetscInt :: s, source_offset, source_index
    type(list_type) :: source_controls
    PetscInt :: num_sources, total_num_sources, num_zone_sources, n_all, i
    PetscErrorCode :: ierr, err
    PetscReal, parameter :: start_time = 0._dp
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    type(tracer_type), allocatable :: tracers(:)
    PetscInt, parameter :: expected_num_sources = 23
    PetscMPIInt :: rank, num_procs
    IS :: source_is
    PetscInt, allocatable :: zone_source(:), isort(:)
    PetscInt, allocatable :: zone_source_sorted(:), zone_source_all(:)
    PetscInt, allocatable :: zone_source_counts(:), zone_source_displacements(:)
    PetscInt, parameter :: expected_zone_source_cells(3) = [0, 4, 8]

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)
    json => fson_parse_mpi(trim(adjustl(data_path)) // "source/test_source.json")

    call thermo%init()
    call eos%init(json, thermo)
    call setup_tracers(json, eos, tracers, err = err)
    call source%init(eos, size(tracers))
    call mesh%init(eos, json)
    call DMCreateLabel(mesh%serial_dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(gravity, json, err = err)
    call DMGetGlobalVector(mesh%dm, fluid_vector, ierr); CHKERRQ(ierr) ! dummy- not used

    call setup_sources(json, mesh%dm, mesh%cell_natural_global, eos, tracers%name, &
         thermo, start_time, fluid_vector, fluid_range_start, source_vector, &
         source_range_start, num_sources, total_num_sources, source_controls, &
         source_is, err = err)
    call test%assert(0, err, "error")

    if (rank == 0) then
      call test%assert(expected_num_sources, total_num_sources, "number of sources")
    end if

    call global_vec_section(source_vector, source_section)
    call VecGetArrayReadF90(source_vector, source_array, ierr); CHKERRQ(ierr)
    allocate(zone_source(num_sources))
    zone_source = -1
    num_zone_sources = 0

    do s = 0, num_sources - 1
       source_offset = global_section_offset(source_section, s, source_range_start)
       call source%assign(source_array, source_offset)
       source_index = nint(source%source_index)
       select case (source_index)
       case (0)
            call source_test(test, source_index, source, &
                 0, 10._dp, 90.e3_dp, 0, 0, [0._dp, 0._dp])
         case (1)
            call source_test(test, source_index, source, &
                 1, 5._dp, 100.e3_dp, 2, 0)
         case (2)
            call source_test(test, source_index, source, &
                 2, 1000._dp, 0._dp, 3, 3)
         case (3)
            call source_test(test, source_index, source, &
                 3, -2._dp, default_source_injection_enthalpy, 1, 0)
         case (4)
            call source_test(test, source_index, source, &
                 4, -3._dp, 200.e3_dp, 1, 0)
         case (5)
            call source_test(test, source_index, source, &
                 5, -5._dp, default_source_injection_enthalpy, 0, 0)
         case (6)
            call source_test(test, source_index, source, &
                 6, -2000._dp, 0._dp, 3, 3)
         case (7)
            call source_test(test, source_index, source, &
                 7, default_source_rate, default_source_injection_enthalpy, 1, 0)
         case (8)
            call source_test(test, source_index, source, &
                 8, default_source_rate, 1000.e3_dp, 2, 0)
         case (9)
            call source_test(test, source_index, source, &
                 0, default_source_rate, 0._dp, 3, 3)
         case (10)
            call source_test(test, source_index, source, &
                 1, 3._dp, 150.e3_dp, 1, 1)
         case (11)
            call source_test(test, source_index, source, &
                 2, default_source_rate, default_source_injection_enthalpy, 1, 1)
         case (12)
            call source_test(test, source_index, source, &
                 3, default_source_rate, 80.e3_dp, 2, 2)
         case (13)
            call source_test(test, source_index, source, &
                 4, default_source_rate, 90.e3_dp, 1, 2)
         case (14)
            call source_test(test, source_index, source, &
                 5, default_source_rate, 500.e3_dp, 2, 3)
         case (15)
            call source_test(test, source_index, source, &
                 6, default_source_rate, 100.e3_dp, default_source_component, 2)
         case (16)
            call source_test(test, source_index, source, &
                 1, 10._dp, 50.e3_dp, default_source_component, 0, &
                 [1.e-3_dp, 1.e-3_dp])
         case (17)
            call source_test(test, source_index, source, &
                 2, 5._dp, 40.e3_dp, default_source_component, 0, &
                 [2.e-3_dp, 3.e-3_dp])
         case (18)
            call source_test(test, source_index, source, &
                 3, 7.5_dp, 30.e3_dp, default_source_component, 0, &
                 [3.e-3_dp, 5.e-3_dp])
         case (19)
            call source_test(test, source_index, source, &
                 4, 3.5_dp, 60.e3_dp, default_source_component, 0, &
                 [0._dp, 2.e-3_dp])
         case (20:)
            num_zone_sources = num_zone_sources + 1
            zone_source(num_zone_sources) = nint(source%natural_cell_index)
         end select
    end do

    ! Test cells in last source, defined on a zone:
    zone_source = pack(zone_source, zone_source >= 0)
    zone_source_counts = get_mpi_int_gather_array()
    zone_source_displacements = get_mpi_int_gather_array()
    call MPI_gather(num_zone_sources, 1, MPI_INTEGER, zone_source_counts, 1, &
         MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       zone_source_displacements = [[0], &
            array_cumulative_sum(zone_source_counts(1: num_procs - 1))]
       n_all = sum(zone_source_counts)
    else
       n_all = 1
    end if
    allocate(zone_source_all(n_all))
    call MPI_gatherv(zone_source, num_zone_sources, MPI_INTEGER, &
         zone_source_all, zone_source_counts, zone_source_displacements, &
         MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       isort = [(i - 1, i = 1, n_all)]
       call PetscSortIntWithPermutation(n_all, &
            zone_source_all, isort, ierr); CHKERRQ(ierr)
       isort = isort + 1 ! convert to 1-based
       allocate(zone_source_sorted(n_all))
       do i = 1, n_all
          zone_source_sorted(i) = zone_source_all(isort(i))
       end do
       call test%assert(expected_zone_source_cells, zone_source_sorted, &
            "zone source cells")
       deallocate(zone_source_sorted, isort)
    end if

    deallocate(zone_source, zone_source_counts, zone_source_displacements, &
         zone_source_all, tracers)
    call ISDestroy(source_is, ierr); CHKERRQ(ierr)
    call source%destroy()
    call VecRestoreArrayReadF90(source_vector, source_array, ierr); CHKERRQ(ierr)
    call VecDestroy(source_vector, ierr); CHKERRQ(ierr)
    call source_controls%destroy()
    call DMRestoreGlobalVector(mesh%dm, fluid_vector, ierr); CHKERRQ(ierr)
    call mesh%destroy_distribution_data()
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  contains

    subroutine source_test(test, source_index, source, index, rate, enthalpy, &
         injection_component, production_component, tracer)
      !! Runs asserts for a single source.

      class(unit_test_type), intent(in out) :: test
      PetscInt, intent(in) :: source_index
      type(source_type), intent(in) :: source
      PetscInt, intent(in) :: index
      PetscInt, intent(in) :: injection_component, production_component
      PetscReal, intent(in) :: rate, enthalpy
      PetscReal, intent(in), optional :: tracer(:)
      ! Locals:
      character(12) :: srcstr
      PetscReal, parameter :: tol = 1.e-6_dp

      write(srcstr, '(a, i2, a)') 'source[', source_index, ']'
      call test%assert(index, nint(source%natural_cell_index), &
           trim(srcstr) // ": natural index")
      call test%assert(rate, source%rate, &
           trim(srcstr) // ": rate")
      call test%assert(enthalpy, source%injection_enthalpy, &
           trim(srcstr) // ": enthalpy")
      call test%assert(injection_component, nint(source%injection_component), &
           trim(srcstr) // ": injection component")
      call test%assert(production_component, nint(source%production_component), &
           trim(srcstr) // ": production component")
      if (present(tracer)) then
         call test%assert(tracer, source%tracer_injection_rate, &
              trim(srcstr) // ": tracer")
      end if

    end subroutine source_test

  end subroutine test_setup_sources

!------------------------------------------------------------------------

  subroutine test_source_index(test)
    ! test source_index

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use source_module
    use IAPWS_module
    use eos_we_module
    use tracer_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(tracer_type), allocatable :: tracers(:)
    Vec :: fluid_vector, source_vector
    PetscInt :: num_sources, total_num_sources
    PetscInt :: fluid_range_start, source_range_start
    type(list_type) :: source_controls
    IS :: source_index
    PetscInt, pointer, contiguous :: source_index_array(:)
    PetscMPIInt :: rank, num_procs
    PetscErrorCode :: err, ierr
    PetscReal, parameter :: start_time = 0._dp
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscInt, parameter :: expected_num_sources = 6

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)
    json => fson_parse_mpi(trim(adjustl(data_path)) // "source/test_source_index.json")

    call thermo%init()
    call eos%init(json, thermo)
    call setup_tracers(json, eos, tracers, err = err)
    call mesh%init(eos, json)
    call DMCreateLabel(mesh%serial_dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(gravity, json, err = err)
    call DMGetGlobalVector(mesh%dm, fluid_vector, ierr); CHKERRQ(ierr) ! dummy- not used

    call setup_sources(json, mesh%dm, mesh%cell_natural_global, eos, tracers%name, &
         thermo, start_time, fluid_vector, fluid_range_start, source_vector, &
         source_range_start, num_sources, total_num_sources, source_controls, &
         source_index, err = err)
    call test%assert(0, err, "error")

    if (rank == 0) then
      call test%assert(expected_num_sources, total_num_sources, "number of sources")
    end if

    call ISGetIndicesF90(source_index, source_index_array, ierr); CHKERRQ(ierr)
    call test%assert(all(source_index_array >= 0), "indices >= 0")
    call ISRestoreIndicesF90(source_index, source_index_array, ierr); CHKERRQ(ierr)

    call ISDestroy(source_index, ierr); CHKERRQ(ierr)
    call VecDestroy(source_vector, ierr); CHKERRQ(ierr)
    call source_controls%destroy()
    call DMRestoreGlobalVector(mesh%dm, fluid_vector, ierr); CHKERRQ(ierr)
    call mesh%destroy_distribution_data()
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_source_index

!------------------------------------------------------------------------

end module source_setup_test
