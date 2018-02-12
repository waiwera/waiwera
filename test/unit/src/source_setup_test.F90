module source_setup_test

  ! Tests for source setup module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fruit
  use source_setup_module
  use eos_test

  implicit none
  private

public :: test_setup_sources

contains

!------------------------------------------------------------------------

  subroutine test_setup_sources

    ! setup_sources() test

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use source_module
    use IAPWS_module
    use eos_module
    use dm_utils_module, only: global_vec_section, global_section_offset

    character(16), parameter :: path = "data/source/"
    type(IAPWS_type) :: thermo
    type(eos_test_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(source_type) :: source
    Vec :: fluid_vector, source_vector
    PetscReal, pointer, contiguous :: source_array(:)
    PetscSection :: source_section
    PetscInt :: fluid_range_start, source_range_start
    PetscInt :: s, source_offset, source_index
    type(list_type) :: source_controls
    PetscInt :: expected_num_sources, num_sources, total_num_sources
    PetscErrorCode :: ierr, err
    PetscReal, parameter :: start_time = 0._dp
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscMPIInt :: rank
    PetscViewer :: viewer

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    json => fson_parse_mpi(trim(path) // "test_source.json")
    viewer = PETSC_NULL_VIEWER

    call thermo%init()
    call eos%init(json, thermo)
    call source%init(eos)
    call mesh%init(json)
    call DMCreateLabel(mesh%original_dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(eos, gravity, json, viewer = viewer, err = err)
    call DMGetGlobalVector(mesh%dm, fluid_vector, ierr); CHKERRQ(ierr) ! dummy- not used

    call setup_sources(json, mesh%dm, mesh%cell_order, eos, thermo, start_time, &
         fluid_vector, fluid_range_start, source_vector, source_range_start, &
         num_sources, source_controls, err = err)
    call assert_equals(0, err, "error")

    expected_num_sources = fson_value_count_mpi(json, "source")
    call MPI_reduce(num_sources, total_num_sources, 1, MPI_INTEGER, MPI_SUM, &
         0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
      call assert_equals(expected_num_sources, total_num_sources, "number of sources")
    end if

    call global_vec_section(source_vector, source_section)
    call VecGetArrayReadF90(source_vector, source_array, ierr); CHKERRQ(ierr)

    do s = 0, num_sources - 1
       call global_section_offset(source_section, s, &
            source_range_start, source_offset, ierr); CHKERRQ(ierr)
       call source%assign(source_array, source_offset)
       source_index = nint(source%source_index)
       select case (source_index)
       case (0)
            call source_test(source_index, source, &
                 0, 10._dp, 90.e3_dp, 0, 0)
         case (1)
            call source_test(source_index, source, &
                 1, 5._dp, 100.e3_dp, 2, 0)
         case (2)
            call source_test(source_index, source, &
                 2, 1000._dp, 0._dp, 3, 3)
         case (3)
            call source_test(source_index, source, &
                 3, -2._dp, default_source_injection_enthalpy, 1, 0)
         case (4)
            call source_test(source_index, source, &
                 4, -3._dp, 200.e3_dp, 1, 0)
         case (5)
            call source_test(source_index, source, &
                 5, -5._dp, default_source_injection_enthalpy, 0, 0)
         case (6)
            call source_test(source_index, source, &
                 6, -2000._dp, 0._dp, 3, 3)
         case (7)
            call source_test(source_index, source, &
                 7, default_source_rate, default_source_injection_enthalpy, 1, 0)
         case (8)
            call source_test(source_index, source, &
                 8, default_source_rate, 1000.e3_dp, 2, 0)
         case (9)
            call source_test(source_index, source, &
                 0, default_source_rate, 0._dp, 3, 3)
         case (10)
            call source_test(source_index, source, &
                 1, 3._dp, 150.e3_dp, 1, 1)
         case (11)
            call source_test(source_index, source, &
                 2, default_source_rate, default_source_injection_enthalpy, 1, 1)
         case (12)
            call source_test(source_index, source, &
                 3, default_source_rate, 80.e3_dp, 2, 2)
         case (13)
            call source_test(source_index, source, &
                 4, default_source_rate, 90.e3_dp, 1, 2)
         case (14)
            call source_test(source_index, source, &
                 5, default_source_rate, 500.e3_dp, 2, 3)
         case (15)
            call source_test(source_index, source, &
                 6, default_source_rate, 100.e3_dp, default_source_component, 2)
         end select
    end do

    call source%destroy()
    call VecRestoreArrayReadF90(source_vector, source_array, ierr); CHKERRQ(ierr)
    call VecDestroy(source_vector, ierr); CHKERRQ(ierr)
    call source_controls%destroy()
    call DMRestoreGlobalVector(mesh%dm, fluid_vector, ierr); CHKERRQ(ierr)
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  contains

    subroutine source_test(source_index, source, index, rate, enthalpy, &
         injection_component, production_component)
      !! Runs asserts for a single source.
      PetscInt, intent(in) :: source_index
      type(source_type), intent(in) :: source
      PetscInt, intent(in) :: index
      PetscInt, intent(in) :: injection_component, production_component
      PetscReal, intent(in) :: rate, enthalpy
      ! Locals:
      character(12) :: srcstr
      PetscReal, parameter :: tol = 1.e-6_dp

      write(srcstr, '(a, i2, a)') 'source[', source_index, ']'
      call assert_equals(index, nint(source%natural_cell_index), &
           trim(srcstr) // ": natural index")
      call assert_equals(rate, source%rate, tol, &
           trim(srcstr) // ": rate")
      call assert_equals(enthalpy, source%injection_enthalpy, tol, &
           trim(srcstr) // ": enthalpy")
      call assert_equals(injection_component, nint(source%injection_component), &
           trim(srcstr) // ": injection component")
      call assert_equals(production_component, nint(source%production_component), &
           trim(srcstr) // ": production component")

    end subroutine source_test

  end subroutine test_setup_sources

!------------------------------------------------------------------------

end module source_setup_test
