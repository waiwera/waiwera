module source_setup_test

  ! Tests for source setup module

  use kinds_module
  use mpi_module
  use fruit
  use source_setup_module
  use eos_test

  implicit none
  private

#include <petsc/finclude/petsc.h90>

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

    character(16), parameter :: path = "data/source/"
    type(IAPWS_type) :: thermo
    type(eos_test_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(list_type) :: sources, source_controls
    PetscInt :: expected_num_sources, num_sources
    PetscErrorCode :: ierr

    json => fson_parse_mpi(trim(path) // "test_source.json")

    call thermo%init()
    call eos%init(json, thermo)
    call mesh%init(json)
    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(eos%primary_variable_names)

    call setup_sources(json, mesh%dm, eos, sources, source_controls)

    expected_num_sources = fson_value_count_mpi(json, "source")
    call MPI_reduce(sources%count, num_sources, 1, MPI_INTEGER, MPI_SUM, &
         mpi%input_rank, mpi%comm, ierr)
    if (mpi%rank == mpi%input_rank) then
      call assert_equals(expected_num_sources, num_sources, "number of sources")
    end if

    call sources%traverse(source_test_iterator)

    call sources%destroy(source_list_node_data_destroy)
    call source_controls%destroy()
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  contains

    subroutine source_test(tag, source, index, rate, enthalpy, &
         injection_component, production_component)
      !! Runs asserts for a single source.
      character(*), intent(in) :: tag
      type(source_type), intent(in) :: source
      PetscInt, intent(in) :: index
      PetscInt, intent(in) :: injection_component, production_component
      PetscReal, intent(in) :: rate, enthalpy
      ! Locals:
      PetscReal, parameter :: tol = 1.e-6_dp

      call assert_equals(index, source%cell_natural_index, &
           trim(tag) // ": natural index")
      call assert_equals(rate, source%rate, tol, &
           trim(tag) // ": rate")
      call assert_equals(enthalpy, source%injection_enthalpy, tol, &
           trim(tag) // ": enthalpy")
      call assert_equals(injection_component, source%injection_component, &
           trim(tag) // ": injection component")
      call assert_equals(production_component, source%production_component, &
           trim(tag) // ": production component")

    end subroutine source_test

    subroutine source_test_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      select type (source => node%data)
      type is (source_type)
         select case (node%tag)
         case ("mass injection 1")
            call source_test(node%tag, source, &
                 0, 10._dp, 90.e3_dp, 0, 0)
         case ("mass injection 2")
            call source_test(node%tag, source, &
                 1, 5._dp, 100.e3_dp, 2, 0)
         case ("heat injection")
            call source_test(node%tag, source, &
                 2, 1000._dp, 0._dp, 3, 3)
         case ("mass component production")
            call source_test(node%tag, source, &
                 3, -2._dp, 0._dp, 1, 0)
         case ("mass component production enthalpy")
            call source_test(node%tag, source, &
                 4, -3._dp, 200.e3_dp, 1, 0)
         case ("mass production")
            call source_test(node%tag, source, &
                 5, -5._dp, 0._dp, 0, 0)
         case ("heat production")
            call source_test(node%tag, source, &
                 6, -2000._dp, 0._dp, 3, 3)
         case ("no rate mass")
            call source_test(node%tag, source, &
                 7, default_source_rate, default_source_injection_enthalpy, 1, 0)
         case ("no rate mass enthalpy")
            call source_test(node%tag, source, &
                 8, default_source_rate, 1000.e3_dp, 2, 0)
         case ("no rate heat")
            call source_test(node%tag, source, &
                 0, default_source_rate, 0._dp, 3, 3)
         case ("production component 1")
            call source_test(node%tag, source, &
                 1, 3._dp, 150.e3_dp, 1, 1)
         case ("production component 2")
            call source_test(node%tag, source, &
                 2, default_source_rate, default_source_injection_enthalpy, 1, 1)
         case ("production component 3")
            call source_test(node%tag, source, &
                 3, default_source_rate, 80.e3_dp, 2, 2)
         case ("production component 4")
            call source_test(node%tag, source, &
                 4, default_source_rate, 90.e3_dp, 1, 2)
         case ("production component 5")
            call source_test(node%tag, source, &
                 5, default_source_rate, 500.e3_dp, 2, 3)
         case ("production component 6")
            call source_test(node%tag, source, &
                 6, default_source_rate, 100.e3_dp, default_source_component, 2)
         end select
      end select

      stopped = PETSC_FALSE

    end subroutine source_test_iterator

    subroutine source_list_node_data_destroy(node)
      type(list_node_type), pointer, intent(in out) :: node

      select type (source => node%data)
      type is (source_type)
         call source%destroy()
      end select
    end subroutine source_list_node_data_destroy

  end subroutine test_setup_sources

!------------------------------------------------------------------------

end module source_setup_test
