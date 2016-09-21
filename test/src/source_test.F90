module source_test

  ! Test for source module

  use kinds_module
  use mpi_module
  use fruit
  use source_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

public :: test_setup_sources, test_source_update_flow

contains

!------------------------------------------------------------------------

  subroutine test_setup_sources

    ! setup_sources() test

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use eos_module, only: max_primary_variable_name_length

    character(16), parameter :: path = "data/source/"
    PetscInt, parameter :: nc = 2
    PetscInt :: np
    character(max_primary_variable_name_length), allocatable :: &
         primary_variable_names(:)
    PetscBool :: isothermal
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(list_type) :: sources
    PetscErrorCode :: ierr

    primary_variable_names = &
         [character(max_primary_variable_name_length) :: &
         "P1", "P2", "T"]
    np = size(primary_variable_names)

    isothermal = (nc == np)

    json => fson_parse_mpi(trim(path) // "test_source.json")

    call mesh%init(json)
    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(primary_variable_names)

    call setup_sources(json, mesh%dm, np, isothermal, sources)
    call sources%traverse(source_test_iterator)

    call sources%destroy(source_list_node_data_destroy)
    call mesh%destroy()
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

  subroutine test_source_update_flow
    ! update_flow() test

    use fluid_module, only: fluid_type

    type(source_type) :: source
    type(fluid_type) :: fluid
    PetscBool :: isothermal = PETSC_FALSE
    PetscInt, parameter :: num_components = 2, num_phases = 2
    PetscInt, parameter :: num_primary = num_components + 1
    PetscInt, parameter :: offset = 1
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal, parameter :: tol = 1.e-6_dp

    if (mpi%rank == mpi%output_rank) then

       call fluid%init(num_components, num_phases)

       allocate(fluid_data(offset - 1 + fluid%dof))
       fluid_data = [2.7e5_dp, 130._dp, 4._dp, 3._dp, &
            935._dp, 1.e-6_dp, 0.8_dp, 0.7_dp, 83.9e3_dp, 5.461e5_dp, 0.7_dp, 0.3_dp, &
            1.5_dp,  2.e-7_dp, 0.2_dp, 0.3_dp, 800.e3_dp, 2.540e6_dp, 0.4_dp, 0.6_dp]

       call fluid%assign(fluid_data, offset)

       call source_flow_test("inject 1", 10._dp, 200.e3_dp, 1, 0, &
            [10._dp, 0._dp, 2.e6_dp])
       call source_flow_test("inject 2", 5._dp, 200.e3_dp, 2, 0, &
            [0._dp, 5._dp, 1.e6_dp])
       call source_flow_test("inject 3", 5._dp, 200.e3_dp, 2, 2, &
            [0._dp, 5._dp, 1.e6_dp])
       call source_flow_test("inject heat", 1000._dp, 0._dp, 3, 0, &
            [0._dp, 0._dp, 1000._dp])

       call source_flow_test("produce all", -5._dp, 0._dp, 0, 0, &
            [-3.4948610582_dp, -1.5051389418_dp, -431766.653977922_dp])
       call source_flow_test("produce 1", -5._dp, 0._dp, 0, 1, &
            [-5._dp, 0._dp, -431766.653977922_dp])
       call source_flow_test("produce heat", -5000._dp, 0._dp, 0, 3, &
            [0._dp, 0._dp, -5000._dp])

       call source_flow_test("no flow 1", 0._dp, 100.e3_dp, 1, 0, &
            [0._dp, 0._dp, 0._dp])
       call source_flow_test("no flow all", 0._dp, 100.e3_dp, 0, 0, &
            [0._dp, 0._dp, 0._dp])

       call fluid%destroy()
       deallocate(fluid_data)

    end if

  contains

    subroutine source_flow_test(tag, rate, enthalpy, injection_component, &
         production_component, flow)
      !! Runs asserts for single flow update_source() test.

      character(*), intent(in) :: tag
      PetscReal, intent(in) :: rate, enthalpy
      PetscInt, intent(in) :: injection_component, production_component
      PetscReal, intent(in) :: flow(:)

      call source%init(0, 0, num_primary, rate, enthalpy, &
           injection_component, production_component)
      call source%update_flow(fluid, isothermal)
      call assert_equals(flow, source%flow, &
           num_primary, tol, "Source update_flow() " // trim(tag))
      call source%destroy()

    end subroutine source_flow_test

  end subroutine test_source_update_flow

!------------------------------------------------------------------------

end module source_test
