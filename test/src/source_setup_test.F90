module source_setup_test

  ! Tests for source setup module

  use kinds_module
  use mpi_module
  use fruit
  use eos_module
  use source_setup_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  type, extends(eos_type) :: eos_test_type
     !! EOS used only for testing source setup.
   contains
     private
     procedure, public :: init => eos_test_init
     procedure, public :: transition => eos_test_transition
     procedure, public :: bulk_properties => eos_test_bulk_properties
     procedure, public :: phase_properties => eos_test_phase_properties
     procedure, public :: primary_variables => eos_test_primary_variables
     procedure, public :: check_primary_variables => eos_test_check_primary_variables
  end type eos_test_type

public :: test_setup_sources

contains

!------------------------------------------------------------------------

  subroutine eos_test_init(self, json, thermo, logfile)

    use fson
    use fson_mpi_module, only: fson_get_mpi
    use logfile_module
    use thermodynamics_module

    class(eos_test_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    class(thermodynamics_type), intent(in), target :: thermo
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscReal, parameter :: default_pressure = 1.0e5_dp
    PetscReal, parameter :: default_temperature = 20._dp
    PetscReal, parameter :: default_gas_partial_pressure = 0._dp

    self%name = "test"
    self%description = "Test EOS"
    self%primary_variable_names = [&
         "pressure                     ", &
         "temperature/vapour_saturation", &
         "gas partial pressure         "]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 2
    self%phase_names = ["liquid", "vapour"]
    self%num_components = 2
    self%component_names = ["water", "gas  "]

    self%default_primary = [default_pressure, default_temperature, &
         default_gas_partial_pressure]
    self%default_region = 1

    self%thermo => thermo

  end subroutine eos_test_init

!------------------------------------------------------------------------

  subroutine eos_test_transition(self, primary, old_fluid, fluid, &
       transition, err)

    use fluid_module, only: fluid_type

    class(eos_test_type), intent(in out) :: self
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in) :: old_fluid
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err

    transition = PETSC_FALSE
    err = 0

  end subroutine eos_test_transition

!------------------------------------------------------------------------

  subroutine eos_test_bulk_properties(self, primary, fluid, err)

    use fluid_module, only: fluid_type

    class(eos_test_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscErrorCode, intent(out) :: err

    err = 0

  end subroutine eos_test_bulk_properties

!------------------------------------------------------------------------

  subroutine eos_test_phase_properties(self, primary, rock, fluid, err)

    use rock_module, only: rock_type
    use fluid_module, only: fluid_type

    class(eos_test_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables)
    type(rock_type), intent(in out) :: rock
    type(fluid_type), intent(in out) :: fluid
    PetscErrorCode, intent(out) :: err

    err = 0

  end subroutine eos_test_phase_properties

!------------------------------------------------------------------------

  subroutine eos_test_primary_variables(self, fluid, primary)

    use fluid_module, only: fluid_type

    class(eos_test_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(out) :: primary(self%num_primary_variables)

    primary = 0._dp

  end subroutine eos_test_primary_variables

!------------------------------------------------------------------------

  PetscErrorCode function eos_test_check_primary_variables(self, fluid, &
       primary) result(err)

    use fluid_module, only: fluid_type

    class(eos_test_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(in) :: primary(self%num_primary_variables)

    err = 0

  end function eos_test_check_primary_variables

!------------------------------------------------------------------------

  subroutine test_setup_sources

    ! setup_sources() test

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use source_module
    use IAPWS_module

    character(16), parameter :: path = "data/source/"
    type(IAPWS_type) :: thermo
    type(eos_test_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(list_type) :: sources, source_controls
    PetscErrorCode :: ierr

    json => fson_parse_mpi(trim(path) // "test_source.json")

    call thermo%init()
    call eos%init(json, thermo)
    call mesh%init(json)
    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(eos%primary_variable_names)

    call setup_sources(json, mesh%dm, eos, sources, source_controls)

    call sources%traverse(source_test_iterator)
    call sources%destroy(source_list_node_data_destroy)

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
