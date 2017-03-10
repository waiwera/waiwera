module eos_test
  !! Tests for eos module.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fruit
  use fson
  use fson_mpi_module
  use eos_module

  implicit none
  private

  type, extends(eos_type), public :: eos_test_type
     !! Dummy EOS used only for testing.
   contains
     private
     procedure, public :: init => eos_test_init
     procedure, public :: transition => eos_test_transition
     procedure, public :: bulk_properties => eos_test_bulk_properties
     procedure, public :: phase_properties => eos_test_phase_properties
     procedure, public :: primary_variables => eos_test_primary_variables
     procedure, public :: check_primary_variables => eos_test_check_primary_variables
  end type eos_test_type

  public :: test_eos_component_index

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

  subroutine eos_test_transition(self, old_primary, primary, &
       old_fluid, fluid, transition, err)

    use fluid_module, only: fluid_type

    class(eos_test_type), intent(in out) :: self
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
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

  subroutine test_eos_component_index
    ! eos component_index() test

    use IAPWS_module

    type(IAPWS_type) :: thermo
    type(eos_test_type) :: eos
    type(fson_value), pointer :: json
    character(max_component_name_length) :: name

    json => fson_parse_mpi(str = "{}")
    call thermo%init()
    call eos%init(json, thermo)

    name = "water"
    call assert_equals(1, eos%component_index(name), name)

    name = "Water"
    call assert_equals(1, eos%component_index(name), name)

    name = "gas"
    call assert_equals(2, eos%component_index(name), name)

    name = "energy"
    call assert_equals(eos%num_primary_variables, &
         eos%component_index(name), name)

    name = "fred"
    call assert_equals(-1, eos%component_index(name), name)

    call eos%destroy()
    call thermo%destroy()

  end subroutine test_eos_component_index

!------------------------------------------------------------------------

end module eos_test
