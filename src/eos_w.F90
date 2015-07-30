module eos_w_module
  !! Isothermal pure water equation of state.

  use kinds_module
  use eos_module
  use fson
  use thermodynamics_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  type, public, extends(eos_type) :: eos_w_type
     !! Isothermal pure water equation of state type.
     private
     PetscReal, public :: temperature  !! Constant temperature
   contains
     private
     procedure :: transition => eos_w_transition
     procedure, public :: init => eos_w_init
     procedure, public :: check_primary => eos_w_check_primary
     procedure, public :: fluid_properties => eos_w_fluid_properties
  end type eos_w_type

contains

!------------------------------------------------------------------------

  subroutine eos_w_init(self, json, thermo)
    !! Initialise isothermal pure water EOS.

    use fson_mpi_module, only: fson_get_mpi

    class(eos_w_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    ! Locals:
    PetscReal, parameter :: default_temperature = 20._dp ! deg C

    self%name = "w"
    self%description = "Isothermal pure water"
    self%primary_variable_names = ["Pressure"]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 1
    self%num_components = 1
    self%isothermal = .true.

    self%thermo => thermo

    call fson_get_mpi(json, "eos.temperature", default_temperature, &
         self%temperature)

  end subroutine eos_w_init

!------------------------------------------------------------------------
  
  subroutine eos_w_transition(self, region1, region2, primary)
    !! Perform transitions between thermodynamic regions for isothermal
    !! pure water

    class(eos_w_type), intent(in) :: self
    PetscInt, intent(in) :: region1, region2
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)

    continue ! no transitions needed

  end subroutine eos_w_transition

!------------------------------------------------------------------------

  subroutine eos_w_check_primary(self, region, primary)
    !! Check primary variables for current region and make
    !! transition if needed for isothermal pure water

    class(eos_w_type), intent(in) :: self
    PetscInt, intent(in out) :: region
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)

    continue ! no checks needed

  end subroutine eos_w_check_primary

!------------------------------------------------------------------------

  subroutine eos_w_fluid_properties(self, region, primary, fluid)
    !! Calculate fluid properties from region and primary variables
    !! for isothermal pure water.

    use fluid_module, only: fluid_type
    class(eos_w_type), intent(in out) :: self
    PetscInt, intent(in) :: region !! Thermodynamic region index
    PetscReal, intent(in), target :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object

    fluid%pressure = primary(1)
    fluid%temperature = self%temperature
    fluid%region = 1

    call fluid%phase_properties(self%thermo)

  end subroutine eos_w_fluid_properties

!------------------------------------------------------------------------

end module eos_w_module
