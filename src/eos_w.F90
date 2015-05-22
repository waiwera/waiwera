module eos_w_module
  !! Isothermal pure water equation of state.

  use kinds_module
  use eos_module
  use thermodynamics_module

  implicit none
  private

#include <petsc-finclude/petsc.h90>

  type, public, extends(eos_type) :: eos_w_type
     !! Isothermal pure water equation of state type.
     private
     PetscReal, public :: temperature = 20._dp !! Constant temperature
   contains
     private
     procedure :: transition => eos_w_transition
     procedure, public :: init => eos_w_init
     procedure, public :: check_primary => eos_w_check_primary
     procedure, public :: fluid_properties => eos_w_fluid_properties
  end type eos_w_type

contains

!------------------------------------------------------------------------

  subroutine eos_w_init(self, thermo)
    !! Initialise isothermal pure water EOS.

    class(eos_w_type), intent(in out) :: self
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object

    self%name = "W"
    self%description = "Isothermal pure water"
    self%primary_variable_names = ["Pressure"]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 1
    self%num_components = 1

    self%thermo => thermo

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
    type(fluid_type), intent(out) :: fluid !! Fluid object
    ! Locals:
    PetscReal :: properties(2)
    PetscInt :: err

    fluid%pressure = primary(1)
    fluid%temperature = self%temperature
    fluid%region = 1

    call self%thermo%water%properties([fluid%pressure, fluid%temperature],\
    properties, err)

    fluid%phase(1)%density = properties(1)
    fluid%phase(1)%internal_energy = properties(2)

    call self%thermo%water%viscosity(fluid%temperature, fluid%pressure, \
    fluid%phase(1)%density, fluid%phase(1)%viscosity)

  end subroutine eos_w_fluid_properties

!------------------------------------------------------------------------

end module eos_w_module
