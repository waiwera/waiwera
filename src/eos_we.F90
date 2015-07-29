module eos_we_module
  !! Pure water and energy equation of state.

  use kinds_module
  use eos_module
  use fson
  use thermodynamics_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  type, public, extends(eos_type) :: eos_we_type
     !! Pure water and energy equation of state type.
     private
   contains
     private
     procedure :: transition => eos_we_transition
     procedure, public :: init => eos_we_init
     procedure, public :: check_primary => eos_we_check_primary
     procedure, public :: fluid_properties => eos_we_fluid_properties
  end type eos_we_type

contains

!------------------------------------------------------------------------

  subroutine eos_we_init(self, json, thermo)
    !! Initialise pure water and energy EOS.

    use fson_mpi_module, only: fson_get_mpi

    class(eos_we_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    ! Locals:
    PetscReal, parameter :: default_temperature = 20._dp ! deg C

    self%name = "we"
    self%description = "Pure water and energy"
    self%primary_variable_names = ["Pressure", "Temperature"]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 2
    self%num_components = 1

    self%thermo => thermo

  end subroutine eos_we_init

!------------------------------------------------------------------------
  
  subroutine eos_we_transition(self, region1, region2, primary)
    !! Perform transitions between thermodynamic regions.

    class(eos_we_type), intent(in) :: self
    PetscInt, intent(in) :: region1, region2
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)

    continue ! TODO

  end subroutine eos_we_transition

!------------------------------------------------------------------------

  subroutine eos_we_check_primary(self, region, primary)
    !! Check primary variables for current region and make
    !! transition if needed.

    class(eos_we_type), intent(in) :: self
    PetscInt, intent(in out) :: region
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)

    continue ! TODO

  end subroutine eos_we_check_primary

!------------------------------------------------------------------------

  subroutine eos_we_fluid_properties(self, region, primary, fluid)
    !! Calculate fluid properties from region and primary variables.

    use fluid_module, only: fluid_type
    class(eos_we_type), intent(in out) :: self
    PetscInt, intent(in) :: region !! Thermodynamic region index
    PetscReal, intent(in), target :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    ! Locals:
    PetscReal :: properties(2)
    PetscInt :: err

    fluid%pressure = primary(1)
    fluid%temperature = self%temperature
    fluid%region = 1

    call self%thermo%water%properties([fluid%pressure, fluid%temperature], &
         properties, err)

    fluid%phase(1)%saturation = 1._dp
    fluid%phase(1)%relative_permeability = 1._dp
    fluid%phase(1)%density = properties(1)
    fluid%phase(1)%internal_energy = properties(2)
    fluid%phase(1)%specific_enthalpy = fluid%phase(1)%internal_energy + &
         fluid%pressure / fluid%phase(1)%density
    fluid%phase(1)%mass_fraction(1) = 1._dp

    call self%thermo%water%viscosity(fluid%temperature, fluid%pressure, &
         fluid%phase(1)%density, fluid%phase(1)%viscosity)

  end subroutine eos_we_fluid_properties

!------------------------------------------------------------------------

end module eos_we_module
