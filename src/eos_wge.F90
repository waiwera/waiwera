module eos_wge_module
  !! Equation of state for non-isothermal water and non-condensible gas.

  use kinds_module
  use eos_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>
#include <petsc/finclude/petscsys.h>

  type, public, abstract, extends(eos_type) :: eos_wge_type
     !! Pure water, non-condensible gas and energy equation of state type.
     private
   contains
     private
     procedure, public :: init => eos_wge_init
  end type eos_wge_type

contains

!------------------------------------------------------------------------

  subroutine eos_wge_init(self, json, thermo, logfile)
    !! Initialise pure water, non-condensible gas and energy EOS.

    use fson
    use fson_mpi_module, only: fson_get_mpi
    use logfile_module
    use thermodynamics_module

    class(eos_wge_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscReal, parameter :: default_pressure = 1.0e5_dp
    PetscReal, parameter :: default_temperature = 20._dp ! deg C
    PetscReal, parameter :: default_gas_partial_pressure = 0._dp

    self%name = "wge"
    self%description = "Water, non-condensible gas and energy"
    self%primary_variable_names = [ &
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

  end subroutine eos_wge_init

!------------------------------------------------------------------------

end module eos_wge_module
