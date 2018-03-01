module eos_wce_module
  !! Equation of state for non-isothermal water and CO2 NCG.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use eos_wge_module
  use ncg_co2_thermodynamics_module

  implicit none
  private

  type, public, extends(eos_wge_type) :: eos_wce_type
     !! Equation of state object for non-isothermal water and CO2 NCG.
   contains
     private
     procedure, public :: init => eos_wce_init
     procedure, public :: destroy => eos_wce_destroy
  end type eos_wce_type

contains

!------------------------------------------------------------------------

  subroutine eos_wce_init(self, json, thermo, logfile)
    !! Initialise eos_wce object.

    use fson
    use fson_mpi_module, only: fson_get_mpi
    use logfile_module
    use thermodynamics_module

    class(eos_wce_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile

    call self%eos_wge_type%init(json, thermo, logfile)

    self%name = "wce"
    self%description = "Water, CO2 NCG and energy"
    self%primary_variable_names(3) = "CO2 partial pressure"
    self%component_names(2) = "CO2"
    self%required_output_fluid_fields = [ &
         "pressure            ", "temperature         ", &
         "region              ", "CO2_partial_pressure", &
         "vapour_saturation   "]
    self%default_output_fluid_fields = [ &
         "pressure              ", "temperature           ", &
         "region                ", "CO2_partial_pressure  ", &
         "vapour_saturation     "]

    allocate(ncg_co2_thermodynamics_type :: self%gas)
    call self%gas%init()

  end subroutine eos_wce_init

!------------------------------------------------------------------------

  subroutine eos_wce_destroy(self)
    !! Destroy eos_wce object.

    class(eos_wce_type), intent(in out) :: self

    deallocate(self%gas)
    call self%eos_wge_type%destroy()

  end subroutine eos_wce_destroy

!------------------------------------------------------------------------

end module eos_wce_module
