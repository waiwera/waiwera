module eos_wsce_module
  !! Equation of state for non-isothermal water, salt and CO2 NCG.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use eos_wsge_module
  use ncg_co2_thermodynamics_module

  implicit none
  private

  type, public, extends(eos_wsge_type) :: eos_wsce_type
     !! Equation of state object for non-isothermal water, salt and CO2 NCG.
   contains
     private
     procedure, public :: init => eos_wsce_init
     procedure, public :: destroy => eos_wsce_destroy
  end type eos_wsce_type

contains

!------------------------------------------------------------------------

  subroutine eos_wsce_init(self, json, thermo, logfile)
    !! Initialise eos_wsce object.

    use fson
    use fson_mpi_module, only: fson_get_mpi, fson_has_mpi
    use logfile_module
    use thermodynamics_module

    class(eos_wsce_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile

    call self%eos_wsge_type%init(json, thermo, logfile)

    self%name = "wsce"
    self%description = "Water, salt, CO2 NCG and energy"
    self%primary_variable_names(4) = "CO2 partial pressure"
    self%component_names(3) = "CO2"
    self%required_output_fluid_fields = [ &
         "pressure                 ", &
         "temperature              ", &
         "region                   ", &
         "vapour_saturation        ", &
         "liquid_salt_mass_fraction", &
         "solid_saturation         ", &
         "CO2_partial_pressure     "]

    self%default_output_fluid_fields = self%required_output_fluid_fields

    allocate(ncg_co2_thermodynamics_type :: self%gas)
    call self%gas%init()

  end subroutine eos_wsce_init

!------------------------------------------------------------------------

  subroutine eos_wsce_destroy(self)
    !! Destroy eos_wsce object.

    class(eos_wsce_type), intent(in out) :: self

    deallocate(self%gas)
    call self%eos_wsge_type%destroy()

  end subroutine eos_wsce_destroy

!------------------------------------------------------------------------

end module eos_wsce_module
