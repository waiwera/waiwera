module eos_sae_module
  !! Equation of state for supercritical water and air NCG.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use eos_sge_module
  use ncg_air_thermodynamics_module

  implicit none
  private

  type, public, extends(eos_sge_type) :: eos_sae_type
     !! Equation of state object for supercritical water and air NCG.
   contains
     private
     procedure, public :: init => eos_sae_init
  end type eos_sae_type

contains

!------------------------------------------------------------------------

  subroutine eos_sae_init(self, json, thermo, logfile)
    !! Initialise eos_sae object.

    use fson
    use fson_mpi_module, only: fson_get_mpi, fson_has_mpi
    use logfile_module
    use thermodynamics_module

    class(eos_sae_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscReal :: air_partial_pressure_scale

    call self%eos_sge_type%init(json, thermo, logfile)

    self%name = "sae"
    self%description = "Supercritical water, air NCG and energy"
    self%primary_variable_names(3) = "air partial pressure"
    self%component_names(2) = "air"
    self%required_output_fluid_fields = [ &
         "pressure             ", "temperature          ", &
         "region               ", "vapour_saturation    ", &
         "liquid_density       ", "vapour_density       ", &
         "supercritical_density", "air_partial_pressure "]
    self%default_output_fluid_fields = [ &
         "pressure             ", "temperature          ", &
         "region               ", "vapour_saturation    ", &
         "liquid_density       ", "vapour_density       ", &
         "supercritical_density", "liquidlike_fraction  ", &
         "air_partial_pressure "]

    if (fson_has_mpi(json, "eos.primary.scale.air_partial_pressure")) then
       call fson_get_mpi(json, "eos.primary.scale.air_partial_pressure", &
         val = air_partial_pressure_scale)
       self%primary_scale(3, :) = air_partial_pressure_scale
    end if

    allocate(ncg_air_thermodynamics_type :: self%gas)
    call self%gas%init()

  end subroutine eos_sae_init

!------------------------------------------------------------------------

end module eos_sae_module
