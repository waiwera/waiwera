module eos_wsae_module
  !! Equation of state for non-isothermal water, salt and air NCG.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use eos_wsge_module
  use ncg_air_thermodynamics_module

  implicit none
  private

  type, public, extends(eos_wsge_type) :: eos_wsae_type
     !! Equation of state object for non-isothermal water, salt and air NCG.
   contains
     private
     procedure, public :: init => eos_wsae_init
     procedure, public :: destroy => eos_wsae_destroy
  end type eos_wsae_type

contains

!------------------------------------------------------------------------

  subroutine eos_wsae_init(self, json, thermo, logfile)
    !! Initialise eos_wsae object.

    use fson
    use fson_mpi_module, only: fson_get_mpi, fson_has_mpi
    use logfile_module
    use thermodynamics_module

    class(eos_wsae_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscReal :: air_partial_pressure_scale

    call self%eos_wsge_type%init(json, thermo, logfile)

    self%name = "wsae"
    self%description = "Water, salt, air NCG and energy"
    self%primary_variable_names(4) = "air partial pressure"
    self%component_names(3) = "air"
    self%required_output_fluid_fields = [ &
         "pressure                 ", &
         "temperature              ", &
         "region                   ", &
         "vapour_saturation        ", &
         "liquid_salt_mass_fraction", &
         "solid_saturation         ", &
         "air_partial_pressure     "]

    self%default_output_fluid_fields = self%required_output_fluid_fields

    if (fson_has_mpi(json, "eos.primary.scale.air_partial_pressure")) then
       call fson_get_mpi(json, "eos.primary.scale.air_partial_pressure", &
         val = air_partial_pressure_scale)
       self%primary_scale(4, :) = air_partial_pressure_scale
    end if

    allocate(ncg_air_thermodynamics_type :: self%gas)
    call self%gas%init()

  end subroutine eos_wsae_init

!------------------------------------------------------------------------

  subroutine eos_wsae_destroy(self)
    !! Destroy eos_wsae object.

    class(eos_wsae_type), intent(in out) :: self

    deallocate(self%gas)
    call self%eos_wsge_type%destroy()

  end subroutine eos_wsae_destroy

!------------------------------------------------------------------------

end module eos_wsae_module
