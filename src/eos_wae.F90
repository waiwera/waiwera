module eos_wae_module
  !! Equation of state for non-isothermal water and air NCG.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use eos_wge_module
  use ncg_air_thermodynamics_module

  implicit none
  private

  type, public, extends(eos_wge_type) :: eos_wae_type
     !! Equation of state object for non-isothermal water and air NCG.
   contains
     private
     procedure, public :: init => eos_wae_init
     procedure, public :: destroy => eos_wae_destroy
  end type eos_wae_type

contains

!------------------------------------------------------------------------

  subroutine eos_wae_init(self, json, thermo, logfile)
    !! Initialise eos_wae object.

    use fson
    use fson_mpi_module, only: fson_get_mpi, fson_has_mpi
    use logfile_module
    use thermodynamics_module

    class(eos_wae_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscReal :: air_partial_pressure_scale

    call self%eos_wge_type%init(json, thermo, logfile)

    self%name = "wae"
    self%description = "Water, air NCG and energy"
    self%primary_variable_names(3) = "air partial pressure"
    self%component_names(2) = "air"
    self%required_output_fluid_fields = [ &
         "pressure            ", "temperature         ", &
         "region              ", "air_partial_pressure", &
         "vapour_saturation   "]
    self%default_output_fluid_fields = [ &
         "pressure              ", "temperature           ", &
         "region                ", "air_partial_pressure  ", &
         "vapour_saturation     "]

    if (fson_has_mpi(json, "eos.primary.scale.air_partial_pressure")) then
       call fson_get_mpi(json, "eos.primary.scale.air_partial_pressure", &
         val = air_partial_pressure_scale)
       self%primary_scale(3, :) = air_partial_pressure_scale
    end if

    allocate(ncg_air_thermodynamics_type :: self%gas)
    call self%gas%init()

  end subroutine eos_wae_init

!------------------------------------------------------------------------

  subroutine eos_wae_destroy(self)
    !! Destroy eos_wae object.

    class(eos_wae_type), intent(in out) :: self

    deallocate(self%gas)
    call self%eos_wge_type%destroy()

  end subroutine eos_wae_destroy

!------------------------------------------------------------------------

end module eos_wae_module
