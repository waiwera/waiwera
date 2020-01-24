!   Copyright 2016 University of Auckland.

!   This file is part of Waiwera.

!   Waiwera is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.

!   Waiwera is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.

!   You should have received a copy of the GNU Lesser General Public License
!   along with Waiwera.  If not, see <http://www.gnu.org/licenses/>.

module eos_w_module
  !! Equation of state for isothermal pure water.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use eos_module

  implicit none
  private

  type, public, extends(eos_type) :: eos_w_type
     !! Isothermal pure water equation of state type.
     private
     PetscReal, public :: temperature  !! Constant temperature
   contains
     private
     procedure, public :: init => eos_w_init
     procedure, public :: transition => eos_w_transition
     procedure, public :: bulk_properties => eos_w_bulk_properties
     procedure, public :: phase_saturations => eos_w_phase_saturations
     procedure, public :: phase_properties => eos_w_phase_properties
     procedure, public :: primary_variables => eos_w_primary_variables
     procedure, public :: check_primary_variables => eos_w_check_primary_variables
     procedure, public :: conductivity => eos_w_conductivity
  end type eos_w_type

contains

!------------------------------------------------------------------------

  subroutine eos_w_init(self, json, thermo, logfile)
    !! Initialise isothermal pure water EOS.

    use fson
    use fson_mpi_module, only: fson_get_mpi
    use logfile_module
    use thermodynamics_module

    class(eos_w_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscReal :: pressure_scale
    PetscReal, parameter :: default_pressure = 1.0e5_dp
    PetscReal, parameter :: default_temperature = 20._dp ! deg C
    PetscReal, parameter :: default_pressure_scale = 1.e6_dp !! Scale factor for non-dimensionalising pressure

    self%name = "w"
    self%description = "Isothermal pure water"
    self%primary_variable_names = ["pressure"]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 1
    self%phase_names = ["liquid"]
    self%num_components = 1
    self%component_names = ["water"]
    self%isothermal = PETSC_TRUE

    self%default_primary = [default_pressure]
    self%default_region = 1
    self%required_output_fluid_fields = ["pressure", "region  "]
    self%default_output_fluid_fields = ["pressure", "region  "]

    call fson_get_mpi(json, "eos.primary.scale.pressure", default_pressure_scale, &
         pressure_scale, logfile)

    allocate(self%primary_scale(1, 2))
    self%primary_scale = reshape([ &
          pressure_scale, &
          pressure_scale], [1, 2])

    self%thermo => thermo

    call fson_get_mpi(json, "eos.temperature", default_temperature, &
         self%temperature, logfile)

  end subroutine eos_w_init

!------------------------------------------------------------------------
  
  subroutine eos_w_transition(self, old_primary, primary, &
       old_fluid, fluid, transition, err)
    !! For eos_w, check primary variables for a cell and make
    !! thermodynamic region transitions if needed

    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in out) :: self
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in) :: old_fluid
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err

    ! no transitions needed
    err = 0
    transition = PETSC_FALSE

  end subroutine eos_w_transition

!------------------------------------------------------------------------

  subroutine eos_w_bulk_properties(self, primary, fluid, err)
    !! Calculate fluid bulk properties from region and primary variables
    !! for isothermal pure water.

    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err

    err = 0
    fluid%pressure = primary(1)
    fluid%temperature = self%temperature

    fluid%phase(1)%saturation = 1._dp
    call self%phase_composition(fluid, err)

    fluid%partial_pressure(1) = fluid%pressure

  end subroutine eos_w_bulk_properties

!------------------------------------------------------------------------

  subroutine eos_w_phase_saturations(self, primary, fluid)
    !! Assigns fluid phase saturations from fluid region and primary variables.

    use fluid_module, only: fluid_type
    class(eos_w_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    ! Locals:
    PetscInt :: region

    region = nint(fluid%region)
    fluid%phase(region)%saturation = 1._dp

  end subroutine eos_w_phase_saturations

!------------------------------------------------------------------------

  subroutine eos_w_phase_properties(self, primary, rock, fluid, err)
    !! Calculate fluid phase properties from region and primary variables
    !! for isothermal pure water.
    !! Bulk properties need to be calculated before calling this routine.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_w_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: p
    PetscReal :: properties(2)

    err = 0
    p = nint(fluid%region)

    associate(region => self%thermo%region(p)%ptr)

      call region%properties([fluid%pressure, fluid%temperature], &
           properties, err)

      if (err == 0) then
         associate(phase => fluid%phase(p))

         phase%density = properties(1)
         phase%internal_energy = properties(2)
         phase%specific_enthalpy =  phase%internal_energy + &
              fluid%pressure / phase%density

         phase%relative_permeability = 1._dp
         phase%capillary_pressure = 0._dp
         phase%mass_fraction(1) = 1._dp

         call region%viscosity(fluid%temperature, fluid%pressure, &
              phase%density, phase%viscosity)

       end associate
    end if

  end associate

end subroutine eos_w_phase_properties

!------------------------------------------------------------------------

  subroutine eos_w_primary_variables(self, fluid, primary)
    !! Determine primary variables from fluid properties.

    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(out) :: primary(self%num_primary_variables)

    primary(1) = fluid%pressure

  end subroutine eos_w_primary_variables

!------------------------------------------------------------------------

  subroutine eos_w_check_primary_variables(self, fluid, &
       primary, changed, err)
    !! Check if primary variables are in acceptable bounds, and return
    !! error code accordingly.

    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    PetscBool, intent(out) :: changed
    PetscErrorCode, intent(out) :: err

    changed = PETSC_FALSE

    associate (p => primary(1))
      if ((p < 0._dp) .or. (p > 100.e6_dp)) then
         err = 1
      else
         err = 0
      end if
    end associate

  end subroutine eos_w_check_primary_variables

!------------------------------------------------------------------------

  PetscReal function eos_w_conductivity(self, rock, fluid) result(cond)
    !! Returns effective rock heat conductivity for given fluid properties.

    use rock_module, only: rock_type
    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in) :: self
    type(rock_type), intent(in) :: rock !! Rock object
    type(fluid_type), intent(in) :: fluid !! Fluid object

    cond = rock%wet_conductivity

  end function eos_w_conductivity

!------------------------------------------------------------------------

end module eos_w_module
