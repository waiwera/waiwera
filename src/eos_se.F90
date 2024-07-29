!   Copyright 2024 University of Auckland.

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

module eos_se_module
  !! Equation of state for non-isothermal pure water, sub- or super-critical.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use eos_module
  use eos_we_module
  use root_finder_module
  use thermodynamics_module

  implicit none
  private

  type, public, extends(eos_we_type) :: eos_se_type
     !! Pure supercritical water and energy equation of state type.
     private
   contains
     private
     procedure, public :: init => eos_se_init
     procedure, public :: destroy => eos_se_destroy
     procedure, public :: transition => eos_se_transition
     procedure, public :: bulk_properties => eos_se_bulk_properties
     procedure, public :: phase_properties => eos_se_phase_properties
     procedure, public :: primary_variables => eos_se_primary_variables
     procedure, public :: phase_saturations => eos_se_phase_saturations
     procedure, public :: check_primary_variables => eos_se_check_primary_variables
  end type eos_se_type

contains

!------------------------------------------------------------------------

  subroutine eos_se_init(self, json, thermo, logfile)
    !! Initialise pure supercritical water and energy EOS.

    use fson
    use fson_mpi_module, only: fson_get_mpi
    use logfile_module
    use thermodynamics_module
    use IAPWS_module, only: critical

    class(eos_se_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    procedure(root_finder_function), pointer :: f
    class(*), pointer :: pinterp
    PetscReal, allocatable :: data(:, :)
    PetscReal :: pressure_scale, temperature_scale, density_scale
    PetscReal, parameter :: default_pressure = 1.0e5_dp
    PetscReal, parameter :: default_temperature = 20._dp ! deg C
    PetscReal, parameter :: default_pressure_scale = 1.e6_dp !! Default scale factor for non-dimensionalising pressure
    PetscReal, parameter :: default_temperature_scale = 1.e2_dp !! Default scale factor for non-dimensionalising temperature
    PetscReal, parameter :: default_density_scale = critical%density !! Default scale factor for non-dimensionalising density

    self%name = "se"
    self%description = "Pure supercritical water and energy"
    self%primary_variable_names = ["pressure/density             ", &
         "temperature/vapour_saturation"]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 3
    self%num_mobile_phases = 3
    self%phase_names = ["liquid       ", "vapour       ", "supercritical"]
    self%num_components = 1
    self%component_names = ["water"]

    self%default_primary = [default_pressure, default_temperature]
    self%default_region = 1
    self%default_tracer_phase = "liquid"
    self%required_output_fluid_fields = [ &
         "pressure             ", "temperature          ", &
         "region               ", "vapour_saturation    ", &
         "supercritical_density"]
    self%default_output_fluid_fields = [ &
         "pressure             ", "temperature          ", &
         "region               ", "vapour_saturation    ", &
         "supercritical_density"]

    call fson_get_mpi(json, "eos.primary.scale.pressure", default_pressure_scale, &
         pressure_scale, logfile)
    call fson_get_mpi(json, "eos.primary.scale.temperature", default_temperature_scale, &
         temperature_scale, logfile)
    call fson_get_mpi(json, "eos.primary.scale.density", default_density_scale, &
         density_scale, logfile)
    allocate(self%primary_scale(2, 4))
    self%primary_scale = reshape([ &
          pressure_scale, temperature_scale, &
          pressure_scale, temperature_scale, &
          density_scale, temperature_scale, &
          pressure_scale, 1._dp], [2, 4])

    self%thermo => thermo

    ! Set up saturation line finder:
    allocate(primary_variable_interpolator_type :: &
         self%primary_variable_interpolator)
    allocate(data(2, 1 + self%num_primary_variables))
    data = 0._dp
    data(:, 1) = [0._dp, 1._dp]
    call self%primary_variable_interpolator%init(data)
    deallocate(data)
    self%primary_variable_interpolator%thermo => self%thermo
    f => eos_we_saturation_difference
    pinterp => self%primary_variable_interpolator
    call self%saturation_line_finder%init(f, context = pinterp)

  end subroutine eos_se_init

!------------------------------------------------------------------------

  subroutine eos_se_destroy(self)
    !! Destroys pure supercritical water and energy EOS.

    class(eos_se_type), intent(in out) :: self

    deallocate(self%primary_variable_names)
    deallocate(self%phase_names, self%component_names)
    deallocate(self%default_primary)
    deallocate(self%primary_scale)
    self%thermo => null()

    call self%saturation_line_finder%destroy()
    call self%primary_variable_interpolator%destroy()
    deallocate(self%primary_variable_interpolator)

  end subroutine eos_se_destroy

!------------------------------------------------------------------------

  subroutine eos_se_transition(self, old_primary, primary, &
       old_fluid, fluid, transition, err)
    !! For eos_se, check primary variables for a cell and make
    !! thermodynamic region transitions if needed.

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in) :: old_fluid
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: old_region
    PetscReal :: saturation_pressure

    err = 0
    transition = PETSC_FALSE
    old_region = nint(old_fluid%region)

    if (old_region == 4) then  ! Two-phase
       associate (vapour_saturation => primary(2))

         if (vapour_saturation < 0._dp) then
            call self%transition_to_single_phase(old_primary, old_fluid, &
                 1, primary, fluid, transition, err)
         else if (vapour_saturation > 1._dp) then
            call self%transition_to_single_phase(old_primary, old_fluid, &
                 2, primary, fluid, transition, err)
         end if

     end associate
    else  ! Single-phase
       associate (pressure => primary(1), temperature => primary(2))

         call self%thermo%saturation%pressure(temperature, &
              saturation_pressure, err)

         if (err == 0) then
            if (((old_region == 1) .and. (pressure < saturation_pressure)) .or. &
                 ((old_region == 2) .and. (pressure > saturation_pressure))) then
               call self%transition_to_two_phase(saturation_pressure, &
                    old_primary, old_fluid, primary, fluid, transition, err)
            end if
         end if

       end associate
    end if

  end subroutine eos_se_transition

!------------------------------------------------------------------------

  subroutine eos_se_bulk_properties(self, primary, fluid, err)
    !! Calculate fluid bulk properties from region and primary variables
    !! for non-isothermal pure supercritical water.

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscInt :: region

    err = 0
    fluid%pressure = primary(1)
    region = nint(fluid%region)

    if (region == 4) then
       ! Two-phase
       call self%thermo%saturation%temperature(fluid%pressure, &
            fluid%temperature, err)
    else
       ! Single-phase
       fluid%temperature = primary(2)
    end if

    if (err == 0) then
       fluid%permeability_factor = 1._dp
       call self%phase_composition(fluid, err)
       if (err == 0) then
          call self%phase_saturations(primary, fluid)
          fluid%partial_pressure(1) = fluid%pressure
       end if
    end if

  end subroutine eos_se_bulk_properties

!------------------------------------------------------------------------

  subroutine eos_se_phase_saturations(self, primary, fluid)
    !! Assigns fluid phase saturations from fluid region and primary variables.

    use fluid_module, only: fluid_type
    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    ! Locals:
    PetscInt :: region
    PetscReal :: saturation_pressure
    PetscErrorCode :: err

    region = nint(fluid%region)

    select case (region)
    case (1)
       call set_saturations(1._dp, 0._dp, 0._dp)
    case (2)
       call set_saturations(0._dp, 1._dp, 0._dp)
    case (3)
       if (fluid%temperature <= self%thermo%critical%temperature) then
         call self%thermo%saturation%pressure(fluid%temperature, &
              saturation_pressure, err)
         if (fluid%pressure > saturation_pressure) then
            call set_saturations(1._dp, 0._dp, 0._dp)
         else
            call set_saturations(0._dp, 1._dp, 0._dp)
         end if
       else
          if (fluid%pressure <= self%thermo%critical%pressure) then
            call set_saturations(0._dp, 1._dp, 0._dp)
          else
            call set_saturations(0._dp, 0._dp, 1._dp)
          end if
       end if
    case (4)
       call set_saturations(1._dp - primary(2), primary(2), 0._dp)
    end select

  contains

    subroutine set_saturations(s1, s2, s3)
      PetscReal, intent(in) :: s1, s2, s3

      fluid%phase(1)%saturation = s1
      fluid%phase(2)%saturation = s2
      fluid%phase(3)%saturation = s3

    end subroutine set_saturations

  end subroutine eos_se_phase_saturations

!------------------------------------------------------------------------

  subroutine eos_se_phase_properties(self, primary, rock, fluid, err)
    !! Calculate fluid phase properties from updated fluid region and primary variables.
    !! Bulk properties need to be calculated before calling this routine.

    use rock_module, only: rock_type
    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscInt :: p, phases
    PetscReal :: properties(2), sl, relative_permeability(2), capillary_pressure(2)

    err = 0
    phases = nint(fluid%phase_composition)

    sl = fluid%phase(1)%saturation
    relative_permeability = rock%relative_permeability%values(sl)
    capillary_pressure = [rock%capillary_pressure%value(sl, fluid%temperature), &
         0._dp]

    do p = 1, self%num_phases
       associate(phase => fluid%phase(p), &
            region => self%thermo%region(p)%ptr)

         if (btest(phases, p - 1)) then

            call region%properties([fluid%pressure, fluid%temperature], &
                 properties, err)

            if (err == 0) then

               phase%density = properties(1)
               phase%internal_energy = properties(2)
               phase%specific_enthalpy = phase%internal_energy + &
                    fluid%pressure / phase%density

               phase%mass_fraction(1) = 1._dp
               phase%relative_permeability = relative_permeability(p)
               phase%capillary_pressure =  capillary_pressure(p)

               call region%viscosity(fluid%temperature, fluid%pressure, &
                    phase%density, phase%viscosity)

            else
               exit
            end if

         else
            phase%density = 0._dp
            phase%internal_energy = 0._dp
            phase%specific_enthalpy = 0._dp
            phase%relative_permeability = 0._dp
            phase%capillary_pressure = 0._dp
            phase%viscosity = 0._dp
            phase%mass_fraction(1) = 0._dp
         end if

       end associate
    end do

  end subroutine eos_se_phase_properties

!------------------------------------------------------------------------

  subroutine eos_se_primary_variables(self, fluid, primary)
    !! Determine primary variables from fluid properties.

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(out) :: primary(self%num_primary_variables)
    ! Locals:
    PetscInt :: region

    primary(1) = fluid%pressure

    region = nint(fluid%region)
    if (region == 4) then
       primary(2) = fluid%phase(2)%saturation
    else
       primary(2) = fluid%temperature
    end if

  end subroutine eos_se_primary_variables

!------------------------------------------------------------------------

   subroutine eos_se_check_primary_variables(self, fluid, &
       primary, changed, err)
    !! Check if primary variables are in acceptable bounds, and return error
    !! code accordingly.

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    PetscBool, intent(out) :: changed
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: region

    changed = PETSC_FALSE
    err = 0

    associate (p => primary(1))
      if ((p < 0._dp) .or. (p > 100.e6_dp)) then
         err = 1
      else
         region = nint(fluid%region)
         if (region == 4) then
            associate (vapour_saturation => primary(2))
              if ((vapour_saturation < -1._dp) .or. &
                   (vapour_saturation > 2._dp)) then
                 err = 1
              end if
            end associate
         else
            associate (t => primary(2))
              if ((t < 0._dp) .or. (t > 800._dp)) then
                 err = 1
              end if
            end associate
         end if
      end if
    end associate

  end subroutine eos_se_check_primary_variables

!------------------------------------------------------------------------

end module eos_se_module
