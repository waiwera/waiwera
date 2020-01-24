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

module eos_we_module
  !! Equation of state for non-isothermal pure water.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use eos_module
  use root_finder_module
  use thermodynamics_module

  implicit none
  private

  type, public, extends(eos_type) :: eos_we_type
     !! Pure water and energy equation of state type.
     private
     type(root_finder_type) :: saturation_line_finder
     type(primary_variable_interpolator_type), pointer :: &
          primary_variable_interpolator
   contains
     private
     procedure, public :: init => eos_we_init
     procedure, public :: destroy => eos_we_destroy
     procedure, public :: transition => eos_we_transition
     procedure, public :: transition_to_single_phase => eos_we_transition_to_single_phase
     procedure, public :: transition_to_two_phase => eos_we_transition_to_two_phase
     procedure, public :: bulk_properties => eos_we_bulk_properties
     procedure, public :: phase_properties => eos_we_phase_properties
     procedure, public :: primary_variables => eos_we_primary_variables
     procedure, public :: phase_saturations => eos_we_phase_saturations
     procedure, public :: check_primary_variables => eos_we_check_primary_variables
  end type eos_we_type

contains

!------------------------------------------------------------------------

  subroutine eos_we_init(self, json, thermo, logfile)
    !! Initialise pure water and energy EOS.

    use fson
    use fson_mpi_module, only: fson_get_mpi
    use logfile_module
    use thermodynamics_module

    class(eos_we_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    procedure(root_finder_function), pointer :: f
    class(*), pointer :: pinterp
    PetscReal, allocatable :: data(:, :)
    PetscReal :: pressure_scale, temperature_scale
    PetscErrorCode :: err
    PetscReal, parameter :: default_pressure = 1.0e5_dp
    PetscReal, parameter :: default_temperature = 20._dp ! deg C
    PetscReal, parameter :: default_pressure_scale = 1.e6_dp !! Default scale factor for non-dimensionalising pressure
    PetscReal, parameter :: default_temperature_scale = 1.e2_dp !! Default scale factor for non-dimensionalising temperature

    self%name = "we"
    self%description = "Pure water and energy"
    self%primary_variable_names = ["pressure                     ", &
         "temperature/vapour_saturation"]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 2
    self%phase_names = ["liquid", "vapour"]
    self%num_components = 1
    self%component_names = ["water"]

    self%default_primary = [default_pressure, default_temperature]
    self%default_region = 1
    self%required_output_fluid_fields = [ &
         "pressure         ", "temperature      ", &
         "region           ", "vapour_saturation"]
    self%default_output_fluid_fields = [ &
         "pressure              ", "temperature           ", &
         "region                ", "vapour_saturation     "]

    call fson_get_mpi(json, "eos.primary.scale.pressure", default_pressure_scale, &
         pressure_scale, logfile)
    call fson_get_mpi(json, "eos.primary.scale.temperature", default_temperature_scale, &
         temperature_scale, logfile)
    allocate(self%primary_scale(2, 4))
    self%primary_scale = reshape([ &
          pressure_scale, temperature_scale, &
          pressure_scale, temperature_scale, &
          0._dp, 0._dp, &
          pressure_scale, 1._dp], [2, 4])

    self%thermo => thermo

    ! Set up saturation line finder:
    allocate(primary_variable_interpolator_type :: &
         self%primary_variable_interpolator)
    allocate(data(2, 1 + self%num_primary_variables))
    data = 0._dp
    data(:, 1) = [0._dp, 1._dp]
    call self%primary_variable_interpolator%init(data, err)
    deallocate(data)
    self%primary_variable_interpolator%thermo => self%thermo
    f => eos_we_saturation_difference
    pinterp => self%primary_variable_interpolator
    call self%saturation_line_finder%init(f, context = pinterp)

  end subroutine eos_we_init

!------------------------------------------------------------------------

  subroutine eos_we_destroy(self)
    !! Destroys pure water and energy EOS.

    class(eos_we_type), intent(in out) :: self

    deallocate(self%primary_variable_names)
    deallocate(self%phase_names, self%component_names)
    deallocate(self%default_primary)
    deallocate(self%primary_scale)
    self%thermo => null()

    call self%saturation_line_finder%destroy()
    call self%primary_variable_interpolator%destroy()
    deallocate(self%primary_variable_interpolator)

  end subroutine eos_we_destroy

!------------------------------------------------------------------------

  subroutine eos_we_transition_to_single_phase(self, old_primary, old_fluid, &
       new_region, primary, fluid, transition, err)
    !! For eos_we, make transition from two-phase to single-phase with
    !! specified region.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in out) :: self
    type(fluid_type), intent(in) :: old_fluid
    PetscInt, intent(in) :: new_region
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: old_saturation_pressure, pressure_factor
    PetscReal :: saturation_bound, xi
    PetscReal :: interpolated_primary(self%num_primary_variables)
    PetscReal, parameter :: small = 1.e-6_dp

    err = 0
    transition = PETSC_FALSE

    if (new_region == 1) then
       saturation_bound = 0._dp
       pressure_factor = 1._dp + small
    else
       saturation_bound = 1._dp
       pressure_factor = 1._dp - small
    end if

    self%primary_variable_interpolator%val(:, 1) = old_primary
    self%primary_variable_interpolator%val(:, 2) = primary
    call self%primary_variable_interpolator%find_component_at_index(&
         saturation_bound, 2, xi, err)

    associate (pressure => primary(1), temperature => primary(2), &
         interpolated_pressure => interpolated_primary(1))

      if (err == 0) then

         interpolated_primary = self%primary_variable_interpolator%interpolate(xi)
         pressure = pressure_factor * interpolated_pressure
         call self%thermo%saturation%temperature(interpolated_pressure, &
              temperature, err)
         if (err == 0) then
            fluid%region = dble(new_region)
            transition = PETSC_TRUE
         end if

      else

         call self%thermo%saturation%pressure(old_fluid%temperature, &
              old_saturation_pressure, err)
         if (err == 0) then
            pressure = pressure_factor * old_saturation_pressure
            temperature = old_fluid%temperature
            fluid%region = dble(new_region)
            transition = PETSC_TRUE
         end if

      end if

    end associate

  end subroutine eos_we_transition_to_single_phase

!------------------------------------------------------------------------

  subroutine eos_we_transition_to_two_phase(self, saturation_pressure, &
       old_primary, old_fluid, primary, fluid, transition, err)
    !! For eos_we, make transition from single-phase to two-phase.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in out) :: self
    PetscReal, intent(in) :: saturation_pressure
    type(fluid_type), intent(in) :: old_fluid
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: old_region
    PetscReal :: interpolated_primary(self%num_primary_variables)
    PetscReal :: xi
    PetscReal, parameter :: small = 1.e-6_dp

    err = 0
    associate (pressure => primary(1), vapour_saturation => primary(2), &
      interpolated_pressure => interpolated_primary(1))

      self%primary_variable_interpolator%val(:, 1) = old_primary
      self%primary_variable_interpolator%val(:, 2) = primary
      call self%saturation_line_finder%find()

      if (self%saturation_line_finder%err == 0) then
         xi = self%saturation_line_finder%root
         interpolated_primary = self%primary_variable_interpolator%interpolate(xi)
         pressure = interpolated_pressure
      else
         pressure = saturation_pressure
      end if

      old_region = nint(old_fluid%region)
      if (old_region == 1) then
         vapour_saturation = small
      else
         vapour_saturation = 1._dp - small
      end if

      fluid%region = dble(4)
      transition = PETSC_TRUE

    end associate

  end subroutine eos_we_transition_to_two_phase

!------------------------------------------------------------------------

  subroutine eos_we_transition(self, old_primary, primary, &
       old_fluid, fluid, transition, err)
    !! For eos_we, check primary variables for a cell and make
    !! thermodynamic region transitions if needed.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in out) :: self
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

  end subroutine eos_we_transition

!------------------------------------------------------------------------

  subroutine eos_we_bulk_properties(self, primary, fluid, err)
    !! Calculate fluid bulk properties from region and primary variables
    !! for non-isothermal pure water.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in out) :: self
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
       call self%phase_composition(fluid, err)
       if (err == 0) then
          call self%phase_saturations(primary, fluid)
          fluid%partial_pressure(1) = fluid%pressure
       end if
    end if

  end subroutine eos_we_bulk_properties

!------------------------------------------------------------------------

  subroutine eos_we_phase_saturations(self, primary, fluid)
    !! Assigns fluid phase saturations from fluid region and primary variables.

    use fluid_module, only: fluid_type
    class(eos_we_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    ! Locals:
    PetscInt :: region

    region = nint(fluid%region)

    select case (region)
    case (1)
       fluid%phase(1)%saturation = 1._dp
       fluid%phase(2)%saturation = 0._dp
    case (2)
       fluid%phase(1)%saturation = 0._dp
       fluid%phase(2)%saturation = 1._dp
    case (4)
       fluid%phase(1)%saturation = 1._dp - primary(2)
       fluid%phase(2)%saturation = primary(2)
    end select

  end subroutine eos_we_phase_saturations

!------------------------------------------------------------------------

  subroutine eos_we_phase_properties(self, primary, rock, fluid, err)
    !! Calculate fluid phase properties from updated fluid region and primary variables.
    !! Bulk properties need to be calculated before calling this routine.

    use rock_module, only: rock_type
    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in out) :: self
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

  end subroutine eos_we_phase_properties

!------------------------------------------------------------------------

  subroutine eos_we_primary_variables(self, fluid, primary)
    !! Determine primary variables from fluid properties.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in) :: self
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

  end subroutine eos_we_primary_variables

!------------------------------------------------------------------------

   subroutine eos_we_check_primary_variables(self, fluid, &
       primary, changed, err)
    !! Check if primary variables are in acceptable bounds, and return error
    !! code accordingly.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in) :: self
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

  end subroutine eos_we_check_primary_variables

!------------------------------------------------------------------------

  PetscReal function eos_we_saturation_difference(x, context) result(dp)
    !! Returns difference between saturation pressure and pressure at
    !! normalised point 0 <= x <= 1 along line between start and end
    !! primary variables.

    PetscReal, intent(in) :: x
    class(*), pointer, intent(in out) :: context
    ! Locals:
    PetscReal, allocatable :: var(:)
    PetscReal :: Ps
    PetscInt :: err

    select type (context)
    type is (primary_variable_interpolator_type)
       allocate(var(context%dim))
       var = context%interpolate_at_index(x)
       associate(P => var(1), T => var(2))
         call context%thermo%saturation%pressure(T, Ps, err)
         dp = P - Ps
       end associate
       deallocate(var)
    end select

  end function eos_we_saturation_difference

!------------------------------------------------------------------------

end module eos_we_module
