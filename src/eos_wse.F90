!   Copyright 2022 University of Auckland.

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

module eos_wse_module
  !! Equation of state for non-isothermal water/salt.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use eos_module
  use eos_we_module
  use root_finder_module
  use thermodynamics_module
  use salt_thermodynamics_module

  implicit none
  private

  type, public, extends(eos_we_type) :: eos_wse_type
     !! Water, salt and energy equation of state type.
     private
     PetscInt, allocatable :: water_region(:) !! Water region from mixture region
     PetscBool, allocatable :: halite(:) !! Halite presence from mixture region
   contains
     private
     procedure, public :: init => eos_wse_init
     procedure, public :: transition => eos_wse_transition
     procedure, public :: transition_to_single_phase => eos_wse_transition_to_single_phase
     procedure, public :: transition_to_two_phase => eos_wse_transition_to_two_phase
     procedure, public :: halite_transition => eos_wse_halite_transition
     procedure, public :: bulk_properties => eos_wse_bulk_properties
     procedure, public :: phase_properties => eos_wse_phase_properties
     procedure, public :: primary_variables => eos_wse_primary_variables
     procedure, public :: phase_saturations => eos_wse_phase_saturations
     procedure, public :: check_primary_variables => eos_wse_check_primary_variables
  end type eos_wse_type

  type, extends(primary_variable_interpolator_type) :: &
       eos_wse_primary_variable_interpolator_type
     !! Interpolator with flag for the presence of halite.
     private
     PetscBool :: halite
  end type eos_wse_primary_variable_interpolator_type

contains

!------------------------------------------------------------------------

  subroutine eos_wse_init(self, json, thermo, logfile)
    !! Initialise water, salt and energy EOS.

    use fson
    use fson_mpi_module, only: fson_get_mpi
    use logfile_module
    use thermodynamics_module

    class(eos_wse_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    procedure(root_finder_function), pointer :: f
    class(*), pointer :: pinterp
    PetscReal, allocatable :: data(:, :)
    PetscReal :: pressure_scale, temperature_scale
    PetscReal, parameter :: default_pressure = 1.0e5_dp
    PetscReal, parameter :: default_temperature = 20._dp ! deg C
    PetscReal, parameter :: default_salt_mass_fraction = 0._dp
    PetscReal, parameter :: default_pressure_scale = 1.e6_dp !! Default scale factor for non-dimensionalising pressure
    PetscReal, parameter :: default_temperature_scale = 1.e2_dp !! Default scale factor for non-dimensionalising temperature

    self%name = "wse"
    self%description = "Water, salt and energy"
    self%primary_variable_names = [ &
         "pressure                           ", &
         "temperature/vapour_saturation      ", &
         "salt_mass_fraction/solid_saturation"]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 3
    self%num_mobile_phases = 2
    self%phase_names = ["liquid", "vapour", "solid " ]
    self%num_components = 2
    self%component_names = ["water", "salt "]

    ! Mixture regions above 4 have halite present:
    self%water_region = [1, 2, 0, 4, 1, 2, 0, 4]
    self%halite = [PETSC_FALSE, PETSC_FALSE, PETSC_FALSE, PETSC_FALSE, &
         PETSC_TRUE, PETSC_TRUE, PETSC_FALSE, PETSC_TRUE]

    self%default_primary = [default_pressure, default_temperature, default_salt_mass_fraction]
    self%default_region = 1
    self%default_tracer_phase = "liquid"
    self%required_output_fluid_fields = [ &
         "pressure                 ", &
         "temperature              ", &
         "region                   ", &
         "vapour_saturation        ", &
         "liquid_salt_mass_fraction", &
         "vapour_salt_mass_fraction", &
         "solid_saturation         "]
    self%default_output_fluid_fields = self%required_output_fluid_fields

    call fson_get_mpi(json, "eos.primary.scale.pressure", default_pressure_scale, &
         pressure_scale, logfile)
    call fson_get_mpi(json, "eos.primary.scale.temperature", default_temperature_scale, &
         temperature_scale, logfile)
    allocate(self%primary_scale(3, 8))
    self%primary_scale = reshape([ &
          pressure_scale, temperature_scale, 1._dp, &
          pressure_scale, temperature_scale, 1._dp, &
          0._dp, 0._dp, 0._dp, &
          pressure_scale, 1._dp, 1._dp, &
          pressure_scale, temperature_scale, 1._dp, &
          pressure_scale, temperature_scale, 1._dp, &
          0._dp, 0._dp, 0._dp, &
          pressure_scale, 1._dp, 1._dp], [3, 8])

    self%thermo => thermo

    ! Set up saturation line finder:
    allocate(eos_wse_primary_variable_interpolator_type :: &
         self%primary_variable_interpolator)
    allocate(data(2, 1 + self%num_primary_variables))
    data = 0._dp
    data(:, 1) = [0._dp, 1._dp]
    call self%primary_variable_interpolator%init(data)
    deallocate(data)
    self%primary_variable_interpolator%thermo => self%thermo
    f => eos_wse_saturation_difference
    pinterp => self%primary_variable_interpolator
    call self%saturation_line_finder%init(f, context = pinterp)

  end subroutine eos_wse_init

!------------------------------------------------------------------------

  subroutine eos_wse_transition_to_single_phase(self, old_primary, old_fluid, &
       new_region, primary, fluid, transition, err)
    !! For eos_wse, make transition from two-phase to single-phase with
    !! specified region.

    use fluid_module, only: fluid_type

    class(eos_wse_type), intent(in out) :: self
    type(fluid_type), intent(in) :: old_fluid
    PetscInt, intent(in) :: new_region
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: old_region, new_water_region
    PetscBool :: old_halite
    PetscReal :: old_saturation_pressure, pressure_factor
    PetscReal :: saturation_bound, xi, solid_saturation, salt_mass_fraction
    PetscReal :: interpolated_primary(self%num_primary_variables)
    PetscReal, parameter :: small = 1.e-6_dp

    err = 0
    transition = PETSC_FALSE

    old_region = nint(old_fluid%region)
    old_halite = self%halite(old_region)
    if (old_halite) then
       solid_saturation = primary(3)
    else
       solid_saturation = 0._dp
    end if
    new_water_region = self%water_region(new_region)
    if (new_water_region == 1) then
       saturation_bound = 0._dp
       pressure_factor = 1._dp + small
    else
       saturation_bound = 1._dp - solid_saturation
       pressure_factor = 1._dp - small
    end if

    self%primary_variable_interpolator%val(:, 1) = old_primary
    self%primary_variable_interpolator%val(:, 2) = primary
    select type (interpolator => self%primary_variable_interpolator)
    type is (eos_wse_primary_variable_interpolator_type)
       interpolator%halite = old_halite
    end select
    call self%primary_variable_interpolator%find_component_at_index( &
         saturation_bound, 2, xi, err)

    associate (pressure => primary(1), temperature => primary(2), &
         salt => primary(3), &
         interpolated_pressure => interpolated_primary(1), &
         interpolated_salt => interpolated_primary(3))

      if (err == 0) then

         interpolated_primary = self%primary_variable_interpolator%interpolate(xi)
         pressure = pressure_factor * interpolated_pressure
         salt = interpolated_salt
         call brine_saturation_temperature(interpolated_pressure, salt_mass_fraction, &
              self%thermo, temperature, err)
         if (err == 0) then
            fluid%region = dble(new_region)
            transition = PETSC_TRUE
         end if

      else

         if (old_halite) then
            call halite_solubility(old_fluid%temperature, salt_mass_fraction, err)
         else
            salt_mass_fraction = old_primary(3)
         end if
         if (err == 0) then
            call brine_saturation_pressure(old_fluid%temperature, salt_mass_fraction, &
                 self%thermo, old_saturation_pressure, err)
            if (err == 0) then
               pressure = pressure_factor * old_saturation_pressure
               temperature = old_fluid%temperature
               fluid%region = dble(new_region)
               transition = PETSC_TRUE
            end if
         end if
      end if

    end associate

  end subroutine eos_wse_transition_to_single_phase

!------------------------------------------------------------------------

  subroutine eos_wse_transition_to_two_phase(self, saturation_pressure, &
       old_primary, old_fluid, primary, fluid, transition, err)
    !! For eos_wse, make transition from single-phase to two-phase.

    use fluid_module, only: fluid_type

    class(eos_wse_type), intent(in out) :: self
    PetscReal, intent(in) :: saturation_pressure
    type(fluid_type), intent(in) :: old_fluid
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: old_region, old_water_region, new_region
    PetscBool :: old_halite
    PetscReal :: interpolated_primary(self%num_primary_variables)
    PetscReal :: xi, solid_saturation
    PetscReal, parameter :: small = 1.e-6_dp

    err = 0
    associate (pressure => primary(1), vapour_saturation => primary(2), &
         salt => primary(3), &
         interpolated_pressure => interpolated_primary(1), &
         interpolated_salt => interpolated_primary(3))

      old_region = nint(old_fluid%region)
      old_water_region = self%water_region(old_region)
      old_halite = self%halite(old_region)

      self%primary_variable_interpolator%val(:, 1) = old_primary
      self%primary_variable_interpolator%val(:, 2) = primary
      select type (interpolator => self%primary_variable_interpolator)
      type is (eos_wse_primary_variable_interpolator_type)
         interpolator%halite = old_halite
      end select
      call self%saturation_line_finder%find()

      if (self%saturation_line_finder%err == 0) then
         xi = self%saturation_line_finder%root
         interpolated_primary = self%primary_variable_interpolator%interpolate(xi)
         pressure = interpolated_pressure
         salt = interpolated_salt
      else
         pressure = saturation_pressure
      end if

      if (old_halite) then
         solid_saturation = primary(3)
         new_region = 8
      else
         solid_saturation = 0._dp
         new_region = 4
      end if
      if (old_water_region == 1) then
         vapour_saturation = small
      else
         vapour_saturation = 1._dp - solid_saturation - small
      end if

      fluid%region = dble(new_region)
      transition = PETSC_TRUE

    end associate

  end subroutine eos_wse_transition_to_two_phase

!------------------------------------------------------------------------

  subroutine eos_wse_halite_transition(self, primary, fluid, transition, err)
    !! For eos_wse, check for halite transitions (e.g. precipitation,
    !! dissolution).

    use fluid_module, only: fluid_type

    class(eos_wse_type), intent(in out) :: self
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(in out) :: transition
    PetscErrorCode, intent(in out) :: err
    ! Locals:
    PetscInt :: region
    PetscReal :: temperature, solubility
    PetscReal :: solid_saturation, salt_mass_fraction
    PetscReal, parameter :: small = 1.e-6_dp

    region = nint(fluid%region)

    select case (region)

    case (1, 4) ! Liquid phase present without halite

       salt_mass_fraction = primary(3)
       if (region == 1) then ! liquid
          temperature = primary(2)
       else ! two-phase
          associate (pressure => primary(1))
            call brine_saturation_temperature(pressure, salt_mass_fraction, &
                 self%thermo, temperature, err)
          end associate
       end if

       if (err == 0) then
          call halite_solubility(temperature, solubility, err)
          if (salt_mass_fraction > solubility) then
             ! halite precipitates out of liquid:
             solid_saturation = small
             primary(3) = solid_saturation
             fluid%region = dble(region + 4)
             transition = PETSC_TRUE
          end if
       end if

    case (2) ! Vapour phase only without halite

       salt_mass_fraction = primary(3)
       if (salt_mass_fraction > 0._dp) then
          ! halite precipitates out of vapour:
          solid_saturation = small
          primary(3) = solid_saturation
          fluid%region = dble(6)
          transition = PETSC_TRUE
       end if

    case (5, 8) ! Liquid phase present with halite

       solid_saturation = primary(3)
       if (solid_saturation < 0._dp) then
          ! halite dissolves into liquid:

          if (region == 5) then ! liquid
             temperature = primary(2)
             call halite_solubility(temperature, solubility, err)
          else ! two-phase
             associate(pressure => primary(1))
               call halite_solubility_two_phase(pressure, self%thermo, &
                    solubility, err)
             end associate
          end if

          if (err == 0) then
             salt_mass_fraction = solubility - small
             primary(3) = salt_mass_fraction
             fluid%region = dble(region - 4)
             transition = PETSC_TRUE
          end if
       end if

    case (6) ! Vapour phase only with halite

       solid_saturation = primary(3)
       if (solid_saturation < 0._dp) then
          ! halite disappears (can't dissolve into vapour phase):
          salt_mass_fraction = 0._dp
          primary(3) = salt_mass_fraction
          fluid%region = dble(2)
          transition = PETSC_TRUE
       end if

    end select

  end subroutine eos_wse_halite_transition

!------------------------------------------------------------------------

  subroutine eos_wse_transition(self, old_primary, primary, &
       old_fluid, fluid, transition, err)
    !! For eos_wse, check primary variables for a cell and make
    !! thermodynamic region transitions if needed.

    use fluid_module, only: fluid_type

    class(eos_wse_type), intent(in out) :: self
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in) :: old_fluid
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: old_region, old_water_region, region_offset, new_region
    PetscBool :: old_halite
    PetscReal :: solid_saturation, salt_mass_fraction, saturation_pressure

    err = 0
    transition = PETSC_FALSE
    old_region = nint(old_fluid%region)
    old_water_region = self%water_region(old_region)
    old_halite = self%halite(old_region)

    if (old_water_region == 4) then  ! Two-phase
       if (old_halite) then
          region_offset = 4
       else
          region_offset = 0
       end if
       associate (vapour_saturation => primary(2))

         if (vapour_saturation < 0._dp) then
            new_region = region_offset + 1
            call self%transition_to_single_phase(old_primary, old_fluid, &
                 new_region, primary, fluid, transition, err)
         else
            if (old_halite) then
               solid_saturation = primary(3)
            else
               solid_saturation = 0._dp
            end if
            if (vapour_saturation > 1._dp - solid_saturation) then
               new_region = region_offset + 2
               call self%transition_to_single_phase(old_primary, old_fluid, &
                    new_region, primary, fluid, transition, err)
            end if
         end if

     end associate
    else  ! Single-phase
       associate (pressure => primary(1), temperature => primary(2))

         if (old_halite) then
            call halite_solubility(temperature, salt_mass_fraction, err)
         else
            salt_mass_fraction = primary(3)
         end if

         if (err == 0) then

            call brine_saturation_pressure(temperature, salt_mass_fraction, &
                 self%thermo, saturation_pressure, err)

            if (err == 0) then
               if (((old_water_region == 1) .and. (pressure < saturation_pressure)) .or. &
                    ((old_water_region == 2) .and. (pressure > saturation_pressure))) then
                  call self%transition_to_two_phase(saturation_pressure, &
                       old_primary, old_fluid, primary, fluid, transition, err)
               end if
            end if

         end if

       end associate
    end if

    if (err == 0) then
       call self%halite_transition(primary, fluid, transition, err)
    end if

  end subroutine eos_wse_transition

!------------------------------------------------------------------------

  subroutine eos_wse_phase_composition(self, fluid, err)
    !! Determines fluid phase composition from bulk properties and
    !! thermodynamic region for non-isothermal water/salt.

    use fluid_module, only: fluid_type

    class(eos_wse_type), intent(in out) :: self
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: region, water_region, phases

    region = nint(fluid%region)
    water_region = self%water_region(region)
    phases = self%thermo%phase_composition(water_region, fluid%pressure, &
         fluid%temperature)
    if (phases > 0) then
       fluid%phase_composition = dble(phases)
       err = 0
    else
       err = 1
    end if

  end subroutine eos_wse_phase_composition

!------------------------------------------------------------------------

  subroutine eos_wse_bulk_properties(self, primary, fluid, err)
    !! Calculate fluid bulk properties from region and primary variables
    !! for non-isothermal water/salt.

    use fluid_module, only: fluid_type

    class(eos_wse_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscInt :: region, water_region
    PetscReal :: salt_mass_fraction

    err = 0
    fluid%pressure = primary(1)
    region = nint(fluid%region)
    water_region = self%water_region(region)

    if (water_region == 4) then ! two-phase
       if (region == 4) then ! without halite
          salt_mass_fraction = primary(3)
       else ! with halite
          call halite_solubility_two_phase(fluid%pressure, self%thermo, &
               salt_mass_fraction, err)
       end if
       call brine_saturation_temperature(fluid%pressure, salt_mass_fraction, &
            self%thermo, fluid%temperature, err)
    else ! single-phase
       fluid%temperature = primary(2)
    end if

    if (err == 0) then
       call self%phase_composition(fluid, err)
       if (err == 0) then
          call self%phase_saturations(primary, fluid)
          fluid%partial_pressure(1) = fluid%pressure
          fluid%partial_pressure(2) = 0._dp
       end if
    end if

  end subroutine eos_wse_bulk_properties

!------------------------------------------------------------------------

  subroutine eos_wse_phase_saturations(self, primary, fluid)
    !! Assigns fluid phase saturations from fluid region and primary variables.

    use fluid_module, only: fluid_type
    class(eos_wse_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    ! Locals:
    PetscInt :: region, water_region
    PetscBool :: halite
    PetscReal :: solid_saturation, fluid_saturation

    region = nint(fluid%region)
    water_region = self%water_region(region)
    halite = self%halite(region)

    if (halite) then
       solid_saturation = primary(3)
    else
       solid_saturation = 0._dp
    end if
    fluid_saturation = 1._dp - solid_saturation

    select case (water_region)
    case (1)
       fluid%phase(1)%saturation = fluid_saturation
       fluid%phase(2)%saturation = 0._dp
    case (2)
       fluid%phase(1)%saturation = 0._dp
       fluid%phase(2)%saturation = fluid_saturation
    case (4)
       fluid%phase(1)%saturation = fluid_saturation - primary(2)
       fluid%phase(2)%saturation = primary(2)
    end select

    fluid%phase(3)%saturation = solid_saturation

  end subroutine eos_wse_phase_saturations

!------------------------------------------------------------------------

  subroutine eos_wse_phase_properties(self, primary, rock, fluid, err)
    !! Calculate fluid phase properties from updated fluid region and primary variables.
    !! Bulk properties need to be calculated before calling this routine.

    use rock_module, only: rock_type
    use fluid_module, only: fluid_type

    class(eos_wse_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscInt :: p, phases, region
    PetscBool :: halite
    PetscReal :: properties(2), sl, salt_mass_fraction, phase_salt_mass_fraction
    PetscReal :: relative_permeability(2), capillary_pressure(2)

    err = 0
    phases = nint(fluid%phase_composition)
    region = nint(fluid%region)
    halite = self%halite(region)
    if (halite) then
       call halite_solubility(fluid%temperature, salt_mass_fraction, err)
    else
       salt_mass_fraction = primary(3)
    end if

    if (err == 0) then

       sl = fluid%phase(1)%saturation
       relative_permeability = rock%relative_permeability%values(sl)
       capillary_pressure = [rock%capillary_pressure%value(sl, fluid%temperature), &
            0._dp]

       do p = 1, self%num_mobile_phases
          associate(phase => fluid%phase(p), region => self%thermo%region(p)%ptr)

            if (btest(phases, p - 1)) then

               if (p == 1) then
                  call brine_properties(fluid%pressure, fluid%temperature, &
                       salt_mass_fraction, self%thermo, properties, err)
                  phase_salt_mass_fraction = salt_mass_fraction
               else
                  call region%properties([fluid%pressure, fluid%temperature], &
                       properties, err)
                  phase_salt_mass_fraction = 0._dp
               end if

               if (err == 0) then

                  phase%density = properties(1)
                  phase%internal_energy = properties(2)
                  phase%specific_enthalpy = phase%internal_energy + &
                       fluid%pressure / phase%density

                  phase%mass_fraction = [1._dp - phase_salt_mass_fraction, &
                       phase_salt_mass_fraction]
                  phase%relative_permeability = relative_permeability(p)
                  phase%capillary_pressure =  capillary_pressure(p)

                  if (p == 1) then
                     call brine_viscosity(fluid%temperature, fluid%pressure, &
                          phase%density, salt_mass_fraction, self%thermo, &
                          phase%viscosity, err)
                  else
                     call region%viscosity(fluid%temperature, fluid%pressure, &
                          phase%density, phase%viscosity)
                  end if

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
               phase%mass_fraction = 0._dp
            end if

          end associate
       end do

       if (err == 0) then
          associate(solid_phase => fluid%phase(3))

            call halite_properties(fluid%pressure, fluid%temperature, &
                 properties, err)
            if (err == 0) then
               solid_phase%density = properties(1)
               solid_phase%internal_energy = properties(2)
               solid_phase%specific_enthalpy = solid_phase%internal_energy + &
                    fluid%pressure / solid_phase%density
               solid_phase%relative_permeability = 0._dp
               solid_phase%capillary_pressure = 0._dp
               solid_phase%viscosity = 0._dp
               solid_phase%mass_fraction = [0._dp, 1._dp]
            end if

          end associate
       end if

    end if

  end subroutine eos_wse_phase_properties

!------------------------------------------------------------------------

  subroutine eos_wse_primary_variables(self, fluid, primary)
    !! Determine primary variables from fluid properties.

    use fluid_module, only: fluid_type

    class(eos_wse_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(out) :: primary(self%num_primary_variables)
    ! Locals:
    PetscInt :: region, water_region
    PetscBool :: halite

    region = nint(fluid%region)
    water_region = self%water_region(region)
    halite = self%halite(region)

    primary(1) = fluid%pressure

    if (water_region == 4) then
       primary(2) = fluid%phase(2)%saturation
    else
       primary(2) = fluid%temperature
    end if

    if (halite) then
       primary(3) = fluid%phase(3)%saturation
    else
       primary(3) = fluid%component_mass_fraction(2)
    end if

  end subroutine eos_wse_primary_variables

!------------------------------------------------------------------------

   subroutine eos_wse_check_primary_variables(self, fluid, &
       primary, changed, err)
    !! Check if primary variables are in acceptable bounds, and return error
    !! code accordingly.

    use fluid_module, only: fluid_type

    class(eos_wse_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    PetscBool, intent(out) :: changed
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: region, water_region

    changed = PETSC_FALSE
    err = 0

    associate (p => primary(1), salt => primary(3))
      if ((p < 0._dp) .or. (p > 100.e6_dp)) then
         err = 1
      else
         if ((salt < -1._dp) .or. (salt > 2._dp)) then
            err = 1
         else
            region = nint(fluid%region)
            water_region = self%water_region(region)
            if (water_region == 4) then
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
      end if
    end associate

  end subroutine eos_wse_check_primary_variables

!------------------------------------------------------------------------

  PetscReal function eos_wse_saturation_difference(x, context) result(dp)
    !! Returns difference between brine saturation pressure and
    !! pressure at normalised point 0 <= x <= 1 along line between
    !! start and end primary variables.

    PetscReal, intent(in) :: x
    class(*), pointer, intent(in out) :: context

    ! Locals:
    PetscReal, allocatable :: var(:)
    PetscReal :: salt_mass_fraction, Ps
    PetscInt :: err

    select type (context)
    type is (eos_wse_primary_variable_interpolator_type)
       allocate(var(context%dim))
       var = context%interpolate_at_index(x)
       associate(P => var(1), T => var(2))
         if (context%halite) then
            call halite_solubility(T, salt_mass_fraction, err)
         else
            salt_mass_fraction = var(3)
         end if
         call brine_saturation_pressure(T, salt_mass_fraction, &
              context%thermo, Ps, err)
         dp = P - Ps
       end associate
       deallocate(var)
    end select

  end function eos_wse_saturation_difference

!------------------------------------------------------------------------

end module eos_wse_module
