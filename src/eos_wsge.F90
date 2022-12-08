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

module eos_wsge_module
  !! Equation of state for non-isothermal water/salt/NCG.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use eos_module
  use eos_wse_module
  use ncg_thermodynamics_module
  use root_finder_module
  use thermodynamics_module
  use salt_thermodynamics_module
  use fluid_module

  implicit none
  private

  type, public, extends(eos_wse_type) :: eos_wsge_type
     !! Water, salt, NCG and energy equation of state type.
     private
     class(ncg_thermodynamics_type), allocatable, public :: gas
   contains
     private
     procedure, public :: init => eos_wsge_init
     procedure, public :: transition => eos_wsge_transition
     procedure, public :: transition_to_single_phase => eos_wsge_transition_to_single_phase
     procedure, public :: transition_to_two_phase => eos_wsge_transition_to_two_phase
     procedure, public :: halite_transition => eos_wsge_halite_transition
     procedure, public :: phase_composition => eos_wsge_phase_composition
     procedure, public :: bulk_properties => eos_wsge_bulk_properties
     procedure, public :: phase_properties => eos_wsge_phase_properties
     procedure, public :: primary_variables => eos_wsge_primary_variables
     procedure, public :: check_primary_variables => eos_wsge_check_primary_variables
  end type eos_wsge_type

contains

!------------------------------------------------------------------------

  subroutine eos_wsge_init(self, json, thermo, logfile)
    !! Initialise water, salt, NCG and energy EOS.

    use fson
    use fson_mpi_module, only: fson_get_mpi, fson_has_mpi, fson_type_mpi
    use fson_value_m, only: TYPE_STRING, TYPE_REAL, TYPE_NULL
    use logfile_module
    use thermodynamics_module
    use utils_module, only: str_to_lower

    class(eos_wsge_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    procedure(root_finder_function), pointer :: f
    class(*), pointer :: pinterp
    PetscReal, allocatable :: data(:, :)
    PetscReal :: pressure_scale, temperature_scale, partial_pressure_scale
    PetscInt :: scale_type
    character(max_fluid_modifier_name_length) :: permeability_modifier_type_name
    type(fson_value), pointer :: perm_json
    PetscReal, parameter :: default_pressure = 1.0e5_dp
    PetscReal, parameter :: default_temperature = 20._dp ! deg C
    PetscReal, parameter :: default_salt_mass_fraction = 0._dp
    PetscReal, parameter :: default_gas_partial_pressure = 0._dp
    PetscReal, parameter :: default_pressure_scale = 1.e6_dp !! Default scale factor for non-dimensionalising pressure
    PetscReal, parameter :: default_temperature_scale = 1.e2_dp !! Default scale factor for non-dimensionalising temperature
    PetscReal, parameter :: default_partial_pressure_scale = 1.e6_dp !! Default scale factor for non-dimensionalising partial pressure
    character(max_fluid_modifier_name_length), parameter :: &
         default_permeability_modifier_type_name = "identity"

    self%name = "wsge"
    self%description = "Water, salt, non-condensible gas and energy"
    self%primary_variable_names = [ &
         "pressure                           ", &
         "temperature/vapour_saturation      ", &
         "salt_mass_fraction/solid_saturation", &
         "gas partial pressure               "]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 3
    self%num_mobile_phases = 2
    self%phase_names = ["liquid", "vapour", "solid " ]
    self%num_components = 3
    self%component_names = ["water", "salt ", "gas  "]

    self%default_primary = [default_pressure, default_temperature, &
         default_salt_mass_fraction, default_gas_partial_pressure]
    self%default_region = 1
    self%default_tracer_phase = "liquid"
    self%required_output_fluid_fields = [ &
         "pressure                 ", &
         "temperature              ", &
         "region                   ", &
         "vapour_saturation        ", &
         "liquid_salt_mass_fraction", &
         "solid_saturation         ", &
         "gas_partial_pressure     "]
    self%default_output_fluid_fields = self%required_output_fluid_fields

    call fson_get_mpi(json, "eos.primary.scale.pressure", default_pressure_scale, &
         pressure_scale, logfile)
    call fson_get_mpi(json, "eos.primary.scale.temperature", default_temperature_scale, &
         temperature_scale, logfile)
    scale_type = fson_type_mpi(json, "eos.primary.scale.partial_pressure")
    select case (scale_type)
    case (TYPE_STRING, TYPE_NULL)
       self%scale => eos_wsge_scale_adaptive
       self%unscale => eos_wsge_unscale_adaptive
       partial_pressure_scale = 0._dp
    case (TYPE_REAL)
       call fson_get_mpi(json, "eos.primary.scale.partial_pressure", &
            default_partial_pressure_scale, partial_pressure_scale, logfile)
    end select

    allocate(self%primary_scale(4, 8))
    self%primary_scale = reshape([ &
          pressure_scale, temperature_scale, 1._dp, partial_pressure_scale, &
          pressure_scale, temperature_scale, 1._dp, partial_pressure_scale, &
          0._dp, 0._dp, 0._dp, 0._dp, &
          pressure_scale, 1._dp, 1._dp, partial_pressure_scale, &
          pressure_scale, temperature_scale, 1._dp, partial_pressure_scale, &
          pressure_scale, temperature_scale, 1._dp, partial_pressure_scale, &
          0._dp, 0._dp, 0._dp, 0._dp, &
          pressure_scale, 1._dp, 1._dp, partial_pressure_scale], [4, 8])

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
    f => eos_wsge_saturation_difference
    pinterp => self%primary_variable_interpolator
    call self%saturation_line_finder%init(f, context = pinterp)

    ! Set up permeability modifier:
    call fson_get_mpi(json, "eos.permeability_modifier.type", &
         default_permeability_modifier_type_name, &
         permeability_modifier_type_name, logfile)
    select case (str_to_lower(permeability_modifier_type_name))
    case ("power")
       allocate(fluid_permeability_factor_power_type :: self%permeability_modifier)
    case ("verma-pruess")
       allocate(fluid_permeability_factor_verma_pruess_type :: self%permeability_modifier)
    case default
       allocate(fluid_permeability_factor_null_type :: self%permeability_modifier)
    end select
    if (fson_has_mpi(json, "eos.permeability_modifier")) then
       call fson_get_mpi(json, "eos.permeability_modifier", perm_json)
    else
       perm_json => null()
    end if
    call self%permeability_modifier%init(perm_json, logfile)

  end subroutine eos_wsge_init

!------------------------------------------------------------------------

  subroutine eos_wsge_transition_to_single_phase(self, old_primary, old_fluid, &
       new_region, primary, fluid, transition, err)
    !! For eos_wsge, make transition from two-phase to single-phase with
    !! specified region.

    class(eos_wsge_type), intent(in out) :: self
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
    PetscReal :: saturation_bound, xi, interpolated_brine_pressure
    PetscReal :: solid_saturation, salt_mass_fraction
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

    associate (pressure => primary(1), temperature => primary(2), &
         salt => primary(3), partial_pressure => primary(4), &
         interpolated_pressure => interpolated_primary(1), &
         interpolated_salt => interpolated_primary(3), &
         interpolated_partial_pressure => interpolated_primary(4))

      partial_pressure = max(0._dp, min(partial_pressure, pressure))
      self%primary_variable_interpolator%val(:, 1) = old_primary
      self%primary_variable_interpolator%val(:, 2) = primary
      select type (interpolator => self%primary_variable_interpolator)
      type is (eos_wse_primary_variable_interpolator_type)
         interpolator%halite = old_halite
      end select
      call self%primary_variable_interpolator%find_component_at_index( &
           saturation_bound, 2, xi, err)

      if (err == 0) then

         interpolated_primary = self%primary_variable_interpolator%interpolate(xi)
         interpolated_brine_pressure = interpolated_pressure - &
              interpolated_partial_pressure
         pressure = pressure_factor * interpolated_brine_pressure + &
              interpolated_partial_pressure
         salt = interpolated_salt
         partial_pressure = interpolated_partial_pressure
         if (old_halite) then
            call halite_solubility_two_phase(interpolated_brine_pressure, &
                 self%thermo, salt_mass_fraction, err)
         else
            salt_mass_fraction = salt
         end if
         if (err == 0) then
            call brine_saturation_temperature(interpolated_brine_pressure, &
                 salt_mass_fraction, self%thermo, temperature, err)
            if (err == 0) then
               fluid%region = dble(new_region)
               transition = PETSC_TRUE
            end if
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
               pressure = pressure_factor * old_saturation_pressure + partial_pressure
               temperature = old_fluid%temperature
               fluid%region = dble(new_region)
               transition = PETSC_TRUE
            end if
         end if
      end if

    end associate

  end subroutine eos_wsge_transition_to_single_phase

!------------------------------------------------------------------------

  subroutine eos_wsge_transition_to_two_phase(self, saturation_pressure, &
       old_primary, old_fluid, primary, fluid, transition, err)
    !! For eos_wsge, make transition from single-phase to two-phase.

    class(eos_wsge_type), intent(in out) :: self
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
         salt => primary(3), partial_pressure => primary(4), &
         interpolated_pressure => interpolated_primary(1), &
         interpolated_salt => interpolated_primary(3), &
         interpolated_partial_pressure => interpolated_primary(4))

      old_region = nint(old_fluid%region)
      old_water_region = self%water_region(old_region)
      old_halite = self%halite(old_region)

      partial_pressure = max(0._dp, min(partial_pressure, pressure))
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
         partial_pressure = interpolated_partial_pressure
      else
         pressure = saturation_pressure + partial_pressure
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

  end subroutine eos_wsge_transition_to_two_phase

!------------------------------------------------------------------------

  subroutine eos_wsge_halite_transition(self, primary, fluid, transition, err)
    !! For eos_wsge, check for halite transitions (e.g. precipitation,
    !! dissolution).

    class(eos_wsge_type), intent(in out) :: self
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(in out) :: transition
    PetscErrorCode, intent(in out) :: err
    ! Locals:
    PetscInt :: region
    PetscReal :: temperature, solubility, brine_pressure
    PetscReal :: solid_saturation, salt_mass_fraction
    PetscReal, parameter :: small = 1.e-6_dp

    region = nint(fluid%region)

    select case (region)

    case (1, 4) ! Liquid phase present without halite

       salt_mass_fraction = primary(3)
       if (region == 1) then ! liquid
          temperature = primary(2)
       else ! two-phase
          associate (pressure => primary(1), partial_pressure => primary(4))
            brine_pressure = pressure - partial_pressure
            call brine_saturation_temperature(brine_pressure, salt_mass_fraction, &
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
             associate(pressure => primary(1), partial_pressure => primary(4))
               brine_pressure = pressure - partial_pressure
               call halite_solubility_two_phase(brine_pressure, self%thermo, &
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

  end subroutine eos_wsge_halite_transition

!------------------------------------------------------------------------

  subroutine eos_wsge_transition(self, old_primary, primary, &
       old_fluid, fluid, transition, err)
    !! For eos_wsge, check primary variables for a cell and make
    !! thermodynamic region transitions if needed.

    class(eos_wsge_type), intent(in out) :: self
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in) :: old_fluid
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: old_region, old_water_region, region_offset, new_region
    PetscBool :: old_halite
    PetscReal :: solid_saturation, salt_mass_fraction
    PetscReal :: saturation_pressure, water_pressure

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
       associate (pressure => primary(1), temperature => primary(2), &
            partial_pressure => primary(3))

         if (old_halite) then
            call halite_solubility(temperature, salt_mass_fraction, err)
         else
            salt_mass_fraction = primary(3)
         end if

         if (err == 0) then

            call brine_saturation_pressure(temperature, salt_mass_fraction, &
                 self%thermo, saturation_pressure, err)

            if (err == 0) then
               water_pressure = pressure - partial_pressure
               if (((old_water_region == 1) .and. (water_pressure < saturation_pressure)) .or. &
                    ((old_water_region == 2) .and. (water_pressure > saturation_pressure))) then
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

  end subroutine eos_wsge_transition

!------------------------------------------------------------------------

  subroutine eos_wsge_phase_composition(self, fluid, err)
    !! Determines fluid phase composition from bulk properties and
    !! thermodynamic region for non-isothermal water/salt and NCG.

    class(eos_wsge_type), intent(in out) :: self
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: region, water_region, phases

    region = nint(fluid%region)
    water_region = self%water_region(region)
    phases = self%thermo%phase_composition(water_region, &
         fluid%partial_pressure(1), fluid%temperature)
    if (phases > 0) then
       fluid%phase_composition = dble(phases)
       err = 0
    else
       err = 1
    end if

  end subroutine eos_wsge_phase_composition

!------------------------------------------------------------------------

  subroutine eos_wsge_bulk_properties(self, primary, fluid, err)
    !! Calculate fluid bulk properties from region and primary variables
    !! for non-isothermal water/salt/NCG.

    class(eos_wsge_type), intent(in out) :: self
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

    associate(partial_pressure => primary(4))
      fluid%partial_pressure(1) = fluid%pressure - partial_pressure
      fluid%partial_pressure(2) = 0._dp
      fluid%partial_pressure(3) = partial_pressure
    end associate

    if (water_region == 4) then ! two-phase
       if (region == 4) then ! without halite
          salt_mass_fraction = primary(3)
       else ! with halite
          call halite_solubility_two_phase(fluid%partial_pressure(1), self%thermo, &
               salt_mass_fraction, err)
       end if
       if (err == 0) then
          call brine_saturation_temperature(fluid%partial_pressure(1), &
               salt_mass_fraction, self%thermo, fluid%temperature, err)
       end if
    else ! single-phase
       fluid%temperature = primary(2)
    end if

    if (err == 0) then
       call self%phase_composition(fluid, err)
       if (err == 0) then
          call self%phase_saturations(primary, fluid)
          call self%permeability_modifier%modify(fluid)
       end if
    end if

  end subroutine eos_wsge_bulk_properties

!------------------------------------------------------------------------

  subroutine eos_wsge_phase_properties(self, primary, rock, fluid, err)
    !! Calculate fluid phase properties from updated fluid region and primary variables.
    !! Bulk properties need to be calculated before calling this routine.

    use rock_module, only: rock_type

    class(eos_wsge_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscInt :: p, phases, region
    PetscBool :: halite
    PetscReal :: henrys_constant
    PetscReal :: properties(2), sl, ss, sl_effective, xg
    PetscReal :: salt_mass_fraction, phase_salt_mass_fraction
    PetscReal :: relative_permeability(2), capillary_pressure(2)
    PetscReal :: brine_pressure(2), energy_solution(2)
    PetscReal :: gas_properties(2), effective_gas_properties(2)
    PetscReal :: b_viscosity, b_enthalpy

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
       ss = fluid%phase(3)%saturation
       sl_effective = sl / (1._dp - ss)
       relative_permeability = rock%relative_permeability%values(sl_effective)
       capillary_pressure = 0._dp
       henrys_constant = 0._dp
       energy_solution = 0._dp

       call self%gas%properties(fluid%partial_pressure(3), fluid%temperature, &
            gas_properties, err)

       if (err == 0) then

          ! effective brine pressure in each phase:
          brine_pressure = [fluid%pressure, fluid%partial_pressure(1)]

          if (btest(phases, 0)) then
             capillary_pressure(1) = rock%capillary_pressure%value(sl_effective, &
                  fluid%temperature)
             call self%gas%henrys_constant(fluid%temperature, henrys_constant, err)
             if (err == 0) then
                call self%gas%energy_solution(fluid%temperature, henrys_constant, &
                     energy_solution(1), err)
             end if
          end if

          if (err == 0) then

             do p = 1, self%num_mobile_phases
                associate(phase => fluid%phase(p), region => self%thermo%region(p)%ptr)

                  if (btest(phases, p - 1)) then

                     if (p == 1) then
                        call brine_properties(brine_pressure(p), fluid%temperature, &
                             salt_mass_fraction, self%thermo, properties, err)
                        phase_salt_mass_fraction = salt_mass_fraction
                     else
                        call region%properties([brine_pressure(p), fluid%temperature], &
                             properties, err)
                        phase_salt_mass_fraction = 0._dp
                     end if

                     if (err == 0) then

                        call self%gas%effective_properties(gas_properties, p, &
                             effective_gas_properties)

                        associate(brine_density => properties(1), &
                             brine_internal_energy => properties(2), &
                             gas_density => effective_gas_properties(1), &
                             gas_enthalpy => effective_gas_properties(2))

                          call self%gas%mass_fraction(fluid%partial_pressure(3), &
                               fluid%temperature, p, gas_density, brine_density, &
                               henrys_constant, xg, err)

                          if (err == 0) then

                             if (p == 1) then
                                call brine_viscosity(fluid%temperature, fluid%pressure, &
                                     salt_mass_fraction, self%thermo, b_viscosity, err)
                             else
                                call region%viscosity(fluid%temperature, fluid%pressure, &
                                     phase%density, b_viscosity)
                             end if

                             call self%gas%mixture_viscosity(b_viscosity, &
                                  fluid%temperature, fluid%partial_pressure(3), xg, p, &
                                  phase%viscosity, err)

                             if (err == 0) then
                                phase%density = brine_density + gas_density
                                phase%mass_fraction = [1._dp - phase_salt_mass_fraction - xg, &
                                     phase_salt_mass_fraction, xg]
                                phase%relative_permeability = relative_permeability(p)
                                phase%capillary_pressure =  capillary_pressure(p)
                                b_enthalpy = brine_internal_energy &
                                     + brine_pressure(p) / brine_density
                                phase%specific_enthalpy = b_enthalpy * (1._dp - xg) &
                                     + (gas_enthalpy + energy_solution(p)) * xg
                                phase%internal_energy = phase%specific_enthalpy &
                                     - fluid%pressure / phase%density
                             else
                                exit
                             end if
                          else
                             exit
                          end if
                        end associate
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
                     solid_phase%mass_fraction = [0._dp, 1._dp, 0._dp]
                  end if

                end associate
             end if

          end if
       end if
    end if

  end subroutine eos_wsge_phase_properties

!------------------------------------------------------------------------

  subroutine eos_wsge_primary_variables(self, fluid, primary)
    !! Determine primary variables from fluid properties.

    class(eos_wsge_type), intent(in) :: self
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

    primary(4) = fluid%partial_pressure(3)

  end subroutine eos_wsge_primary_variables

!------------------------------------------------------------------------

   subroutine eos_wsge_check_primary_variables(self, fluid, &
       primary, changed, err)
    !! Check if primary variables are in acceptable bounds, and return error
    !! code accordingly.

    class(eos_wsge_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    PetscBool, intent(out) :: changed
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: region, water_region
    PetscReal :: p
    PetscReal, parameter :: small = 1.e-6_dp

    changed = PETSC_FALSE
    err = 0

    associate (total_pressure => primary(1), salt => primary(3), &
         partial_pressure => primary(4))

      if (total_pressure > 0._dp) then

         associate(max_partial_pressure => (1._dp - small) * total_pressure)
           if (partial_pressure > max_partial_pressure) then
              partial_pressure = max_partial_pressure
              changed = PETSC_TRUE
           else if (partial_pressure < 0._dp) then
              partial_pressure = 0._dp
              changed = PETSC_TRUE
           end if
         end associate

         if (salt < 0._dp) then
            salt = 0._dp
            changed = PETSC_TRUE
         else if (salt > 1._dp) then
            err = 1
         end if

         p = total_pressure - partial_pressure
         if (p > 100.e6_dp) then
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

      else
         err = 1
      end if

    end associate

  end subroutine eos_wsge_check_primary_variables

!------------------------------------------------------------------------

  function eos_wsge_scale_adaptive(self, primary, region) result(scaled_primary)
    !! Non-dimensionalise eos_wsge primary variables by scaling. The
    !! first three variables (pressure, temperature or saturation and
    !! salt mass fraction or solid saturation) are scaled by fixed
    !! constants. The fourth variable, NCG partial pressure, is scaled
    !! adaptively by total pressure in the cell.

    class(eos_type), intent(in) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables)
    PetscInt, intent(in) :: region
    PetscReal :: scaled_primary(self%num_primary_variables)

    scaled_primary(1:3) = primary(1:3) / self%primary_scale(1:3, region)
    associate(scaled_partial_pressure => scaled_primary(4), &
         pressure => primary(1), partial_pressure => primary(4))
      scaled_partial_pressure = partial_pressure / pressure
    end associate

  end function eos_wsge_scale_adaptive

!------------------------------------------------------------------------

  function eos_wsge_unscale_adaptive(self, scaled_primary, region) result(primary)
    !! Re-dimensionalise eos_wgse scaled primary variables.

    class(eos_type), intent(in) :: self
    PetscReal, intent(in) :: scaled_primary(self%num_primary_variables)
    PetscInt, intent(in) :: region
    PetscReal :: primary(self%num_primary_variables)

    primary(1:3) = scaled_primary(1:3) * self%primary_scale(1:3, region)
    associate(scaled_partial_pressure => scaled_primary(4), &
         pressure => primary(1), partial_pressure => primary(4))
      partial_pressure = scaled_partial_pressure * pressure
    end associate

  end function eos_wsge_unscale_adaptive

!------------------------------------------------------------------------

  PetscReal function eos_wsge_saturation_difference(x, context) result(dp)
    !! Returns difference between brine saturation pressure and
    !! brine pressure at normalised point 0 <= x <= 1 along line between
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
       associate(P => var(1), T => var(2), Pg => var(4))
         if (context%halite) then
            call halite_solubility(T, salt_mass_fraction, err)
         else
            salt_mass_fraction = var(3)
         end if
         call brine_saturation_pressure(T, salt_mass_fraction, &
              context%thermo, Ps, err)
         dp = P - Pg - Ps
       end associate
       deallocate(var)
    end select

  end function eos_wsge_saturation_difference

!------------------------------------------------------------------------

end module eos_wsge_module
