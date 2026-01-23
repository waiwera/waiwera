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
  use IAPWS_module
  use fluid_module

  implicit none
  private

  type, public, extends(eos_we_type) :: eos_se_type
     !! Pure supercritical water and energy equation of state type.
     private
     PetscInt, public :: region3_phase(4) = [1, 2, -1, 3] !! Map phase composition to phase index in region 3
     PetscBool, public :: pressure_conditions !! Option to allow region 3 initial and boundary conditions specified with pressure instead of density
     class(fluid_modifier_type), allocatable, public :: relative_permeability_modifier !! Modifies effective relative permeability for temperature effects
   contains
     private
     procedure, public :: init => eos_se_init
     procedure, public :: init_relative_permeability_modifier => &
          eos_se_init_relative_permeability_modifier
     procedure, public :: destroy => eos_se_destroy
     procedure, public :: region_1_transitions => eos_se_region_1_transitions
     procedure, public :: region_2_transitions => eos_se_region_2_transitions
     procedure, public :: region_3_transitions => eos_se_region_3_transitions
     procedure, public :: region_4_transitions => eos_se_region_4_transitions
     procedure, public :: transition => eos_se_transition
     procedure, public :: transition_to_single_phase => eos_se_transition_to_single_phase
     procedure, public :: transition_single_phase_to_region3 => eos_se_transition_single_phase_to_region3
     procedure, public :: transition_region4_to_supercritical => eos_se_transition_region4_to_supercritical
     procedure, public :: fluid_properties => eos_se_fluid_properties
     procedure, public :: region_1_fluid_properties => eos_se_region_1_fluid_properties
     procedure, public :: region_2_fluid_properties => eos_se_region_2_fluid_properties
     procedure, public :: region_3_fluid_properties => eos_se_region_3_fluid_properties
     procedure, public :: region_4_fluid_properties => eos_se_region_4_fluid_properties
     procedure, public :: primary_variables => eos_se_primary_variables
     procedure, public :: phase_saturations => eos_se_phase_saturations
     procedure, public :: check_primary_variables => eos_se_check_primary_variables
     procedure, public :: convert_fluid => eos_se_convert_fluid
     procedure, public :: process_conditions => eos_se_process_conditions
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
    use utils_module, only: str_to_lower

    class(eos_se_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    procedure(root_finder_routine), pointer :: fs
    PetscReal :: pressure_scale, temperature_scale, density_scale
    character(10) :: conditions
    PetscReal, parameter :: default_pressure = 1.0e5_dp
    PetscReal, parameter :: default_temperature = 20._dp ! deg C
    PetscReal, parameter :: default_pressure_scale = 1.e6_dp !! Default scale factor for non-dimensionalising pressure
    PetscReal, parameter :: default_temperature_scale = 1.e2_dp !! Default scale factor for non-dimensionalising temperature
    PetscReal, parameter :: default_density_scale = critical%density !! Default scale factor for non-dimensionalising density
    character(10), parameter :: default_conditions = "density"

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
         "liquid_density       ", "vapour_density       ", &
         "supercritical_density"]
    self%default_output_fluid_fields = [ &
         "pressure             ", "temperature          ", &
         "region               ", "vapour_saturation    ", &
         "liquid_density       ", "vapour_density       ", &
         "supercritical_density", "liquidlike_fraction  ", &
         "supercritical_phases "]

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

    fs => eos_we_saturation_difference
    allocate(primary_variable_interpolator_type :: self%primary_variable_interpolator)
    call self%init_line_finder(self%saturation_line_finder, &
         self%primary_variable_interpolator, fs, init_interpolator = PETSC_TRUE)

    call fson_get_mpi(json, "eos.conditions", default_conditions, &
         conditions, logfile)
    self%pressure_conditions = (str_to_lower(conditions) == "pressure")

    call self%init_relative_permeability_modifier(json, logfile)

  end subroutine eos_se_init

!------------------------------------------------------------------------

  subroutine eos_se_init_relative_permeability_modifier(self, json, logfile)
    !! Initialise relative permeability modifier from JSON.

    use fson
    use fson_mpi_module, only: fson_get_mpi, fson_has_mpi, fson_type_mpi
    use fson_value_m, only: TYPE_OBJECT, TYPE_NULL
    use utils_module, only: str_to_lower
    use logfile_module

    class(eos_se_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    character(max_fluid_modifier_name_length), parameter :: &
         default_relative_permeability_modifier_type_name = "linear"
    PetscInt :: modifier_type
    character(max_fluid_modifier_name_length) :: relative_permeability_modifier_type_name
    type(fson_value), pointer :: rperm_json
    PetscReal :: default_min_temperature

    modifier_type = fson_type_mpi(json, "eos.relative_permeability_modifier")
    select case (modifier_type)
    case (TYPE_OBJECT)
       call fson_get_mpi(json, "eos.relative_permeability_modifier.type", &
            default_relative_permeability_modifier_type_name, &
            relative_permeability_modifier_type_name, logfile)
       select case (str_to_lower(relative_permeability_modifier_type_name))
       case ("linear")
          call init_linear_modifier()
       case default ! null modifier
          allocate(fluid_modifier_type :: self%relative_permeability_modifier)
       end select
    case (TYPE_NULL)
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
               str_key = "eos.relative_permeability_modifier.type", &
               str_value = default_relative_permeability_modifier_type_name)
       end if
       call init_linear_modifier()
    end select
    if (fson_has_mpi(json, "eos.relative_permeability_modifier")) then
       call fson_get_mpi(json, "eos.relative_permeability_modifier", rperm_json)
    else
       rperm_json => null()
    end if
    call self%relative_permeability_modifier%init(rperm_json, logfile)

  contains

    subroutine init_linear_modifier()

      allocate(fluid_relative_permeability_linear_temperature_type :: &
           self%relative_permeability_modifier)
      select type (modifier => self%relative_permeability_modifier)
      type is (fluid_relative_permeability_linear_temperature_type)
         select type (thermo => self%thermo)
         type is (IAPWS_type)
            default_min_temperature = thermo%temperature_bdy_1_3
         end select
         call fson_get_mpi(json, &
              "eos.relative_permeability_modifier.minimum_temperature", &
              default_min_temperature, modifier%min_temperature, logfile)
         modifier%max_temperature = self%thermo%critical%temperature
      end select

    end subroutine init_linear_modifier

  end subroutine eos_se_init_relative_permeability_modifier

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

    call self%relative_permeability_modifier%destroy()

  end subroutine eos_se_destroy

!------------------------------------------------------------------------

  subroutine eos_se_transition_single_phase_to_region3(self, primary, &
       fluid, transition, err)
    !! For eos_se, carry out transition from single-phase (region 1 or
    !! 2) to region 3.

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: water_pressure, density
    PetscInt :: region

    err = 0
    call self%enforce_consistency(primary)
    region = nint(fluid%region)
    call self%water_pressure(primary, region, PETSC_FALSE, water_pressure, err)
    if (err == 0) then
       select type (region3 => self%thermo%region(3)%ptr)
       type is (IAPWS_region3_type)
          associate (temperature => primary(2))
            call region3%density([water_pressure, temperature], density, &
                 err, polish = PETSC_TRUE)
          end associate
       end select
    end if

    if (err == 0) then
       fluid%region = dble(3)
       primary(1) = density
       transition = PETSC_TRUE
    end if

  end subroutine eos_se_transition_single_phase_to_region3

!------------------------------------------------------------------------

  subroutine eos_se_transition_region4_to_supercritical(self, primary, fluid, &
       transition, err)
      !! For eos_se, make transition from region 4 to supercritical
      !! region 3. The pressure (> critical pressure) is retained,
      !! while temperature is interpolated between the bounds of the
      !! Widom delta at that pressure, according to the old vapour
      !! saturation (effectively used as an estimate of 1 - liquidlike
      !! fraction).

    use fluid_module, only: fluid_type
    use utils_module, only: hermite_spline_inv_00

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: delta(2), xi, t_delta, density_delta, water_pressure
    PetscReal, parameter :: small = 1.e-3_dp

    err = 0
    select type (region3 => self%thermo%region(3)%ptr)
    type is (IAPWS_region3_type)

       call self%water_pressure(primary, 4, PETSC_FALSE, water_pressure, err)
       associate (Sl => 1._dp - primary(2))
         select type (thermo => self%thermo)
         type is (IAPWS_type)
            call thermo%widom_delta(water_pressure, delta, err)
         end select
         if (err == 0) then
            xi = hermite_spline_inv_00(Sl)
            t_delta = (1._dp - xi) * delta(1) + xi * delta(2)
            call region3%density([water_pressure, t_delta], density_delta, &
                 err, polish = PETSC_TRUE)
         end if
       end associate

       associate (density => primary(1), temperature => primary(2))
         if (err == 0) then
            density = density_delta
            temperature = t_delta
         else ! fallback
            density = self%thermo%critical%density
            temperature = (1._dp + small) * self%thermo%critical%temperature
            err = 0
         end if
         fluid%region = dble(3)
         transition = PETSC_TRUE
       end associate

    end select

  end subroutine eos_se_transition_region4_to_supercritical

!------------------------------------------------------------------------

  subroutine eos_se_transition_to_single_phase(self, old_primary, old_fluid, &
       new_region, primary, fluid, transition, err)
    !! For eos_se, make transition from two-phase to single-phase with
    !! specified region (this may be modified to region 3 for high
    !! pressures).

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in out) :: self
    type(fluid_type), intent(in) :: old_fluid
    PetscInt, intent(in) :: new_region !! Default new region (1, 2)
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: old_saturation_pressure, pressure_factor
    PetscReal :: saturation_bound, xi
    PetscReal :: interpolated_primary(self%num_primary_variables)
    PetscReal :: old_fluid_primary(self%num_primary_variables)
    PetscReal :: interpolated_water_pressure, water_pressure
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

    select type (thermo => self%thermo)
    type is (IAPWS_type)

       self%primary_variable_interpolator%val(:, 1) = old_primary
       self%primary_variable_interpolator%val(:, 2) = primary
       call self%primary_variable_interpolator%set_index(1)
       call self%primary_variable_interpolator%find_component_at_index(&
            saturation_bound, 2, xi, err)

       if (err == 0) then

          interpolated_primary = self%primary_variable_interpolator%interpolate(xi)
          call self%water_pressure(interpolated_primary, 4, PETSC_FALSE, &
               interpolated_water_pressure, err)
          if (err == 0) then
             associate (interpolated_pressure => interpolated_primary(1))

               if (interpolated_water_pressure > thermo%critical%pressure) then

                  call self%transition_region4_to_supercritical(primary, fluid, &
                       transition, err)

               else

                  primary = interpolated_primary
                  associate (temperature => primary(2))
                    call thermo%saturation%temperature(interpolated_water_pressure, &
                         temperature, err)
                  end associate

                  if (err == 0) then
                     if (interpolated_water_pressure > thermo%saturation_pressure_bdy_1_3) then
                        call region4_above_bdy_1_3_transitions()
                     else
                        water_pressure = pressure_factor * interpolated_water_pressure
                        call self%set_water_pressure(water_pressure, primary)
                        fluid%region = dble(new_region)
                        transition = PETSC_TRUE
                     end if
                  end if
               end if
             end associate
          end if
       end if

       if (err > 0) then ! fallback

         call self%primary_variables(old_fluid, old_fluid_primary)
         call self%saturation_pressure(old_fluid_primary, nint(old_fluid%region), &
              old_saturation_pressure, err)

         if (err == 0) then
            associate(temperature => primary(2))
              temperature = old_fluid%temperature
            end associate
            if (old_saturation_pressure <= thermo%saturation_pressure_bdy_1_3) then
               water_pressure = pressure_factor * old_saturation_pressure
               call self%set_water_pressure(water_pressure, primary)
               fluid%region = dble(new_region)
               transition = PETSC_TRUE
            else
               call region4_above_bdy_1_3_transitions()
            end if
         end if

       end if

    end select

  contains

!........................................................................

    subroutine region4_above_bdy_1_3_transitions()
      !! Transitions from region 4 to subcritical regions 1, 2 or
      !! 3. It is assumed the primary array contains pressures and
      !! temperatures.

      ! Locals:
      PetscReal :: old_density, interpolated_temperature, water_param(2)
      PetscReal :: interpolated_density, old_temperature
      PetscReal :: old_component_density(old_fluid%num_components)
      PetscReal :: boundary_pressure, vapour_props(2), water_pressure
      PetscBool :: liquid
      PetscReal, parameter :: small = 1.e-6_dp

      liquid = (new_region == 1)
      old_component_density = old_fluid%component_density()
      old_density = old_component_density(1)
      old_temperature = old_fluid%temperature
      call self%thermo%saturation%temperature(interpolated_water_pressure, &
           interpolated_temperature, err)
      if (err == 0) then
         water_param = [interpolated_water_pressure, interpolated_temperature]
         select type (region3 => self%thermo%region(3)%ptr)
         type is (IAPWS_region3_type)
            call region3%saturation_density(water_param, liquid, &
                 interpolated_density, err, polish = PETSC_TRUE)
         end select
         if (err == 0) then
            associate (density => primary(1), temperature => primary(2))

              density = interpolated_density + small * &
                   (interpolated_density - old_density)
              temperature = interpolated_temperature + small * &
                   (interpolated_temperature - old_temperature)

              if (liquid) then
                 fluid%region = dble(3)
                 transition = PETSC_TRUE
              else
                 select type (thermo => self%thermo)
                 type is (IAPWS_type)
                    call thermo%boundary23%pressure(temperature, boundary_pressure)
                    call thermo%region(2)%ptr%properties([boundary_pressure, &
                         temperature], vapour_props, err)
                 end select
                 if (err == 0) then
                    associate (boundary_23_density => vapour_props(1))
                      if (density < boundary_23_density) then
                         select type (region2 => self%thermo%region(2)%ptr)
                         type is (IAPWS_region2_type)
                            water_pressure = boundary_pressure
                            call region2%pressure(primary, water_pressure, err)
                         end select
                         if (err == 0) then
                            call self%set_water_pressure(water_pressure, primary)
                            fluid%region = dble(2)
                            transition = PETSC_TRUE
                         end if
                      else
                         fluid%region = dble(3)
                         transition = PETSC_TRUE
                      end if
                    end associate
                 end if

              end if
            end associate
         end if
      end if

    end subroutine region4_above_bdy_1_3_transitions

  end subroutine eos_se_transition_to_single_phase

!------------------------------------------------------------------------

  subroutine eos_se_region_1_transitions(self, old_primary, primary, &
       old_fluid, fluid, transition, err)
    !! For eos_se, carry out phase transitions from region 1 to 3 or
    !! 4.

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in) :: old_fluid
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: water_pressure, saturation_pressure, xi
    PetscInt :: old_region

    old_region = nint(old_fluid%region)
    call self%water_pressure(primary, 1, PETSC_FALSE, water_pressure, err)
    if (err == 0) then

       associate (temperature => primary(2))
         select type (thermo => self%thermo)
         type is (IAPWS_type)

            if (temperature <= thermo%critical%temperature) then

               call self%saturation_pressure(primary, old_region, &
                    saturation_pressure, err)
               if (err == 0) then
                  if (water_pressure < saturation_pressure) then
                     call self%transition_to_two_phase(saturation_pressure, &
                          old_primary, old_fluid, primary, fluid, transition, err)
                  else if (temperature > thermo%temperature_bdy_1_3) then
                     call self%transition_single_phase_to_region3(primary, &
                          fluid, transition, err)
                  end if
               end if

            else

               if (water_pressure > thermo%critical%pressure) then
                  call self%transition_single_phase_to_region3(primary, &
                      fluid, transition, err)
               else

                  self%primary_variable_interpolator%val(:, 1) = old_primary
                  self%primary_variable_interpolator%val(:, 2) = primary
                  call self%primary_variable_interpolator%set_index(1)
                  call self%primary_variable_interpolator%find_component_at_index(&
                       thermo%critical%temperature, 2, xi, err)
                  if (err == 0) then

                     primary = self%primary_variable_interpolator%interpolate(xi)
                     call self%water_pressure(primary, 1, PETSC_FALSE, &
                          water_pressure, err)
                     if (err == 0) then
                        if (water_pressure < thermo%critical%pressure) then
                           call self%transition_to_two_phase(water_pressure, &
                                old_primary, old_fluid, primary, fluid, transition, err)
                        else
                           call self%transition_single_phase_to_region3(primary, &
                                fluid, transition, err)
                        end if
                     end if

                  end if
               end if
            end if

         end select
       end associate
    end if

  end subroutine eos_se_region_1_transitions

!------------------------------------------------------------------------

  subroutine eos_se_region_2_transitions(self, old_primary, primary, &
       old_fluid, fluid, transition, err)
    !! For eos_se, carry out phase transitions from region 2 to 3 or
    !! 4.

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in) :: old_fluid
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: water_pressure, saturation_pressure, pressure_bdy_2_3
    PetscInt :: old_region

    old_region = nint(old_fluid%region)
    call self%water_pressure(primary, 2, PETSC_FALSE, water_pressure, err)
    if (err == 0) then

       associate (temperature => primary(2))
         select type (thermo => self%thermo)
         type is (IAPWS_type)

            if (temperature <= thermo%critical%temperature) then

               call self%saturation_pressure(primary, old_region, &
                    saturation_pressure, err)
               if (err == 0) then

                  if (water_pressure > saturation_pressure) then

                     call self%transition_to_two_phase(saturation_pressure, &
                          old_primary, old_fluid, primary, fluid, transition, err)

                  else if (temperature > thermo%temperature_bdy_1_3) then

                     call thermo%boundary23%pressure(temperature, pressure_bdy_2_3)
                     if (water_pressure > pressure_bdy_2_3) then
                        call self%transition_single_phase_to_region3(primary, &
                             fluid, transition, err)
                     end if

                  end if
               end if

            else
               select type (region3 => thermo%region(3)%ptr)
               type is (IAPWS_region3_type)

                  call thermo%boundary23%pressure(temperature, pressure_bdy_2_3)

                  if (water_pressure > pressure_bdy_2_3) then

                     if (water_pressure > thermo%critical%pressure) then
                        call self%transition_single_phase_to_region3(primary, &
                             fluid, transition, err)
                     else
                        call self%transition_single_phase_to_region3(primary, &
                             fluid, transition, err)
                     end if

                  end if

               end select
            end if

         end select
       end associate

    end if

  end subroutine eos_se_region_2_transitions

!------------------------------------------------------------------------

  subroutine eos_se_region_3_transitions(self, old_primary, primary, &
       old_fluid, fluid, transition, err)
    !! For eos_se, carry out phase transitions from region 3 to 1, 2
    !! or 4.

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in) :: old_fluid
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err

    associate (temperature => primary(2))
      select type (thermo => self%thermo)
      type is (IAPWS_type)
         if (temperature < thermo%temperature_bdy_1_3) then
            call region3_to_below_bdy_1_3_transitions()
         else
            call region3_to_above_bdy_1_3_transitions()
         end if
      end select
    end associate

  contains

!........................................................................

    subroutine region3_to_below_bdy_1_3_transitions()
      !! Transitions from region 3 to temperatures below the region
      !! 1/3 boundary, in region 1, 2 or 4.

      ! Locals:
      PetscReal :: liquid_props(2), vapour_props(2)
      PetscReal :: saturation_pressure, Sv
      PetscInt :: old_region

      old_region = nint(old_fluid%region)

      associate (density => primary(1), temperature => primary(2), &
           liquid_density => liquid_props(1), vapour_density => vapour_props(1))

        call self%saturation_pressure(primary, old_region, saturation_pressure, err)
        if (err == 0) then

           call self%thermo%region(1)%ptr%properties([saturation_pressure, &
                temperature], liquid_props, err)
           if (err == 0) then
              call self%thermo%region(2)%ptr%properties([saturation_pressure, &
                   temperature], vapour_props, err)
              if (err == 0) then

                 if (density > liquid_density) then

                    call region3_to_single_phase_transitions(1)

                 else if (density > vapour_density) then

                    fluid%region = dble(4)
                    Sv = (liquid_density - density) / (liquid_density - vapour_density)
                    call self%set_water_pressure(saturation_pressure, primary)
                    primary(2) = Sv
                    transition = PETSC_TRUE

                 else

                    call region3_to_single_phase_transitions(2)

                 end if

              end if
           end if
        end if
      end associate

    end subroutine region3_to_below_bdy_1_3_transitions

!........................................................................

    subroutine region3_to_two_phase_above_bdy_1_3_transitions()
      !! Transitions from region 3 to two-phase, at temperatures above
      !! the region 1/3 boundary.

      ! Locals:
      PetscReal :: saturation_pressure, liquid_density, vapour_density
      PetscReal :: density_difference, Sv
      PetscReal, parameter :: small = 1.e-6_dp

      associate (density => primary(1), temperature => primary(2))
        select type (thermo => self%thermo)
        type is (IAPWS_type)

           if ((temperature <= thermo%critical%temperature) .and. &
                (density <= thermo%min_liquid_density_bdy_1_3)) then

              call self%saturation_pressure([old_fluid%pressure, temperature], &
                   nint(old_fluid%region), saturation_pressure, err)

              if (err == 0) then

                 select type (region3 => thermo%region(3)%ptr)
                 type is (IAPWS_region3_type)

                    call region3%saturation_density([saturation_pressure, &
                         temperature], PETSC_TRUE, liquid_density, &
                         err, polish = PETSC_TRUE)
                    if (err == 0) then
                       if (density < liquid_density) then
                          call region3%saturation_density([saturation_pressure, &
                               temperature], PETSC_FALSE, vapour_density, &
                               err, polish = PETSC_TRUE)
                          if (err == 0) then
                             if (density > vapour_density) then
                                density_difference = liquid_density - vapour_density
                                if (density_difference > small) then
                                   Sv = (liquid_density - density) / density_difference
                                else
                                   Sv = 0.5_dp
                                end if
                                fluid%region = dble(4)
                                call self%set_water_pressure(saturation_pressure, primary)
                                primary(2) = Sv
                                transition = PETSC_TRUE
                             end if
                          end if
                       end if
                    end if
                 end select

              end if

           end if

        end select
      end associate

    end subroutine region3_to_two_phase_above_bdy_1_3_transitions

!........................................................................

    subroutine region3_to_above_bdy_1_3_transitions()
      !! Transitions from region 3 to temperatures above the region
      !! 1/3 boundary, in region 2 or 4.

      ! Locals:
      PetscReal :: pressure_bdy_2_3, props(2)
      PetscReal, parameter :: eps = 1.e-6_dp

      associate (density => primary(1), temperature => primary(2))
        select type (thermo => self%thermo)
        type is (IAPWS_type)

           call thermo%boundary23%pressure(temperature, pressure_bdy_2_3)
           pressure_bdy_2_3 = min(pressure_bdy_2_3, thermo%max_pressure)
           call thermo%region(2)%ptr%properties([pressure_bdy_2_3, temperature], &
                props, err)
           if (err == 0) then
              associate(density_bdy_2_3 => props(1))

                if (density < density_bdy_2_3) then

                   call region3_to_single_phase_transitions(2)

                   if (err > 0) then
                      density = (1._dp - eps) * density_bdy_2_3
                      call region3_to_single_phase_transitions(2)
                   end if

                else
                   call region3_to_two_phase_above_bdy_1_3_transitions()
                end if
              end associate

           end if

        end select
      end associate

    end subroutine region3_to_above_bdy_1_3_transitions

!------------------------------------------------------------------------

    subroutine region3_to_single_phase_transitions(new_region)
      !! Transitions from region 3 to the single-phase new_region (1 or
      !! 2).

      PetscInt, intent(in) :: new_region
      ! Locals:
      PetscReal :: water_pressure
      PetscReal :: old_fluid_primary(self%num_primary_variables)

      err = 0
      select type (region => self%thermo%region(new_region)%ptr)
      class is (IAPWS_region_type)
         old_fluid_primary = old_primary
         old_fluid_primary(1) = old_fluid%pressure
         call self%water_pressure(old_fluid_primary, new_region, &
              PETSC_FALSE, water_pressure, err)
         if (err == 0) then
            call region%pressure(primary(1:2), water_pressure, err)
         end if
      end select

      if (err == 0) then
         fluid%region = dble(new_region)
         call self%set_water_pressure(water_pressure, primary)
         transition = PETSC_TRUE
      end if

    end subroutine region3_to_single_phase_transitions

  end subroutine eos_se_region_3_transitions

!------------------------------------------------------------------------

  subroutine eos_se_region_4_transitions(self, old_primary, primary, &
       old_fluid, fluid, transition, err)
    !! For eos_se, carry out phase transitions from region 4 to 1, 2
    !! or 3.

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in) :: old_fluid
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: water_pressure

    associate (vapour_saturation => primary(2))
      call self%water_pressure(primary, 4, PETSC_FALSE, water_pressure, err)
      if (err == 0) then
         if (vapour_saturation < 0._dp) then
            call self%transition_to_single_phase(old_primary, old_fluid, &
                 1, primary, fluid, transition, err)
         else if (vapour_saturation > 1._dp) then
            call self%transition_to_single_phase(old_primary, old_fluid, &
                 2, primary, fluid, transition, err)
         else if (water_pressure > self%thermo%critical%pressure) then
            call self%transition_region4_to_supercritical(primary, fluid, &
                 transition, err)
         end if
      end if
    end associate

  end subroutine eos_se_region_4_transitions

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

    err = 0
    transition = PETSC_FALSE
    old_region = nint(old_fluid%region)

    select case (old_region)
    case (1)
       call self%region_1_transitions(old_primary, primary, &
            old_fluid, fluid, transition, err)
    case (2)
       call self%region_2_transitions(old_primary, primary, &
            old_fluid, fluid, transition, err)
    case (3)
       call self%region_3_transitions(old_primary, primary, &
            old_fluid, fluid, transition, err)
    case (4)
       call self%region_4_transitions(old_primary, primary, &
            old_fluid, fluid, transition, err)
    end select

  end subroutine eos_se_transition

!........................................................................

  subroutine eos_se_region_1_fluid_properties(self, primary, rock, fluid, err)
    !! Calculate region 1 fluid properties from region and primary
    !! variables for pure supercritical water and energy EOS.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err

    call self%eos_we_type%fluid_properties(primary, rock, fluid, err)
    call fluid%phase(3)%zero()

  end subroutine eos_se_region_1_fluid_properties

!------------------------------------------------------------------------

  subroutine eos_se_region_2_fluid_properties(self, primary, rock, fluid, err)
    !! Calculate region 2 fluid properties from region and primary
    !! variables for pure supercritical water and energy EOS.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err

    err = 0
    call self%eos_we_type%bulk_properties(primary, fluid, err)

    if (err == 0) then
       if (fluid%is_supercritical()) then
          call region_2_supercritical_phase_properties()
       else
          call self%eos_we_type%phase_properties(primary, rock, fluid, err)
          call fluid%phase(3)%zero()
       end if
    end if

  contains

!........................................................................

    subroutine region_2_supercritical_phase_properties()
    !! Calculate region 2 supercritical phase properties from region
    !! and primary variables for pure supercritical water and energy
    !! EOS.

    ! Locals:
    PetscInt :: p, pseudo_phases
    PetscReal :: water_primary(2), properties(2), pi_liq
    PetscReal, parameter :: density = 0._dp ! not used

    err = 0
    do p = 1, 2
       call fluid%phase(p)%zero()
    end do

    associate (region => self%thermo%region(2)%ptr, phase => fluid%phase(3), &
         water_pressure => water_primary(1), water_temperature => water_primary(2))

      call self%water_pressure(primary, 2, PETSC_FALSE, water_pressure, err)
      if (err == 0) then

         water_temperature = fluid%temperature

         call region%properties(water_primary, properties, err)
         if (err == 0) then

            phase%saturation = 1._dp
            phase%density = properties(1)
            phase%internal_energy = properties(2)
            phase%specific_enthalpy = phase%internal_energy + &
                 fluid%pressure / phase%density

            phase%mass_fraction(1) = 1._dp
            phase%relative_permeability = 1._dp
            phase%capillary_pressure = 0._dp

            call region%viscosity(fluid%temperature, fluid%pressure, &
                 phase%density, phase%viscosity)

            select type (thermo => self%thermo)
            type is (IAPWS_type)
               call thermo%pi_liquidlike(water_pressure, water_temperature, density, &
                    pi_liq, pseudo_phases, err)
            end select
            if (err == 0) then
               fluid%liquidlike_fraction = pi_liq
               fluid%supercritical_phases = dble(pseudo_phases)
            end if

         end if
      end if
    end associate

  end subroutine region_2_supercritical_phase_properties

  end subroutine eos_se_region_2_fluid_properties

!------------------------------------------------------------------------

  subroutine eos_se_region_3_fluid_properties(self, primary, rock, fluid, err)
    !! Calculate region 3 fluid properties from region and primary
    !! variables for pure supercritical water and energy EOS.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: p, pseudo_phases
    PetscReal :: properties(2), PT(2), sl, pi_liq
    PetscReal :: relative_permeability(2), capillary_pressure(2)

    err = 0

    associate(density => primary(1), temperature => primary(2))

      select type (region => self%thermo%region(3)%ptr)
      type is (IAPWS_region3_type)

         fluid%temperature = temperature
         call region%properties(primary, properties, err)

         if (err == 0) then

            associate(pressure => properties(1), internal_energy => properties(2))

              fluid%pressure = pressure
              PT = [pressure, temperature]
              call self%partial_pressures(PT, 3, fluid%partial_pressure, err)
              if (err == 0) then

                 fluid%permeability_factor = 1._dp

                 call self%phase_composition(fluid, err)
                 if (err == 0) then

                    do p = 1, self%num_phases
                       call fluid%phase(p)%zero()
                    end do

                    call self%phase_saturations(primary, fluid)

                    p = self%region3_phase(nint(fluid%phase_composition))
                    associate(phase => fluid%phase(p))
                      phase%saturation = 1._dp
                      phase%density = density
                      phase%internal_energy = internal_energy
                      phase%specific_enthalpy = phase%internal_energy + &
                           fluid%pressure / phase%density
                      phase%mass_fraction(1) = 1._dp
                      call region%viscosity(fluid%temperature, fluid%pressure, &
                           phase%density, phase%viscosity)

                      if (fluid%temperature <= self%thermo%critical%temperature) then
                         sl = fluid%phase(1)%saturation
                         relative_permeability = rock%relative_permeability%values(sl)
                         capillary_pressure = [rock%capillary_pressure%value(sl, &
                              fluid%temperature), 0._dp]
                         phase%relative_permeability = relative_permeability(p)
                         phase%capillary_pressure =  capillary_pressure(p)
                      else
                         phase%relative_permeability = 1._dp
                         phase%capillary_pressure =  0._dp
                      end if
                    end associate

                    select type (thermo => self%thermo)
                    type is (IAPWS_type)
                       call thermo%pi_liquidlike(pressure, temperature, density, &
                            pi_liq, pseudo_phases, err)
                    end select
                    if (err == 0) then
                       fluid%liquidlike_fraction = pi_liq
                       fluid%supercritical_phases = dble(pseudo_phases)
                    end if

                 end if
              end if
            end associate
         end if

      end select
    end associate

  end subroutine eos_se_region_3_fluid_properties

!------------------------------------------------------------------------

  subroutine eos_se_region_4_fluid_properties(self, primary, rock, fluid, err)
    !! Calculate region 4 fluid properties from region and primary
    !! variables for pure supercritical water and energy EOS.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err

    call self%eos_we_type%bulk_properties(primary, fluid, err)

    associate(pressure => primary(1))
      select type (thermo => self%thermo)
      type is (IAPWS_type)
         if (pressure <= thermo%saturation_pressure_bdy_1_3) then ! T <= 350:
            call self%eos_we_type%phase_properties(primary, rock, fluid, err)
         else
            call region4_above_bdy_1_3_phase_properties()
         end if
      end select
    end associate
    call fluid%phase(3)%zero()
    call self%relative_permeability_modifier%modify(fluid)

  contains

!........................................................................

    subroutine region4_above_bdy_1_3_phase_properties()
      !! Calculate region 4 phase properties from region and primary
      !! variables for temperatures above the region 1/3 boundary.

      ! Locals:
      PetscInt :: p, phases
      PetscReal :: density, sl, properties(2)
      PetscReal :: relative_permeability(2), capillary_pressure(2)
      PetscBool :: liquid

      err = 0

      select type (region3 => self%thermo%region(3)%ptr)
      type is (IAPWS_region3_type)

         phases = nint(fluid%phase_composition)
         sl = fluid%phase(1)%saturation
         relative_permeability = rock%relative_permeability%values(sl)
         capillary_pressure = [rock%capillary_pressure%value(sl, &
              fluid%temperature), 0._dp]

         do p = 1, 2
            associate(phase => fluid%phase(p))

              if (btest(phases, p - 1)) then

                 liquid = (p == 1)
                 call region3%saturation_density([fluid%pressure, &
                      fluid%temperature], liquid, density, err, &
                      polish = PETSC_TRUE)

                 if (err == 0) then

                    call region3%properties([density, fluid%temperature], &
                         properties, err)

                    if (err == 0) then

                       phase%density = density
                       phase%internal_energy = properties(2)
                       phase%specific_enthalpy = phase%internal_energy + &
                            fluid%pressure / phase%density

                       phase%mass_fraction(1) = 1._dp
                       phase%relative_permeability = relative_permeability(p)
                       phase%capillary_pressure = capillary_pressure(p)

                       call region3%viscosity(fluid%temperature, fluid%pressure, &
                            phase%density, phase%viscosity)

                    else
                       exit
                    end if

                 else
                    exit
                 end if

              else
                 call phase%zero()
              end if

            end associate
         end do

      end select

    end subroutine region4_above_bdy_1_3_phase_properties

  end subroutine eos_se_region_4_fluid_properties

!------------------------------------------------------------------------

  subroutine eos_se_fluid_properties(self, primary, rock, fluid, err)
    !! Calculate fluid properties from region and primary variables
    !! for pure supercritical water and energy EOS.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: region

    err = 0
    region = nint(fluid%region)

    select case(region)
    case (1)
       call self%region_1_fluid_properties(primary, rock, fluid, err)
    case (2)
       call self%region_2_fluid_properties(primary, rock, fluid, err)
    case (3)
       call self%region_3_fluid_properties(primary, rock, fluid, err)
    case (4)
       call self%region_4_fluid_properties(primary, rock, fluid, err)
    end select

  end subroutine eos_se_fluid_properties

!------------------------------------------------------------------------

  subroutine eos_se_phase_saturations(self, primary, fluid)
    !! Assigns fluid phase saturations from fluid region and primary variables.

    use fluid_module, only: fluid_type
    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    ! Locals:
    PetscInt :: region, phases, p
    PetscReal :: s(3)

    region = nint(fluid%region)
    select case (region)
    case (1)
       s = [1._dp, 0._dp, 0._dp]
    case (2)
       s = 0._dp
       if (fluid%is_supercritical()) then
          p = 3
       else
          p = 2
       end if
       s(p) = 1._dp
    case (3)
       s = 0._dp
       phases = nint(fluid%phase_composition)
       p = self%region3_phase(phases)
       s(p) = 1._dp
    case (4)
       s = [1._dp - primary(2), primary(2), 0._dp]
    end select

    do p = 1, 3
       fluid%phase(p)%saturation = s(p)
    end do

  end subroutine eos_se_phase_saturations

!------------------------------------------------------------------------

  subroutine eos_se_primary_variables(self, fluid, primary)
    !! Determine primary variables from fluid properties.

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(out) :: primary(self%num_primary_variables)
    ! Locals:
    PetscInt :: region, phases, p

    region = nint(fluid%region)
    select case (region)
    case (1, 2)
       primary(1) = fluid%pressure
       primary(2) = fluid%temperature
    case (3)
       primary(2) = fluid%temperature
       phases = self%thermo%phase_composition(region, fluid%pressure, &
            fluid%temperature)
       p = self%region3_phase(phases)
       primary(1) = fluid%phase(p)%density
    case (4)
       primary(1) = fluid%pressure
       primary(2) = fluid%phase(2)%saturation
    end select

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
    PetscReal :: p, props(2)

    changed = PETSC_FALSE
    err = 0

    region = nint(fluid%region)
    if (region == 3) then
       call self%thermo%region(region)%ptr%properties(primary, props, err)
       if (err == 0) p = props(1)
    else
       p = primary(1)
    end if

    if (err == 0) then
      if ((p < 0._dp) .or. (p > 100.e6_dp)) then
         err = 1
      else
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
   end if

  end subroutine eos_se_check_primary_variables

!------------------------------------------------------------------------

  subroutine eos_se_convert_fluid(self, fluid1, fluid2)

    !! For fluid objects on face, convert supercritical fluid to
    !! equivalent two-phase fluid for the flux calculation. This
    !! enables the flux between sub-critical and supercritical cells
    !! to be computed twice, once treating the supercritical fluid as
    !! liquid and once treating it as vapour. The effective flux is
    !! calculated as a weighted sum of the two, according to the
    !! liquidlike fraction of the supercritical fluid.

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in) :: self
    type(fluid_type), intent(in out) :: fluid1, fluid2 !! Fluid objects

    call convert_scf(fluid1)
    call convert_scf(fluid2)

  contains

    subroutine convert_scf(fluid)
      ! If fluid is supercritical, convert to equivalent two-phase
      ! representation.

      type(fluid_type), target, intent(in out) :: fluid
      ! Locals:
      PetscInt :: super_phases, p
      type(fluid_type) :: tmp

      if (fluid%is_supercritical()) then

         call tmp%init(fluid%num_components, fluid%num_phases)
         call tmp%assign(fluid%internal_data, 1)

         super_phases = nint(fluid%supercritical_phases)

         select case (super_phases)
         case (int(b'001'))

            call tmp%phase(1)%copy(fluid%phase(3))
            call tmp%phase(2)%zero()

         case (int(b'010'))

            call tmp%phase(1)%zero()
            call tmp%phase(2)%copy(fluid%phase(3))

         case (int(b'011'))

            call tmp%phase(1)%copy(fluid%phase(3))
            call tmp%phase(2)%copy(fluid%phase(3))
            tmp%phase(1)%saturation = fluid%liquidlike_fraction
            tmp%phase(2)%saturation = 1._dp - tmp%phase(1)%saturation
            do p = 1, 2
               tmp%phase(p)%relative_permeability = tmp%phase(p)%saturation
            end do

         end select

         tmp%pressure = fluid%pressure
         tmp%temperature = fluid%temperature
         tmp%phase_composition = fluid%supercritical_phases
         call tmp%phase(3)%zero()

         call fluid%assign_internal()
         call tmp%destroy()

      end if

    end subroutine convert_scf

  end subroutine eos_se_convert_fluid

!------------------------------------------------------------------------

  subroutine eos_se_process_conditions(self, primary, region, err)

    !! Carry out processing of initial and boundary conditions -
    !! allowing initialisation of region 3 with pressure and
    !! temperature.

    class(eos_se_type), intent(in) :: self
    PetscReal, intent(in out) :: primary(self%num_primary_variables) !! Primary variables
    PetscInt, intent(in) :: region !! Thermodynamic region
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: water_pressure, density
    err = 0

    if ((self%pressure_conditions) .and. (region == 3)) then
       select type (region3 => self%thermo%region(3)%ptr)
       type is (IAPWS_region3_type)
          call self%water_pressure(primary, region, PETSC_FALSE, &
               water_pressure, err)
          if (err == 0) then
             associate (temperature => primary(2))
               call region3%density([water_pressure, temperature], &
                    density, err, polish = PETSC_TRUE)
             end associate
          end if
          if (err == 0) then
             primary(1) = density
          end if
       end select
    end if

  end subroutine eos_se_process_conditions

!------------------------------------------------------------------------

end module eos_se_module
