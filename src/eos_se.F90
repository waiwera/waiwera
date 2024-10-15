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

  implicit none
  private

  type, public, extends(primary_variable_interpolator_type) :: widom_delta_interpolator_type
     private
     PetscInt, public :: bdy_index !! 1 for liquid-like boundary, 2 for vapour-like
  end type widom_delta_interpolator_type

  type, public, extends(eos_we_type) :: eos_se_type
     !! Pure supercritical water and energy equation of state type.
     private
     PetscInt :: region3_phase(4) = [1, 2, -1, 3] !! Map phase composition to phase index in region 3
     type(root_finder_type), public :: widom_delta_finder
     class(primary_variable_interpolator_type), pointer, public :: widom_delta_interpolator
     type(root_finder_type), public :: two_phase_finder
   contains
     private
     procedure, public :: init => eos_se_init
     procedure, public :: destroy => eos_se_destroy
     procedure, public :: transition => eos_se_transition
     procedure, public :: transition_to_single_phase => eos_se_transition_to_single_phase
     procedure, public :: transition_single_phase_to_region3 => eos_se_transition_single_phase_to_region3
     procedure, public :: transition_region3_to_single_phase => eos_se_transition_region3_to_single_phase
     procedure, public :: transition_region3_to_two_phase => eos_se_transition_region3_to_two_phase
     procedure, public :: transition_region4_to_supercritical => eos_se_transition_region4_to_supercritical
     procedure :: set_delta_interpolator_bdy => eos_se_set_delta_interpolator_bdy
     procedure, public :: fluid_properties => eos_se_fluid_properties
     procedure :: region2_supercritical_fluid_properties => eos_se_region2_supercritical_fluid_properties
     procedure :: region2_fluid_properties => eos_se_region2_fluid_properties
     procedure :: region3_fluid_properties => eos_se_region3_fluid_properties
     procedure :: region4_above_bdy_1_3_fluid_properties => eos_se_region4_above_bdy_1_3_fluid_properties
     procedure, public :: primary_variables => eos_se_primary_variables
     procedure, public :: phase_saturations => eos_se_phase_saturations
     procedure, public :: check_primary_variables => eos_se_check_primary_variables
     procedure, public :: convert_fluid => eos_se_convert_fluid
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
    procedure(root_finder_routine), pointer :: fs, fw, ft
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
    call init_line_finder(self%saturation_line_finder, &
         self%primary_variable_interpolator, fs, init_interpolator = PETSC_TRUE)
    fw => eos_se_widom_delta_difference
    allocate(widom_delta_interpolator_type :: self%widom_delta_interpolator)
    call init_line_finder(self%widom_delta_finder, &
         self%widom_delta_interpolator, fw, init_interpolator = PETSC_TRUE)
    ft => eos_se_two_phase_difference
    call init_line_finder(self%two_phase_finder, &
         self%primary_variable_interpolator, ft, init_interpolator = PETSC_FALSE)

  contains

    subroutine init_line_finder(finder, interpolator, f, init_interpolator)
      !! Initialises line finder (and optionally interpolator) for
      !! interpolating onto saturation line or Widom delta boundaries.

      type(root_finder_type), intent(in out) :: finder
      class(primary_variable_interpolator_type), pointer, &
           intent(in out) :: interpolator
      procedure(root_finder_routine), pointer, intent(in out) :: f
      PetscBool, intent(in) :: init_interpolator
      ! Locals:
      PetscReal, allocatable :: data(:, :)
      class(*), pointer :: pinterp

      if (init_interpolator) then
         allocate(data(2, 1 + self%num_primary_variables))
         data = 0._dp
         data(:, 1) = [0._dp, 1._dp]
         call interpolator%init(data)
         deallocate(data)
      end if

      interpolator%thermo => self%thermo
      pinterp => interpolator
      call finder%init(f, context = pinterp)

    end subroutine init_line_finder

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

    call self%widom_delta_finder%destroy()
    call self%widom_delta_interpolator%destroy()
    deallocate(self%widom_delta_interpolator)

  end subroutine eos_se_destroy

!------------------------------------------------------------------------

  subroutine eos_se_transition_single_phase_to_region3(self, primary, &
       fluid, transition, err)
    !! For eos_se, carry out transition from single phase to region 3.

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: density

    err = 0
    select type (region3 => self%thermo%region(3)%ptr)
    type is (IAPWS_region3_type)
       call region3%density(primary, density, err, polish = PETSC_TRUE)
    end select

    if (err == 0) then
       fluid%region = dble(3)
       primary(1) = density
       transition = PETSC_TRUE
    end if

  end subroutine eos_se_transition_single_phase_to_region3

!------------------------------------------------------------------------

  subroutine eos_se_transition_region3_to_single_phase(self, old_fluid, &
       new_region, primary, fluid, transition, err)
    !! For eos_se, carry out transition from region 3 to the
    !! single-phase new_region (1 or 2).

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in out) :: self
    type(fluid_type), intent(in) :: old_fluid
    PetscInt, intent(in) :: new_region
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: pressure

    err = 0
    select type (region => self%thermo%region(new_region)%ptr)
    class is (IAPWS_region_type)
       pressure = old_fluid%pressure
       call region%pressure(primary, pressure, err)
    end select

    if (err == 0) then
       fluid%region = dble(new_region)
       primary(1) = pressure
       transition = PETSC_TRUE
    end if

  end subroutine eos_se_transition_region3_to_single_phase

!------------------------------------------------------------------------

  subroutine eos_se_transition_region3_to_two_phase(self, primary, &
       fluid, transition, err)
    !! For eos_se, carry out transition from region 3 to two-phase
    !! (region 4).

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: props(2), pi_liq, Sv
    PetscInt :: phases
    PetscReal, parameter :: small = 1.e-6_dp

    err = 0
    select type (region3 => self%thermo%region(3)%ptr)
    type is (IAPWS_region3_type)
       call region3%properties(primary, props, err)
       if (err == 0) then
          associate(pressure => props(1), density => primary(1), &
               temperature => primary(2))
            call region3%pi_liquidlike(pressure, temperature, density, &
                 pi_liq, phases, err)
            if (err == 0) then
               fluid%region = dble(4)
               primary(1) = pressure
               Sv = 1._dp - pi_liq
               Sv = min(max(Sv, small), 1._dp - small)
               primary(2) = Sv
               transition = PETSC_TRUE
            end if
          end associate
       end if
    end select

  end subroutine eos_se_transition_region3_to_two_phase

!------------------------------------------------------------------------

  subroutine eos_se_set_delta_interpolator_bdy(self, bdy_index)
    !! Sets Widom delta interpolator boundary type.

    class(eos_se_type), intent(in out) :: self
    PetscInt, intent(in) :: bdy_index

    select type (interpolator => self%widom_delta_interpolator)
    type is (widom_delta_interpolator_type)
       interpolator%bdy_index = bdy_index
    end select

  end subroutine eos_se_set_delta_interpolator_bdy

!------------------------------------------------------------------------

  subroutine eos_se_transition_region4_to_supercritical(self, primary, fluid, &
       transition)
      !! For eos_se, make transition from region 4 to supercritical region 3.

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    ! Locals:
    PetscReal, parameter :: small = 1.e-6_dp

    associate (density => primary(1), temperature => primary(2))
      fluid%region = dble(3)
      density = self%thermo%critical%density
      temperature = (1._dp + small) * self%thermo%critical%temperature
      transition = PETSC_TRUE
    end associate

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
       call self%primary_variable_interpolator%find_component_at_index(&
            saturation_bound, 2, xi, err)

       if (err == 0) then

          interpolated_primary = self%primary_variable_interpolator%interpolate(xi)
          associate (interpolated_pressure => interpolated_primary(1), &
               temperature => primary(2))

            if (interpolated_pressure > thermo%critical%pressure) then
               call self%transition_region4_to_supercritical(primary, fluid, &
                    transition)
            else
               associate (pressure => primary(1))
                 pressure = interpolated_pressure

                 call thermo%saturation%temperature(pressure, temperature, err)
                 if (err == 0) then

                    if (pressure <= thermo%saturation_pressure_bdy_1_3) then
                       pressure = pressure_factor * pressure
                       fluid%region = dble(new_region)
                       transition = PETSC_TRUE
                    else
                       call region4_above_bdy_1_3_transitions()
                    end if

                 end if
               end associate
            end if
          end associate

       end if

       if (err > 0) then

          call thermo%saturation%pressure(old_fluid%temperature, &
               old_saturation_pressure, err)
          if (err == 0) then
             associate(pressure => primary(1), temperature => primary(2))

               pressure = pressure_factor * old_saturation_pressure
               temperature = old_fluid%temperature

               if (pressure <= thermo%saturation_pressure_bdy_1_3) then
                  fluid%region = dble(new_region)
                  transition = PETSC_TRUE
               else
                  call region4_above_bdy_1_3_transitions()
               end if

             end associate
          end if

       end if

    end select

  contains

!........................................................................

    subroutine region4_above_bdy_1_3_transitions()
      !! Transitions from region 4 to subcritical regions 1, 2 or 3.

      ! Locals:
      PetscBool :: liquid
      PetscReal :: factor, bdy_34_density, bdy_23_pressure, props(2)
      PetscReal, parameter :: small = 1.e-6_dp

      select type (thermo => self%thermo)
      type is (IAPWS_type)

         liquid = (new_region == 1)
         if (liquid) then
            factor = 1._dp + small
         else
            factor = 1._dp - small
         end if

         select type (region3 => thermo%region(3)%ptr)
         type is (IAPWS_region3_type)
            call region3%saturation_density(primary, liquid, bdy_34_density, &
                 err, polish = PETSC_TRUE)
         end select

         if (err == 0) then
            associate (density => primary(1), temperature => primary(2))

              density = factor * bdy_34_density
              fluid%region = dble(3)
              transition = PETSC_TRUE

              if (new_region == 2) then
                 call thermo%boundary23%pressure(temperature, bdy_23_pressure)
                 call thermo%region(2)%ptr%properties([bdy_23_pressure, &
                      temperature], props, err)
                 if (err == 0) then
                    associate(pressure => primary(1), bdy_23_density => &
                         props(1))
                      if (density < bdy_23_density) then
                         fluid%region = dble(2)
                         pressure = factor * bdy_23_pressure
                      end if
                    end associate
                 end if
              end if

            end associate
         end if

      end select

    end subroutine region4_above_bdy_1_3_transitions

  end subroutine eos_se_transition_to_single_phase

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
       call region_1_transitions()
    case (2)
       call region_2_transitions()
    case (3)
       call region_3_transitions()
    case (4)
       call region_4_transitions()
    end select

  contains

!........................................................................

    subroutine region_1_to_supercritical_transitions()
      !! Transitions from region 1 to supercritical region 3

      ! Locals:
      PetscReal :: delta(2), xi

      associate (pressure => primary(1), temperature => primary(2))
        select type (region3 => self%thermo%region(3)%ptr)
        type is (IAPWS_region3_type)

           call region3%widom_delta(pressure, delta, err)

           if (err == 0) then

              if (temperature > delta(1)) then

                 call self%set_delta_interpolator_bdy(WIDOM_DELTA_BDY_LIQUID)
                 self%widom_delta_interpolator%val(:, 1) = old_primary
                 self%widom_delta_interpolator%val(:, 2) = primary
                 call self%widom_delta_finder%find()

                 if (self%widom_delta_finder%err == 0) then
                    xi = self%widom_delta_finder%root
                    primary = self%widom_delta_interpolator%interpolate(xi)
                    call self%transition_single_phase_to_region3(primary, &
                         fluid, transition, err)
                 else
                    err = 1
                 end if

              else
                 call self%transition_single_phase_to_region3(primary, &
                      fluid, transition, err)
              end if

           end if
        end select
      end associate

    end subroutine region_1_to_supercritical_transitions

!........................................................................

    subroutine region_1_transitions()
      !! Transitions from region 1 to 3 or 4

      ! Locals:
      PetscReal :: saturation_pressure

      associate (pressure => primary(1), temperature => primary(2))
        select type (thermo => self%thermo)
        type is (IAPWS_type)

           if (temperature <= thermo%temperature_bdy_1_3) then

              call thermo%saturation%pressure(temperature, &
                   saturation_pressure, err)
              if (err == 0) then
                 if (pressure < saturation_pressure) then
                    call self%transition_to_two_phase(saturation_pressure, &
                         old_primary, old_fluid, primary, fluid, transition, err)
                 end if
              end if

           else

              if (pressure > thermo%critical%pressure) then
                 call region_1_to_supercritical_transitions()
              else

                 call thermo%saturation%pressure(temperature, &
                      saturation_pressure, err)
                 if (err == 0) then

                    if (pressure < saturation_pressure) then
                       call self%transition_to_two_phase(saturation_pressure, &
                            old_primary, old_fluid, primary, fluid, transition, err)
                    else
                       call self%transition_single_phase_to_region3(primary, &
                            fluid, transition, err)

                    end if
                 end if
              end if

           end if

        end select
      end associate

    end subroutine region_1_transitions

!........................................................................

    subroutine region2_to_supercritical_transitions()
      !! Transitions from region 2 to supercritical region 3

      ! Locals:
      PetscReal :: delta(2), xi

      associate (pressure => primary(1), temperature => primary(2))
        select type (region3 => self%thermo%region(3)%ptr)
        type is (IAPWS_region3_type)

           call region3%widom_delta(pressure, delta, err)

           if (err == 0) then

              if (temperature < delta(2)) then

                 call self%set_delta_interpolator_bdy(WIDOM_DELTA_BDY_VAPOUR)
                 self%widom_delta_interpolator%val(:, 1) = old_primary
                 self%widom_delta_interpolator%val(:, 2) = primary
                 call self%widom_delta_finder%find()

                 if (self%widom_delta_finder%err == 0) then
                    xi = self%widom_delta_finder%root
                    primary = self%widom_delta_interpolator%interpolate(xi)
                    call self%transition_single_phase_to_region3(primary, &
                         fluid, transition, err)
                 else
                    err = 1
                 end if
              else
                 call self%transition_single_phase_to_region3(primary, &
                      fluid, transition, err)
              end if
           end if

        end select
      end associate

    end subroutine region2_to_supercritical_transitions

!........................................................................

    subroutine region_2_transitions()
      !! Transitions from region 2 to 3 or 4

      ! Locals:
      PetscReal :: saturation_pressure, pressure_bdy_2_3

      associate (pressure => primary(1), temperature => primary(2))
        select type (thermo => self%thermo)
        type is (IAPWS_type)

           if (temperature <= thermo%critical%temperature) then

              call thermo%saturation%pressure(temperature, &
                   saturation_pressure, err)
              if (err == 0) then

                 if (pressure > saturation_pressure) then

                    call self%transition_to_two_phase(saturation_pressure, &
                         old_primary, old_fluid, primary, fluid, transition, err)

                 else if (temperature > thermo%temperature_bdy_1_3) then

                    call thermo%boundary23%pressure(temperature, pressure_bdy_2_3)
                    if (pressure > pressure_bdy_2_3) then
                       call self%transition_single_phase_to_region3(primary, &
                            fluid, transition, err)
                    end if

                 end if
              end if

           else
              select type (region3 => thermo%region(3)%ptr)
              type is (IAPWS_region3_type)

                 call thermo%boundary23%pressure(temperature, pressure_bdy_2_3)

                 if (pressure > pressure_bdy_2_3) then

                    if (pressure > thermo%critical%pressure) then
                       call region2_to_supercritical_transitions()
                    else
                       call self%transition_single_phase_to_region3(primary, &
                            fluid, transition, err)
                    end if

                 end if

              end select
           end if

        end select
      end associate

    end subroutine region_2_transitions

!........................................................................

    subroutine interpolate_region3_two_phase_intersection(old_primary, primary, &
         fallback_density, err)
      !! Interpolates crossing point between old_primary in region 3
      !! and primary in region 4, on region 3/4 boundary. The result
      !! overwrites the primary variable.

      PetscReal, intent(in) :: old_primary(:)
      PetscReal, intent(in out) :: primary(:)
      PetscReal, intent(in) :: fallback_density !! Density to use if interpolation fails
      PetscErrorCode, intent(out) :: err
      ! Locals:
      PetscReal :: xi, temperature_bdy_3_4

      err = 0

      self%primary_variable_interpolator%val(:, 1) = old_primary
      self%primary_variable_interpolator%val(:, 2) = primary

      select type (thermo => self%thermo)
      type is (IAPWS_type)

         associate(start_primary => self%primary_variable_interpolator%val(:, 1), &
              start_density => self%primary_variable_interpolator%val(1, 1))
           if (start_density > thermo%min_liquid_density_bdy_1_3) then
              call self%primary_variable_interpolator%find_component_at_index(&
                   thermo%min_liquid_density_bdy_1_3, 1, xi, err)
              if (err == 0) then
                 start_primary = self%primary_variable_interpolator%interpolate(xi)
              end if
           end if
         end associate

         if (err == 0) then
            call self%two_phase_finder%find()
            if (self%two_phase_finder%err == 0) then
               xi = self%two_phase_finder%root
               primary = self%primary_variable_interpolator%interpolate(xi)
            else
               err = 1
            end if
         end if

         if (err > 0) then
            call thermo%boundary34%temperature(fallback_density, &
                 temperature_bdy_3_4, err, polish = PETSC_TRUE)
            if (err == 0) then
               primary = [fallback_density, temperature_bdy_3_4]
            end if
         end if

      end select

    end subroutine interpolate_region3_two_phase_intersection

!........................................................................

    subroutine region3_to_below_bdy_1_3_transitions()
      !! Transitions from region 3 to temperatures below the region
      !! 1/3 boundary, in region 1, 2 or 4.

      ! Locals:
      PetscReal :: xi, bdy_primary(self%num_primary_variables)
      PetscReal, parameter :: eps = 1.e-6_dp

      associate (density => primary(1), temperature => primary(2))
        select type (thermo => self%thermo)
        type is (IAPWS_type)

           ! Interpolate primary at region 1/3 boundary:
           self%primary_variable_interpolator%val(:, 1) = old_primary
           self%primary_variable_interpolator%val(:, 2) = primary
           call self%primary_variable_interpolator%find_component_at_index(&
                thermo%temperature_bdy_1_3, 2, xi, err)
           bdy_primary = self%primary_variable_interpolator%interpolate(xi)

           associate (density_bdy_1_3 => bdy_primary(1))

             if (density_bdy_1_3 >= thermo%min_liquid_density_bdy_1_3) then

                call self%transition_region3_to_single_phase(old_fluid, &
                     1, primary, fluid, transition, err)

                if (err > 0) then
                   primary = bdy_primary
                   temperature = (1._dp - eps) * temperature
                   call self%transition_region3_to_single_phase(old_fluid, &
                        1, primary, fluid, transition, err)
                end if

             else if (density_bdy_1_3 <= thermo%max_vapour_density_bdy_1_3) then

                call self%transition_region3_to_single_phase(old_fluid, &
                     2, primary, fluid, transition, err)

                if (err > 0) then
                   primary = bdy_primary
                   temperature = (1._dp - eps) * temperature
                   call self%transition_region3_to_single_phase(old_fluid, &
                        2, primary, fluid, transition, err)
                end if

             else

                call interpolate_region3_two_phase_intersection(old_primary, &
                     primary, density_bdy_1_3, err)
                if (err == 0) then
                   call self%transition_region3_to_two_phase(primary, &
                        fluid, transition, err)
                end if

             end if
           end associate

        end select
      end associate

    end subroutine region3_to_below_bdy_1_3_transitions

!........................................................................

    subroutine region3_to_two_phase_above_bdy_1_3_transitions()
      !! Transitions from region 3 to two-phase, at temperatures above
      !! the region 1/3 boundary.

      ! Locals:
      PetscReal :: temperature_bdy_3_4

      associate (density => primary(1), temperature => primary(2))
        select type (thermo => self%thermo)
        type is (IAPWS_type)

           if ((temperature <= thermo%critical%temperature) .and. &
                (density <= thermo%min_liquid_density_bdy_1_3)) then

              call thermo%boundary34%temperature(density, &
                   temperature_bdy_3_4, err, polish = PETSC_TRUE)
              if (err == 0) then

                 if (temperature < temperature_bdy_3_4) then

                    call interpolate_region3_two_phase_intersection(old_primary, &
                         primary, density, err)
                    if (err == 0) then
                       call self%transition_region3_to_two_phase(primary, &
                            fluid, transition, err)
                    end if

                 end if

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

                   call self%transition_region3_to_single_phase(old_fluid, &
                        2, primary, fluid, transition, err)

                   if (err > 0) then
                      primary = [density_bdy_2_3, temperature]
                      density = (1._dp - eps) * density
                      call self%transition_region3_to_single_phase(old_fluid, &
                           2, primary, fluid, transition, err)
                   end if

                else
                   call region3_to_two_phase_above_bdy_1_3_transitions()
                end if
              end associate

           end if

        end select

      end associate
    end subroutine region3_to_above_bdy_1_3_transitions

!........................................................................

    subroutine region_3_transitions()
      !! Transitions from region 3 to 1, 2 or 4

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

    end subroutine region_3_transitions

!........................................................................

    subroutine region_4_transitions()
      !! Transitions from region 4 to 1, 2 or 3

      associate (pressure => primary(1), vapour_saturation => primary(2))
        if (vapour_saturation < 0._dp) then
           call self%transition_to_single_phase(old_primary, old_fluid, &
                1, primary, fluid, transition, err)
        else if (vapour_saturation > 1._dp) then
           call self%transition_to_single_phase(old_primary, old_fluid, &
                2, primary, fluid, transition, err)
        else if (pressure > self%thermo%critical%pressure) then
           call self%transition_region4_to_supercritical(primary, fluid, &
                transition)
        end if
      end associate

    end subroutine region_4_transitions

  end subroutine eos_se_transition

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
       call self%eos_we_type%fluid_properties(primary, rock, fluid, err)
       call fluid%phase(3)%zero()
    case (2)
       call self%region2_fluid_properties(primary, rock, fluid, err)
    case (3)
       call self%region3_fluid_properties(primary, rock, fluid, err)
    case (4)
       associate(pressure => primary(1))
         select type (thermo => self%thermo)
         type is (IAPWS_type)
            if (pressure <= thermo%saturation_pressure_bdy_1_3) then ! T <= 350:
               call self%eos_we_type%fluid_properties(primary, rock, fluid, err)
               call fluid%phase(3)%zero()
            else
               call self%region4_above_bdy_1_3_fluid_properties(primary, &
                    rock, fluid, err)
            end if
         end select
       end associate
    end select

  end subroutine eos_se_fluid_properties

!------------------------------------------------------------------------

  subroutine eos_se_region2_supercritical_fluid_properties(self, primary, &
       rock, fluid, err)
    !! Calculate region 2 supercritical fluid properties from region
    !! and primary variables for pure supercritical water and energy
    !! EOS.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: p
    PetscReal :: properties(2)

    err = 0

    do p = 1, 2
       call fluid%phase(p)%zero()
    end do

    fluid%permeability_factor = 1._dp
    call self%phase_saturations(primary, fluid)
    fluid%partial_pressure(1) = fluid%pressure
    fluid%liquidlike_fraction = 0._dp
    fluid%supercritical_phases = 2._dp

    associate (region => self%thermo%region(2)%ptr, phase => fluid%phase(3))
      call region%properties(primary, properties, err)
      if (err == 0) then

         phase%density = properties(1)
         phase%internal_energy = properties(2)
         phase%specific_enthalpy = phase%internal_energy + &
              fluid%pressure / phase%density

         phase%mass_fraction(1) = 1._dp
         phase%relative_permeability = 1._dp
         phase%capillary_pressure = 0._dp

         call region%viscosity(fluid%temperature, fluid%pressure, &
              phase%density, phase%viscosity)
      end if
    end associate

  end subroutine eos_se_region2_supercritical_fluid_properties

!------------------------------------------------------------------------

  subroutine eos_se_region2_fluid_properties(self, primary, rock, fluid, err)
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
    fluid%pressure = primary(1)
    fluid%temperature = primary(2)
    call self%phase_composition(fluid, err)

    if (err == 0) then
       if (fluid%is_supercritical()) then
          call self%region2_supercritical_fluid_properties(primary, rock, fluid, err)
       else
          call self%eos_we_type%fluid_properties(primary, rock, fluid, err)
          call fluid%phase(3)%zero()
       end if
    end if

  end subroutine eos_se_region2_fluid_properties

!------------------------------------------------------------------------

  subroutine eos_se_region3_fluid_properties(self, primary, rock, fluid, err)
    !! Calculate region 3 fluid properties from region and primary variables
    !! for pure supercritical water and energy EOS.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: p, pseudo_phases
    PetscReal :: properties(2), sl, pi_liq
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
              fluid%partial_pressure(1) = fluid%pressure
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
                      capillary_pressure = [rock%capillary_pressure%value(sl, fluid%temperature), &
                           0._dp]
                      phase%relative_permeability = relative_permeability(p)
                      phase%capillary_pressure =  capillary_pressure(p)
                   else
                      phase%relative_permeability = 1._dp
                      phase%capillary_pressure =  0._dp
                   end if

                 end associate

                 call region%pi_liquidlike(pressure, temperature, density, &
                      pi_liq, pseudo_phases, err)
                 if (err == 0) then
                    fluid%liquidlike_fraction = pi_liq
                    fluid%supercritical_phases = dble(pseudo_phases)
                 end if

              end if
            end associate
         end if

      end select
    end associate

  end subroutine eos_se_region3_fluid_properties

!------------------------------------------------------------------------

  subroutine eos_se_region4_above_bdy_1_3_fluid_properties(self, primary, &
       rock, fluid, err)
    !! Calculate region 4 fluid properties from region and primary variables
    !! for temperatures above the region 1/3 boundary.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_se_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: p, phases
    PetscReal :: density, sl, properties(2)
    PetscReal :: relative_permeability(2), capillary_pressure(2)
    PetscBool :: liquid

    err = 0

    associate(pressure => primary(1), vapour_saturation => primary(2))

      select type (region3 => self%thermo%region(3)%ptr)
      type is (IAPWS_region3_type)

         fluid%pressure = pressure
         call self%thermo%saturation%temperature(fluid%pressure, &
              fluid%temperature, err)
         if (err == 0) then

            fluid%partial_pressure(1) = fluid%pressure
            fluid%permeability_factor = 1._dp

            call self%phase_composition(fluid, err)
            if (err == 0) then
               call self%phase_saturations(primary, fluid)

               phases = nint(fluid%phase_composition)

               sl = fluid%phase(1)%saturation
               relative_permeability = rock%relative_permeability%values(sl)
               capillary_pressure = [rock%capillary_pressure%value(sl, &
                    fluid%temperature), 0._dp]

               fluid%liquidlike_fraction = sl
               fluid%supercritical_phases = 0._dp

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

               call fluid%phase(3)%zero()

            end if
         end if

      end select

    end associate

  end subroutine eos_se_region4_above_bdy_1_3_fluid_properties

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

    fluid%phase(1)%saturation = s(1)
    fluid%phase(2)%saturation = s(2)
    fluid%phase(3)%saturation = s(3)

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

    !! For fluid objects on face between sub- and super-critical
    !! cells, convert single-phase supercritical fluid to two-phase
    !! for the flux calculation.

    use fluid_module, only: fluid_type

    class(eos_se_type), intent(in) :: self
    type(fluid_type), intent(in out) :: fluid1, fluid2 !! Fluid objects
    ! Locals:
    PetscInt :: num_sc, sc_phases
    type(fluid_type), pointer :: scf
    type(fluid_type) :: tmp

    num_sc = 0
    call check(fluid1)
    call check(fluid2)

    if (num_sc == 1) then

       ! Create a temporary fluid to modify the scf internal data,
       ! while maintaining access to its original data:
       call tmp%init(scf%num_components, scf%num_phases)
       call tmp%assign(scf%internal_data, 1)

       tmp%pressure = scf%pressure
       tmp%temperature = scf%temperature

       sc_phases = nint(scf%supercritical_phases)
       call copy_phase(1)
       call copy_phase(2)
       call tmp%phase(3)%zero()

       tmp%phase_composition = scf%supercritical_phases
       tmp%phase(1)%saturation = scf%liquidlike_fraction
       tmp%phase(2)%saturation = 1._dp - tmp%phase(1)%saturation

       call scf%assign_internal()
       call tmp%destroy()

    end if

  contains

    subroutine check(fluid)
      ! If fluid is supercritical, assign scf pointer to it and
      ! increment num_sc.

      type(fluid_type), target, intent(in) :: fluid

      if (fluid%is_supercritical()) then
         scf => fluid
         num_sc = num_sc + 1
      end if

    end subroutine check

    subroutine copy_phase(i)
      ! Copy supercritical phase to phase i if it is present,
      ! otherwise zero it out.

      PetscInt, intent(in) :: i

      if (btest(sc_phases, i - 1)) then
          call tmp%phase(i)%copy(scf%phase(3))
       else
          call tmp%phase(i)%zero()
       end if

     end subroutine copy_phase

  end subroutine eos_se_convert_fluid

!------------------------------------------------------------------------

  subroutine eos_se_widom_delta_difference(x, context, f, err)
    !! Returns difference between Widom delta boundary temperature and
    !! temperature at normalised point 0 <= x <= 1 along line between
    !! start and end primary variables. Either the liquid-like or
    !! vapour-like boundary is used, based on the context%bdy_index
    !! variable.

    PetscReal, intent(in) :: x
    class(*), pointer, intent(in out) :: context
    PetscReal, intent(out) :: f
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal, allocatable :: var(:)
    PetscReal :: delta(2)

    err = 0
    select type (context)
    type is (widom_delta_interpolator_type)
       allocate(var(context%dim))
       var = context%interpolate_at_index(x)
       associate(P => var(1), T => var(2))
         select type (region3 => context%thermo%region(3)%ptr)
         type is (IAPWS_region3_type)
            call region3%widom_delta(P, delta, err)
            if (err == 0) then
               f = T - delta(context%bdy_index)
            end if
         end select
       end associate
       deallocate(var)
    end select

  end subroutine eos_se_widom_delta_difference

!------------------------------------------------------------------------

  subroutine eos_se_two_phase_difference(x, context, f, err)
    !! Returns temperature difference from region 3/4 boundary at
    !! normalised point 0 <= x <= 1 along line between start and end
    !! primary variables.

    PetscReal, intent(in) :: x
    class(*), pointer, intent(in out) :: context
    PetscReal, intent(out) :: f
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal, allocatable :: var(:)
    PetscReal :: T_bdy

    err = 0
    select type (context)
    type is (primary_variable_interpolator_type)
       allocate(var(context%dim))
       var = context%interpolate_at_index(x)
       associate(rho => var(1), T => var(2))
         select type (thermo => context%thermo)
         type is (IAPWS_type)
            call thermo%boundary34%temperature(rho, T_bdy, err, &
                 polish = PETSC_TRUE)
            if (err == 0) then
               f = T - T_bdy
            end if
         end select
       end associate
       deallocate(var)
    end select

  end subroutine eos_se_two_phase_difference

!------------------------------------------------------------------------

end module eos_se_module
