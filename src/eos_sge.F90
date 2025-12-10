!   Copyright 2025 University of Auckland.

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

module eos_sge_module
  !! Equation of state for non-isothermal water and non-condensible
  !! gas, sub- or super-critical.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use eos_module
  use eos_se_module
  use eos_wge_module
  use root_finder_module
  use thermodynamics_module
  use IAPWS_module
  use fluid_module
  use ncg_thermodynamics_module

  implicit none
  private

  type, public, extends(eos_se_type) :: eos_sge_type
     !! Supercritical water, non-condensible gas and energy equation
     !! of state type.
     private
     class(ncg_thermodynamics_type), allocatable, public :: gas
     class(eos_wge_type), allocatable, public :: eos_wge !! To emulate multiple inheritance from eos_se_type and eos_wge_type
   contains
     private
     procedure, public :: init => eos_sge_init
     procedure, public :: destroy => eos_sge_destroy
     procedure, public :: water_pressure => eos_sge_water_pressure
     procedure, public :: set_water_pressure => eos_sge_set_water_pressure
     procedure, public :: partial_pressure_coefficient => eos_sge_partial_pressure_coefficient
     procedure, public :: partial_pressures => eos_sge_partial_pressures
     procedure, public :: effective_gas_properties => eos_sge_effective_gas_properties
     procedure, public :: enforce_consistency => eos_sge_enforce_consistency
     procedure, public :: region_1_fluid_properties => eos_sge_region_1_fluid_properties
     procedure, public :: region_2_fluid_properties => eos_sge_region_2_fluid_properties
     procedure, public :: region_3_fluid_properties => eos_sge_region_3_fluid_properties
     procedure, public :: region_4_fluid_properties => eos_sge_region_4_fluid_properties
     procedure, public :: primary_variables => eos_sge_primary_variables
     procedure, public :: check_primary_variables => eos_sge_check_primary_variables
  end type eos_sge_type

contains

!------------------------------------------------------------------------

  subroutine eos_sge_init(self, json, thermo, logfile)
    !! Initialise supercritical water, NCG and energy EOS.

    use fson
    use fson_mpi_module, only: fson_get_mpi, fson_has_mpi, fson_type_mpi
    use fson_value_m, only: TYPE_STRING, TYPE_REAL, TYPE_NULL, TYPE_OBJECT
    use logfile_module
    use thermodynamics_module
    use IAPWS_module, only: critical
    use utils_module, only: str_to_lower

    class(eos_sge_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    procedure(root_finder_routine), pointer :: fs
    PetscReal :: pressure_scale, temperature_scale, density_scale, partial_pressure_scale
    PetscInt :: scale_type
    character(10) :: conditions
    PetscReal, parameter :: default_pressure = 1.0e5_dp
    PetscReal, parameter :: default_temperature = 20._dp ! deg C
    PetscReal, parameter :: default_gas_partial_pressure = 0._dp
    PetscReal, parameter :: default_pressure_scale = 1.e6_dp !! Default scale factor for non-dimensionalising pressure
    PetscReal, parameter :: default_temperature_scale = 1.e2_dp !! Default scale factor for non-dimensionalising temperature
    PetscReal, parameter :: default_density_scale = critical%density !! Default scale factor for non-dimensionalising density
    PetscReal, parameter :: default_partial_pressure_scale = 1.e6_dp !! Default scale factor for non-dimensionalising partial pressure
    character(10), parameter :: default_conditions = "density"
    character(max_fluid_modifier_name_length), parameter :: &
         default_relative_permeability_modifier_type_name = "linear"

    self%name = "sge"
    self%description = "Supercritical water, non-condensible gas and energy"
    self%primary_variable_names = [ &
         "pressure/density             ", &
         "temperature/vapour_saturation", &
         "gas partial pressure         "]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 3
    self%num_mobile_phases = 3
    self%phase_names = ["liquid       ", "vapour       ", "supercritical"]
    self%num_components = 2
    self%component_names = ["water", "gas  "]

    self%default_primary = [default_pressure, default_temperature, &
         default_gas_partial_pressure]
    self%default_region = 1
    self%default_tracer_phase = "liquid"
    self%required_output_fluid_fields = [ &
         "pressure             ", "temperature          ", &
         "region               ", "vapour_saturation    ", &
         "liquid_density       ", "vapour_density       ", &
         "supercritical_density", "gas_partial_pressure "]
    self%default_output_fluid_fields = [ &
         "pressure             ", "temperature          ", &
         "region               ", "vapour_saturation    ", &
         "liquid_density       ", "vapour_density       ", &
         "supercritical_density", "liquidlike_fraction  ", &
         "gas_partial_pressure "]

    call fson_get_mpi(json, "eos.primary.scale.pressure", default_pressure_scale, &
         pressure_scale, logfile)
    call fson_get_mpi(json, "eos.primary.scale.temperature", default_temperature_scale, &
         temperature_scale, logfile)
    call fson_get_mpi(json, "eos.primary.scale.density", default_density_scale, &
         density_scale, logfile)
    scale_type = fson_type_mpi(json, "eos.primary.scale.partial_pressure")
    select case (scale_type)
    case (TYPE_STRING, TYPE_NULL)
       self%scale => eos_sge_scale_adaptive
       self%unscale => eos_sge_unscale_adaptive
       partial_pressure_scale = default_partial_pressure_scale
    case (TYPE_REAL)
       call fson_get_mpi(json, "eos.primary.scale.partial_pressure", &
            default_partial_pressure_scale, partial_pressure_scale, logfile)
    end select
    allocate(self%primary_scale(3, 4))
    self%primary_scale = reshape([ &
          pressure_scale, temperature_scale, partial_pressure_scale, &
          pressure_scale, temperature_scale, partial_pressure_scale, &
          density_scale, temperature_scale, partial_pressure_scale, &
          pressure_scale, 1._dp, partial_pressure_scale], [3, 4])

    self%thermo => thermo

    fs => eos_wge_saturation_difference
    allocate(primary_variable_interpolator_type :: self%primary_variable_interpolator)
    call self%init_line_finder(self%saturation_line_finder, &
         self%primary_variable_interpolator, fs, init_interpolator = PETSC_TRUE)

    call fson_get_mpi(json, "eos.conditions", default_conditions, &
         conditions, logfile)
    self%pressure_conditions = (str_to_lower(conditions) == "pressure")

    call self%init_relative_permeability_modifier(json, logfile)

  end subroutine eos_sge_init

!------------------------------------------------------------------------

  subroutine eos_sge_destroy(self)
    !! Destroy supercritical water, non-condensible gas and energy EOS.

    class(eos_sge_type), intent(in out) :: self

    call self%eos_se_type%destroy()

    if (allocated(self%gas)) then
       call self%gas%destroy()
       deallocate(self%gas)
    end if

    if (allocated(self%eos_wge)) then
       call self%eos_wge%destroy()
       deallocate(self%eos_wge)
    end if

  end subroutine eos_sge_destroy

!------------------------------------------------------------------------

  subroutine eos_sge_water_pressure(self, primary, region, &
       liquid, water_pressure, err)
    !! For eos_sge, return water pressure from primary variables, for
    !! regions 1, 2 or 4.

    class(eos_sge_type), intent(in) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables)
    PetscInt, intent(in) :: region
    PetscBool, intent(in) :: liquid
    PetscReal, intent(out) :: water_pressure
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: Pw, t, xi
    PetscReal :: pressure, partial_pressure

    err = 0
    pressure = primary(1)
    partial_pressure = primary(3)

    if (liquid) then
       select case (region)
       case (1, 2)
          associate (temperature => primary(2))
            xi = self%partial_pressure_coefficient(temperature, liquid)
          end associate
          water_pressure = pressure - xi * partial_pressure
       case (4)
          select type (thermo => self%thermo)
          type is (IAPWS_type)
             if (pressure < thermo%saturation_pressure_bdy_1_3) then
                water_pressure = pressure
             else
                Pw = pressure - partial_pressure
                call self%thermo%saturation%temperature(Pw, t, err)
                if (err == 0) then
                   xi = self%partial_pressure_coefficient(t, liquid)
                   water_pressure = pressure - xi * partial_pressure
                end if
             end if
          end select
       end select
    else
       water_pressure = pressure - partial_pressure
    end if

  end subroutine eos_sge_water_pressure

!------------------------------------------------------------------------

  subroutine eos_sge_set_water_pressure(self, water_pressure, primary)
    !! For eos_sge, update primary variables for specified water
    !! pressure (regions 1, 2, 4).

    class(eos_sge_type), intent(in) :: self
    PetscReal, intent(in) :: water_pressure
    PetscReal, intent(in out) :: primary(self%num_primary_variables)

    call self%eos_wge%set_water_pressure(water_pressure, primary)

  end subroutine eos_sge_set_water_pressure

!------------------------------------------------------------------------

  PetscReal function eos_sge_partial_pressure_coefficient(self, temperature, &
       liquid) result(coef)
    !! For eos_sge, return coefficient determining how much of the gas
    !! partial pressure is subtracted from the total pressure to
    !! calculate the effective liquid water pressure (used to
    !! calculate water properties). Specify liquid = true for
    !! sub-critical liquid phase.

    use utils_module, only: hermite_spline_01

    class(eos_sge_type), intent(in) :: self
    PetscReal, intent(in) :: temperature
    PetscBool, intent(in) :: liquid
    ! Locals:
    PetscReal :: xi

    if (liquid) then
       select type (thermo => self%thermo)
       type is (IAPWS_type)
          if (temperature < thermo%temperature_bdy_1_3) then
             coef = 0._dp
          else if (temperature < thermo%critical%temperature) then
             xi = (temperature - thermo%temperature_bdy_1_3) / &
                  (thermo%critical%temperature - thermo%temperature_bdy_1_3)
             coef = hermite_spline_01(xi)
          else
             coef = 1._dp
          end if
       end select
    else
       coef = 1._dp
    end if

  end function eos_sge_partial_pressure_coefficient

!------------------------------------------------------------------------

  subroutine eos_sge_partial_pressures(self, primary, region, &
       partial_pressures, err)
    !! Set partial pressures from primary variables.

    class(eos_sge_type), intent(in) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    PetscInt, intent(in) :: region !! Fluid region
    PetscReal, intent(out) :: partial_pressures(self%num_components) !! Partial pressures
    PetscErrorCode, intent(out) :: err !! Error code

    call self%eos_wge%partial_pressures(primary, region, &
         partial_pressures, err)

  end subroutine eos_sge_partial_pressures

!------------------------------------------------------------------------

  subroutine eos_sge_enforce_consistency(self, primary)
    !! Check internal consistency of primary variables and adjust if
    !! necessary.

    class(eos_sge_type), intent(in) :: self
    PetscReal, intent(in out) :: primary(self%num_primary_variables)

    call self%eos_wge%enforce_consistency(primary)

  end subroutine eos_sge_enforce_consistency

!........................................................................

  subroutine eos_sge_region_1_fluid_properties(self, primary, rock, fluid, err)
    !! Calculate region 1 fluid properties from region and primary
    !! variables for supercritical water, NCG and energy EOS.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_sge_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err

    call self%eos_wge%fluid_properties(primary, rock, fluid, err)
    call fluid%phase(3)%zero()

  end subroutine eos_sge_region_1_fluid_properties

!------------------------------------------------------------------------

  subroutine eos_sge_region_2_fluid_properties(self, primary, rock, fluid, err)
    !! Calculate region 2 fluid properties from region and primary
    !! variables for supercritical water, NCG and energy EOS.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_sge_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err

    err = 0
    call self%eos_wge%bulk_properties(primary, fluid, err)

    if (err == 0) then
       if (fluid%is_supercritical()) then
          call region_2_supercritical_phase_properties()
       else
          call self%eos_wge%phase_properties(primary, rock, fluid, err)
          call fluid%phase(3)%zero()
          fluid%supercritical_phases = 0._dp
       end if
    end if

  contains

!........................................................................

    subroutine region_2_supercritical_phase_properties()
      !! Calculate region 2 supercritical phase properties from region
      !! and primary variables for supercritical water, NCG and energy
      !! EOS. Region 2 supercritical fluid is assumed to be
      !! vapour-like.

      ! Locals:
      PetscInt :: p, pseudo_phases
      PetscReal :: water_primary(2), water_properties(2)
      PetscReal :: water_viscosity, water_enthalpy
      PetscReal :: gas_properties(2), xg, pi_liq
      PetscReal, parameter :: density = 0._dp ! not used

      err = 0
      do p = 1, 2
         call fluid%phase(p)%zero()
      end do

      associate (water_pressure => water_primary(1), &
           water_temperature => water_primary(2), &
           region => self%thermo%region(2)%ptr, phase => fluid%phase(3))

        call self%water_pressure(primary, 2, PETSC_FALSE, water_pressure, err)
        if (err == 0) then

           water_temperature = fluid%temperature

           call self%gas%properties(fluid%partial_pressure(2), fluid%temperature, &
                gas_properties, err)

           if (err == 0) then

              call region%properties(water_primary, water_properties, err)
              if (err == 0) then

                 associate (water_density => water_properties(1), &
                      water_internal_energy => water_properties(2), &
                      gas_density => gas_properties(1), gas_enthalpy => gas_properties(2))

                   call self%gas%mass_fraction(fluid%partial_pressure(2), &
                        fluid%temperature, 2, gas_density, water_density, &
                        0._dp, xg, err)

                   if (err == 0) then

                      call region%viscosity(water_temperature, water_pressure, &
                           water_density, water_viscosity)
                      call self%gas%mixture_viscosity(water_viscosity, &
                           fluid%temperature, fluid%partial_pressure(2), xg, 2, &
                           phase%viscosity, err)

                      if (err == 0) then

                         phase%saturation = 1._dp
                         phase%density = water_density + gas_density
                         phase%mass_fraction = [1._dp - xg, xg]
                         phase%relative_permeability = 1._dp
                         phase%capillary_pressure = 0._dp
                         water_enthalpy = water_internal_energy &
                              + water_pressure / water_density
                         phase%specific_enthalpy = water_enthalpy * (1._dp - xg) &
                              + gas_enthalpy * xg
                         phase%internal_energy = phase%specific_enthalpy &
                              - fluid%pressure / phase%density

                         select type (thermo => self%thermo)
                         type is (IAPWS_type)
                            call thermo%pi_liquidlike(water_pressure, &
                                 water_temperature, density, pi_liq, pseudo_phases, err)
                         end select
                         if (err == 0) then
                            fluid%liquidlike_fraction = pi_liq
                            fluid%supercritical_phases = dble(pseudo_phases)
                         end if

                      end if

                   end if
                 end associate
              end if
           end if
        end if
      end associate

    end subroutine region_2_supercritical_phase_properties

  end subroutine eos_sge_region_2_fluid_properties

!------------------------------------------------------------------------

  subroutine eos_sge_effective_gas_properties(self, p, fluid, &
       water_pressure, water_density, gas_density, &
       energy_solution, effective_gas_density, gas_mass_fraction, &
       mixture_viscosity, err)
    !! Calculate effective gas and mixture properties for phase p. For
    !! sub-critical liquid these are interpolated between liquid and
    !! vapour values so that liquid and vapour phase properties are
    !! equal at the critical point.

    use fluid_module, only: fluid_type

    class(eos_sge_type), intent(in out) :: self
    PetscInt, intent(in) :: p !! phase index
    type(fluid_type), intent(in) :: fluid !! Fluid object
    PetscReal, intent(in) :: water_pressure, water_density, gas_density
    PetscReal, intent(out) :: energy_solution, effective_gas_density, &
         gas_mass_fraction, mixture_viscosity
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: xi, xg(2), visc(2), water_viscosity
    PetscReal :: henrys_constant, constituent_henrys_constant(self%gas%num_constituents)
    PetscInt :: pp

    err = 0

    call self%thermo%region(3)%ptr%viscosity(fluid%temperature, &
         water_pressure, water_density, water_viscosity)

    if ((p == 1) .and. (fluid%temperature < self%thermo%critical%temperature)) then

       xi = self%partial_pressure_coefficient(fluid%temperature, &
            liquid = PETSC_TRUE)

       call self%gas%henrys_constant(fluid%temperature, &
            henrys_constant, constituent_henrys_constant, err)
       if (err == 0) then

          effective_gas_density = xi * gas_density

          call self%gas%energy_solution(fluid%temperature, &
               constituent_henrys_constant, energy_solution, err)
          energy_solution = (1._dp - xi) * energy_solution

          do pp = 1, 2
             call self%gas%mass_fraction(fluid%partial_pressure(2), &
                  fluid%temperature, pp, effective_gas_density, &
                  water_density, henrys_constant, xg(pp), err)
             if (err > 0) exit
          end do
          if (err == 0) then

             gas_mass_fraction = (1._dp - xi) * xg(1) + xi * xg(2)

             do pp = 1, 2
                call self%gas%mixture_viscosity(water_viscosity, &
                     fluid%temperature, fluid%partial_pressure(2), &
                     gas_mass_fraction, pp, visc(pp), err)
                if (err > 0) exit
             end do
             if (err == 0) then
                mixture_viscosity = (1._dp - xi) * visc(1) + xi * visc(2)
             end if

          end if
       end if

    else
       effective_gas_density = gas_density
       henrys_constant = 0._dp
       energy_solution = 0._dp
       call self%gas%mass_fraction(fluid%partial_pressure(2), &
            fluid%temperature, 2, effective_gas_density, water_density, &
            henrys_constant, gas_mass_fraction, err)
       call self%gas%mixture_viscosity(water_viscosity, &
            fluid%temperature, fluid%partial_pressure(2), &
            gas_mass_fraction, 2, mixture_viscosity, err)
    end if

  end subroutine eos_sge_effective_gas_properties

!------------------------------------------------------------------------

  subroutine eos_sge_region_3_fluid_properties(self, primary, rock, fluid, err)
    !! Calculate region 3 fluid properties from region and primary
    !! variables for supercritical water, NCG and energy EOS.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_sge_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: p, pp, phases, pseudo_phases, effective_phases
    PetscReal :: water_properties(2), pi_pseudo_phase(2), xg
    PetscReal :: gas_properties(2)
    PetscReal :: viscosity, energy_solution, water_enthalpy
    PetscReal :: effective_water_pressure, effective_water_density
    PetscReal :: effective_water_internal_energy, effective_gas_density

    err = 0

    associate(water_density => primary(1), temperature => primary(2))

      select type (region => self%thermo%region(3)%ptr)
      type is (IAPWS_region3_type)

         fluid%temperature = temperature
         fluid%permeability_factor = 1._dp

         call region%properties(primary(1:2), water_properties, err)

         if (err == 0) then

            associate (water_pressure => water_properties(1), &
                 water_internal_energy => water_properties(2), &
                 partial_pressure => primary(3), &
                 gas_density => gas_properties(1), &
                 gas_enthalpy => gas_properties(2))

              fluid%pressure = water_pressure + partial_pressure
              fluid%partial_pressure = [water_pressure, partial_pressure]
              call self%phase_composition(fluid, err)

              if (err == 0) then

                 do p = 1, self%num_phases
                    call fluid%phase(p)%zero()
                 end do

                 call self%phase_saturations(primary, fluid)

                 call self%gas%properties(fluid%partial_pressure(2), fluid%temperature, &
                      gas_properties, err)

                 if (err == 0) then

                    phases = nint(fluid%phase_composition)
                    p = self%region3_phase(phases)
                    associate (phase => fluid%phase(p))

                      phase%saturation = 1._dp
                      phase%relative_permeability = 1._dp
                      phase%capillary_pressure =  0._dp

                      select type (thermo => self%thermo)
                      type is (IAPWS_type)
                         call thermo%pi_liquidlike(water_pressure, temperature, &
                              water_density, pi_pseudo_phase(1), pseudo_phases, err)
                      end select

                      if (err == 0) then

                         pi_pseudo_phase(2) = 1._dp - pi_pseudo_phase(1)
                         fluid%liquidlike_fraction = pi_pseudo_phase(1)
                         fluid%supercritical_phases = dble(pseudo_phases)
                         if (fluid%is_supercritical()) then
                            effective_phases = pseudo_phases
                         else
                            effective_phases = phases
                         end if

                         ! Loop over pseudo-phases:
                         do pp = 1, 2
                            if (btest(effective_phases, pp - 1)) then

                               call effective_liquid_properties(pp, fluid%pressure, &
                                    temperature, partial_pressure, water_pressure, &
                                    water_density, water_internal_energy, &
                                    effective_water_pressure, effective_water_density, &
                                    effective_water_internal_energy, err)

                               if (err == 0) then

                                  call self%effective_gas_properties(pp, fluid, &
                                       effective_water_pressure, effective_water_density, &
                                       gas_density, energy_solution, effective_gas_density, &
                                       xg, viscosity, err)

                                  if (err == 0) then
                                     phase%density = phase%density + pi_pseudo_phase(pp) * &
                                          (effective_water_density + effective_gas_density)
                                     phase%mass_fraction = phase%mass_fraction + &
                                          pi_pseudo_phase(pp) * [1._dp - xg, xg]
                                     water_enthalpy = effective_water_internal_energy &
                                          + effective_water_pressure / effective_water_density
                                     phase%specific_enthalpy = phase%specific_enthalpy + &
                                          pi_pseudo_phase(pp) * (water_enthalpy * (1._dp - xg) &
                                          + (gas_enthalpy + energy_solution) * xg)
                                     phase%viscosity = phase%viscosity + &
                                          pi_pseudo_phase(pp) * viscosity
                                  else
                                     exit
                                  end if
                               else
                                  exit
                               end if
                            end if
                         end do
                         if (err == 0) then
                            phase%internal_energy = phase%specific_enthalpy &
                                 - fluid%pressure / phase%density
                         end if
                      end if
                    end associate
                 end if
              end if
            end associate
         end if
      end select
    end associate

  contains

    subroutine effective_liquid_properties(p, pressure, temperature, &
         partial_pressure, water_pressure, water_density, &
         water_internal_energy, effective_water_pressure, &
         effective_water_density, effective_water_internal_energy, err)

      !! Return effective water properties for phase p. For
      !! sub-critical liquid these are interpolated between liquid and
      !! vapour values so that liquid and vapour phase properties are
      !! equal at the critical point.

      PetscInt, intent(in) :: p !! phase index
      PetscReal, intent(in) :: pressure, temperature, partial_pressure
      PetscReal, intent(in) :: water_pressure, water_density, &
           water_internal_energy
      PetscReal, intent(out) :: effective_water_pressure, &
           effective_water_density, effective_water_internal_energy
      PetscErrorCode, intent(out) :: err
      ! Locals:
      PetscReal :: xi, props(2)

      err = 0

      if ((p == 1) .and. (temperature < self%thermo%critical%temperature)) then

         xi = self%partial_pressure_coefficient(temperature, liquid = PETSC_TRUE)
         effective_water_pressure = pressure - xi * partial_pressure

         select type (region => self%thermo%region(3)%ptr)
         type is (IAPWS_region3_type)
            call region%density([effective_water_pressure, temperature], &
                 effective_water_density, err, polish = PETSC_TRUE, &
                 phases = SUBREGION_PHASES_LIQUID)
            if (err == 0) then
               call region%properties([effective_water_density, temperature], &
                    props, err)
               if (err == 0) then
                  effective_water_internal_energy = props(2)
               end if
            end if
         end select

      else
         effective_water_pressure = water_pressure
         effective_water_density = water_density
         effective_water_internal_energy = water_internal_energy
      end if

    end subroutine effective_liquid_properties

  end subroutine eos_sge_region_3_fluid_properties

!------------------------------------------------------------------------

  subroutine eos_sge_region_4_fluid_properties(self, primary, rock, fluid, err)
    !! Calculate region 4 fluid properties from region and primary
    !! variables for supercritical water, NCG and energy EOS.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_sge_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err

    call self%eos_wge%bulk_properties(primary, fluid, err)

    associate(pressure => primary(1))
      select type (thermo => self%thermo)
      type is (IAPWS_type)
         if (pressure <= thermo%saturation_pressure_bdy_1_3) then ! T <= 350:
            call self%eos_wge%phase_properties(primary, rock, fluid, err)
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
      PetscReal :: water_properties(2), water_density, water_enthalpy
      PetscReal :: effective_water_pressure, effective_gas_density
      PetscReal :: energy_solution, sl, xg
      PetscReal :: relative_permeability(2), effective_capillary_pressure
      PetscReal :: gas_properties(2)

      err = 0

      select type (region3 => self%thermo%region(3)%ptr)
      type is (IAPWS_region3_type)

         phases = nint(fluid%phase_composition)
         sl = fluid%phase(1)%saturation
         relative_permeability = rock%relative_permeability%values(sl)

         call self%gas%properties(fluid%partial_pressure(2), fluid%temperature, &
              gas_properties, err)

         if (err == 0) then

            do p = 1, 2
               associate(phase => fluid%phase(p), &
                    water_pressure => fluid%partial_pressure(1), &
                    water_internal_energy => water_properties(2), &
                    gas_density => gas_properties(1), gas_enthalpy => gas_properties(2))

                 if (btest(phases, p - 1)) then

                    call effective_liquid_properties(p, fluid, rock, &
                         effective_water_pressure, effective_capillary_pressure, err)

                    if (err == 0) then

                       call region3%saturation_density([effective_water_pressure, &
                            fluid%temperature], p == 1, water_density, err, &
                            polish = PETSC_TRUE)

                       if (err == 0) then

                          call region3%properties([water_density, fluid%temperature], &
                               water_properties, err)

                          if (err == 0) then

                             call self%effective_gas_properties(p, fluid, &
                                  effective_water_pressure, water_density, &
                                  gas_density, energy_solution, effective_gas_density, &
                                  xg, phase%viscosity, err)

                             if (err == 0) then
                                phase%density = water_density + effective_gas_density
                                phase%mass_fraction = [1._dp - xg, xg]
                                phase%relative_permeability = relative_permeability(p)
                                phase%capillary_pressure = effective_capillary_pressure
                                water_enthalpy = water_internal_energy &
                                     + effective_water_pressure / water_density
                                phase%specific_enthalpy = water_enthalpy * (1._dp - xg) &
                                     + (gas_enthalpy + energy_solution) * xg
                                phase%internal_energy = phase%specific_enthalpy &
                                     - fluid%pressure / phase%density
                             else
                                exit
                             end if
                          else
                             exit
                          end if

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
         end if

      end select

    end subroutine region4_above_bdy_1_3_phase_properties

    subroutine effective_liquid_properties(p, fluid, rock, &
         effective_water_pressure, effective_capillary_pressure, err)

      !! Return effective two-phase water component properties for
      !! phase p. For sub-critical liquid these are interpolated
      !! between liquid and vapour values so that liquid and vapour
      !! phase properties are equal at the critical point.

      PetscInt, intent(in) :: p !! phase index
      type(fluid_type), intent(in) :: fluid !! Fluid object
      type(rock_type), intent(in out) :: rock !! Rock object
      PetscReal, intent(out) :: effective_water_pressure, &
           effective_capillary_pressure
      PetscErrorCode, intent(out) :: err
      ! Locals:
      PetscReal :: sl, cp, xi

      err = 0

      if (p == 1) then
         sl = fluid%phase(1)%saturation
         call self%water_pressure(primary, 4, PETSC_TRUE, &
              effective_water_pressure, err)
         if (err == 0) then
            cp = rock%capillary_pressure%value(sl, fluid%temperature)
            xi = self%partial_pressure_coefficient(fluid%temperature, &
                 liquid = PETSC_TRUE)
            effective_capillary_pressure = (1._dp - xi) * cp
         end if
      else
         effective_water_pressure = fluid%partial_pressure(1)
         effective_capillary_pressure = 0._dp
      end if

    end subroutine effective_liquid_properties

  end subroutine eos_sge_region_4_fluid_properties

!------------------------------------------------------------------------

  subroutine eos_sge_primary_variables(self, fluid, primary)
    !! Determine primary variables from fluid properties for
    !! supercritical water, NCG and energy EOS.

    use fluid_module, only: fluid_type

    class(eos_sge_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(out) :: primary(self%num_primary_variables)

    call self%eos_se_type%primary_variables(fluid, primary)

    primary(3) = fluid%partial_pressure(2)

  end subroutine eos_sge_primary_variables

!------------------------------------------------------------------------

   subroutine eos_sge_check_primary_variables(self, fluid, &
       primary, changed, err)
    !! Check if primary variables are in acceptable bounds, and return error
    !! code accordingly.

    use fluid_module, only: fluid_type

    class(eos_sge_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    PetscBool, intent(out) :: changed
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: region
    PetscReal :: total_pressure, water_pressure, max_partial_pressure
    PetscReal :: props(2), temperature
    PetscReal, parameter :: small = 1.e-6_dp

    changed = PETSC_FALSE
    err = 0

    region = nint(fluid%region)

    if (region == 3) then
       temperature = primary(2)
       call self%thermo%region(region)%ptr%properties(primary(1:2), props, err)
       if (err == 0) water_pressure = props(1)
    else
       call self%water_pressure(primary, region, PETSC_FALSE, water_pressure, err)
       if (err == 0) then
          if (region == 4) then
             call self%thermo%saturation%temperature(water_pressure, &
                  temperature, err)
          else
             temperature = primary(2)
          end if
       end if
    end if

    if (err == 0) then
      if ((water_pressure < 0._dp) .or. (water_pressure > 100.e6_dp)) then
         err = 1
      else

         if (err == 0) then

            associate (partial_pressure => primary(3))
              total_pressure = water_pressure + partial_pressure
              max_partial_pressure = (1._dp - small) * total_pressure
              if (partial_pressure > max_partial_pressure) then
                 partial_pressure = max_partial_pressure
                 changed = PETSC_TRUE
              else if (partial_pressure < 0._dp) then
                 partial_pressure = 0._dp
                 changed = PETSC_TRUE
              end if
            end associate

            if (region == 4) then
               associate (vapour_saturation => primary(2))
                 if ((vapour_saturation < -1._dp) .or. &
                      (vapour_saturation > 2._dp)) then
                    err = 1
                 end if
               end associate
            else
               if ((temperature < 0._dp) .or. (temperature > 800._dp)) then
                  err = 1
               end if
            end if

         end if
      end if
   end if

  end subroutine eos_sge_check_primary_variables

!------------------------------------------------------------------------

  function eos_sge_scale_adaptive(self, primary, region) result(scaled_primary)
    !! Non-dimensionalise eos_sge primary variables by scaling. The
    !! first two variables (pressure or density and temperature or
    !! saturation) are scaled by fixed constants. The third variable,
    !! NCG partial pressure, is scaled adaptively by total pressure in
    !! the cell, except in region 3 (where the first variable is
    !! density rather than total pressure).

    class(eos_type), intent(in) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables)
    PetscInt, intent(in) :: region
    PetscReal :: scaled_primary(self%num_primary_variables)

    if (region == 3) then
       scaled_primary = primary / self%primary_scale(:, region)
    else
       scaled_primary(1:2) = primary(1:2) / self%primary_scale(1:2, region)
       associate(scaled_partial_pressure => scaled_primary(3), &
            pressure => primary(1), partial_pressure => primary(3))
         scaled_partial_pressure = partial_pressure / pressure
       end associate
    end if

  end function eos_sge_scale_adaptive

!------------------------------------------------------------------------

  function eos_sge_unscale_adaptive(self, scaled_primary, region) result(primary)
    !! Re-dimensionalise eos_sge scaled primary variables.

    class(eos_type), intent(in) :: self
    PetscReal, intent(in) :: scaled_primary(self%num_primary_variables)
    PetscInt, intent(in) :: region
    PetscReal :: primary(self%num_primary_variables)

    if (region == 3) then
       primary = scaled_primary * self%primary_scale(:, region)
    else
       primary(1:2) = scaled_primary(1:2) * self%primary_scale(1:2, region)
       associate(scaled_partial_pressure => scaled_primary(3), &
            pressure => primary(1), partial_pressure => primary(3))
         partial_pressure = scaled_partial_pressure * pressure
       end associate
    end if

  end function eos_sge_unscale_adaptive

!------------------------------------------------------------------------

end module eos_sge_module
