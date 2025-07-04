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
     procedure, public :: partial_pressures => eos_sge_partial_pressures
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
    procedure(root_finder_routine), pointer :: fs, fw, ft
    PetscReal :: pressure_scale, temperature_scale, density_scale, partial_pressure_scale
    PetscInt :: scale_type, modifier_type
    character(max_fluid_modifier_name_length) :: relative_permeability_modifier_type_name
    type(fson_value), pointer :: rperm_json
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
       partial_pressure_scale = 0._dp
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
    call init_line_finder(self%saturation_line_finder, &
         self%primary_variable_interpolator, fs, init_interpolator = PETSC_TRUE)
    fw => eos_sge_widom_delta_difference
    allocate(widom_delta_interpolator_type :: self%widom_delta_interpolator)
    call init_line_finder(self%widom_delta_finder, &
         self%widom_delta_interpolator, fw, init_interpolator = PETSC_TRUE)

    call fson_get_mpi(json, "eos.conditions", default_conditions, &
         conditions, logfile)
    self%pressure_conditions = (str_to_lower(conditions) == "pressure")

    ! Set up relative permeability modifier:
    modifier_type = fson_type_mpi(json, "eos.relative_permeability_modifier")
    select case (modifier_type)
    case (TYPE_OBJECT)
       call fson_get_mpi(json, "eos.relative_permeability_modifier.type", &
            default_relative_permeability_modifier_type_name, &
            relative_permeability_modifier_type_name, logfile)
       select case (str_to_lower(relative_permeability_modifier_type_name))
       case ("linear")
          allocate(fluid_relative_permeability_linear_temperature_type :: &
               self%relative_permeability_modifier)
          select type (modifier => self%relative_permeability_modifier)
          type is (fluid_relative_permeability_linear_temperature_type)
             modifier%critical_temperature = self%thermo%critical%temperature
          end select
       case default ! null modifier
          allocate(fluid_modifier_type :: self%relative_permeability_modifier)
       end select
    case (TYPE_NULL)
       allocate(fluid_modifier_type :: self%relative_permeability_modifier)
    end select
    if (fson_has_mpi(json, "eos.relative_permeability_modifier")) then
       call fson_get_mpi(json, "eos.relative_permeability_modifier", rperm_json)
    else
       rperm_json => null()
    end if
    call self%relative_permeability_modifier%init(rperm_json, logfile)

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

  PetscReal function eos_sge_water_pressure(self, primary) result(water_pressure)
    !! For eos_sge, return water pressure from primary variables
    !! (regions 1, 2, 4).

    class(eos_sge_type), intent(in) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables)

    water_pressure = self%eos_wge%water_pressure(primary)

  end function eos_sge_water_pressure

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

  function eos_sge_partial_pressures(self, primary) result (partial_pressures)
    !! Set partial pressures from primary variables.

    class(eos_sge_type), intent(in) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    PetscReal :: partial_pressures(self%num_components)

    partial_pressures = self%eos_wge%partial_pressures(primary)

  end function eos_sge_partial_pressures

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

      associate (water_pressure => water_primary(1), water_temperature => water_primary(2), &
           region => self%thermo%region(2)%ptr, phase => fluid%phase(3))

        water_pressure = fluid%partial_pressure(1)
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
      end associate

    end subroutine region_2_supercritical_phase_properties

  end subroutine eos_sge_region_2_fluid_properties

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
    PetscReal :: henrys_constant, constituent_henrys_constant(self%gas%num_constituents)
    PetscReal :: gas_properties(2), effective_gas_properties(2)
    PetscReal :: viscosity, energy_solution, water_enthalpy, water_viscosity

    err = 0

    associate(water_density => primary(1), temperature => primary(2))

      select type (region => self%thermo%region(3)%ptr)
      type is (IAPWS_region3_type)

         fluid%temperature = temperature
         call region%properties(primary, water_properties, err)

         if (err == 0) then

            associate (water_pressure => water_properties(1), &
                 water_internal_energy => water_properties(2), &
                 partial_pressure => primary(3), &
                 gas_density => effective_gas_properties(1), &
                 gas_enthalpy => effective_gas_properties(2))

              fluid%pressure = water_pressure + partial_pressure
              fluid%partial_pressure = [water_pressure, partial_pressure]
              fluid%permeability_factor = 1._dp

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
                         call thermo%pi_liquidlike(water_pressure, temperature, water_density, &
                              pi_pseudo_phase(1), pseudo_phases, err)
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
                               if (pp == 1) then
                                  water_pressure = fluid%pressure
                                  call self%gas%henrys_constant(fluid%temperature, &
                                       henrys_constant, constituent_henrys_constant, err)
                                  if (err == 0) then
                                     call self%gas%energy_solution(fluid%temperature, &
                                          constituent_henrys_constant, energy_solution, err)
                                  end if
                               else
                                  water_pressure = fluid%partial_pressure(1)
                                  henrys_constant = 0._dp
                                  energy_solution = 0._dp
                               end if
                               if (err == 0) then
                                  call self%gas%effective_properties(gas_properties, pp, &
                                       effective_gas_properties)
                                  call self%gas%mass_fraction(fluid%partial_pressure(2), &
                                       fluid%temperature, pp, gas_density, water_density, &
                                       henrys_constant, xg, err)
                                  if (err == 0) then
                                     call region%viscosity(fluid%temperature, fluid%pressure, &
                                          water_density, water_viscosity)
                                     call self%gas%mixture_viscosity(water_viscosity, &
                                          fluid%temperature, fluid%partial_pressure(2), xg, pp, &
                                          viscosity, err)
                                     if (err == 0) then
                                        phase%density = phase%density + pi_pseudo_phase(pp) * &
                                             (water_density + gas_density)
                                        phase%mass_fraction = phase%mass_fraction + &
                                             pi_pseudo_phase(pp) * [1._dp - xg, xg]
                                        water_enthalpy = water_internal_energy &
                                             + water_pressure / water_density
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
      PetscReal :: water_density, water_properties(2), water_pressure, water_enthalpy
      PetscReal :: water_viscosity, energy_solution, sl, xg
      PetscReal :: relative_permeability(2), capillary_pressure
      PetscBool :: liquid
      PetscReal :: gas_properties(2), effective_gas_properties(2)
      PetscReal :: henrys_constant, constituent_henrys_constant(self%gas%num_constituents)

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
               associate(phase => fluid%phase(p))

                 if (btest(phases, p - 1)) then

                    liquid = (p == 1)

                    if (liquid) then
                       water_pressure = fluid%pressure
                       capillary_pressure = rock%capillary_pressure%value(sl, &
                            fluid%temperature)
                       call self%gas%henrys_constant(fluid%temperature, henrys_constant, &
                            constituent_henrys_constant, err)
                       if (err == 0) then
                          call self%gas%energy_solution(fluid%temperature, &
                               constituent_henrys_constant, energy_solution, err)
                       end if
                    else
                       water_pressure = fluid%partial_pressure(1)
                       capillary_pressure = 0._dp
                       henrys_constant = 0._dp
                       energy_solution = 0._dp
                    end if

                    call region3%saturation_density([water_pressure, &
                         fluid%temperature], liquid, water_density, err, &
                         polish = PETSC_TRUE)

                    if (err == 0) then

                       call region3%properties([water_density, fluid%temperature], &
                            water_properties, err)

                       if (err == 0) then

                          call self%gas%effective_properties(gas_properties, p, &
                               effective_gas_properties)

                          associate(water_internal_energy => water_properties(2), &
                               gas_density => effective_gas_properties(1), &
                               gas_enthalpy => effective_gas_properties(2))

                            call self%gas%mass_fraction(fluid%partial_pressure(2), &
                                 fluid%temperature, p, gas_density, water_density, &
                                 henrys_constant, xg, err)

                            if (err == 0) then

                               call region3%viscosity(fluid%temperature, water_pressure, &
                                    water_density, water_viscosity)
                               call self%gas%mixture_viscosity(water_viscosity, &
                                    fluid%temperature, fluid%partial_pressure(2), xg, p, &
                                    phase%viscosity, err)

                               if (err == 0) then
                                  phase%density = water_density + gas_density
                                  phase%mass_fraction = [1._dp - xg, xg]
                                  phase%relative_permeability = relative_permeability(p)
                                  phase%capillary_pressure = capillary_pressure
                                  water_enthalpy = water_internal_energy &
                                       + water_pressure / water_density
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
                          end associate

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
    PetscReal :: total_pressure, water_pressure, max_partial_pressure, props(2)
    PetscReal, parameter :: small = 1.e-6_dp

    changed = PETSC_FALSE
    err = 0

    region = nint(fluid%region)
    if (region == 3) then
       call self%thermo%region(region)%ptr%properties(primary(1:2), props, err)
       if (err == 0) water_pressure = props(1)
    else
       water_pressure = self%water_pressure(primary)
    end if

    if (err == 0) then
      if ((water_pressure < 0._dp) .or. (water_pressure > 100.e6_dp)) then
         err = 1
      else

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
            associate (t => primary(2))
              if ((t < 0._dp) .or. (t > 800._dp)) then
                 err = 1
              end if
            end associate
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

  subroutine eos_sge_widom_delta_difference(x, context, f, err)
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
    PetscReal, allocatable :: var(:), Pw
    PetscReal :: delta(2)

    err = 0
    select type (context)
    type is (widom_delta_interpolator_type)
       allocate(var(context%dim))
       var = context%interpolate_at_index(x)
       associate(P => var(1), T => var(2), Pg => var(3))
         select type (thermo => context%thermo)
         type is (IAPWS_type)
            Pw = P - Pg
            call thermo%widom_delta(Pw, delta, err)
            if (err == 0) then
               f = T - delta(context%bdy_index)
            end if
         end select
       end associate
       deallocate(var)
    end select

  end subroutine eos_sge_widom_delta_difference

!------------------------------------------------------------------------

end module eos_sge_module
