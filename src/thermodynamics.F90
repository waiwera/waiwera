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

module thermodynamics_module
  !! Thermodynamics constants and abstract types.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module

  implicit none
  private

  PetscInt, parameter, public :: max_thermodynamics_name_length = 16
  PetscInt, parameter, public :: max_thermodynamic_region_name_length = 16

!------------------------------------------------------------------------
! Physical constants
!------------------------------------------------------------------------

  PetscReal, parameter, public :: tc_k = 273.15_dp         !! Conversion from Celsius to Kelvin
  PetscReal, parameter, public :: water_molecular_weight = 18.01528_dp !! Molecular weight of water (g/mol)
  PetscReal, parameter, public :: gas_constant = 8.3144598_dp    !! Gas constant R
  PetscReal, parameter, public :: specific_gas_constant = 0.461526e3_dp !! Specific gas constant R for water
  PetscReal, parameter, public :: ttriple = 0.01_dp              !! Triple point of water

!------------------------------------------------------------------------
! Critical point type
!------------------------------------------------------------------------

  type, public :: critical_point_type
     !! Type for critical point data.
     private
     PetscReal, public :: temperature_k, temperature, pressure, density
  end type critical_point_type

!------------------------------------------------------------------------
! Saturation curve type
!------------------------------------------------------------------------

  type, public, abstract :: saturation_type
     !! Saturation curve type.
     class(thermodynamics_type), pointer, public :: thermo
   contains
     private
       procedure, public :: init => saturation_init
       procedure(saturation_temperature), public, deferred :: temperature
       procedure(saturation_pressure), public, deferred :: pressure
  end type saturation_type

!------------------------------------------------------------------------
! Thermodynamic region type
!------------------------------------------------------------------------

  type, public, abstract :: region_type
     !! Thermodynamic region type.
     character(max_thermodynamic_region_name_length), public :: name
     class(thermodynamics_type), pointer, public :: thermo
   contains
     private
     procedure(region_init), public, deferred :: init
     procedure(region_destroy), public, deferred :: destroy
     procedure(region_properties), public, deferred :: properties
     procedure(region_viscosity), public, deferred :: viscosity
  end type region_type

  ! Pointer to region:
  type, public :: pregion_type
     !! Pointer to thermodynamic region.
     class(region_type), pointer, public :: ptr
   contains
     procedure, public :: set => pregion_set
  end type pregion_type

!------------------------------------------------------------------------
! Thermodynamics type
!------------------------------------------------------------------------

  type, public, abstract :: thermodynamics_type
     !! Thermodynamics type.
     private
     character(max_thermodynamics_name_length), public :: name
     type(critical_point_type), public :: critical !! Critical point
     class(saturation_type), allocatable, public :: saturation !! Saturation curve
     PetscInt, public :: num_regions  !! Number of thermodynamic regions
     class(region_type), allocatable, public :: water !! Pure water region
     class(region_type), allocatable, public :: steam !! Steam region
     class(region_type), allocatable, public :: supercritical !! Supercritical region
     type(pregion_type), allocatable, public :: region(:) !! Array of region pointers
   contains
     private
     procedure(thermodynamics_init_procedure), public, deferred :: init
     procedure(thermodynamics_destroy_procedure), public, deferred :: destroy
     procedure(thermodynamics_phase_composition_procedure), public, deferred :: &
          phase_composition
  end type thermodynamics_type

!------------------------------------------------------------------------

  abstract interface

     subroutine saturation_temperature(self, p, t, err)
       !! Calculates saturation temperature as a function of pressure.
       import :: saturation_type, dp
       class(saturation_type), intent(in) :: self
       PetscReal, intent(in) :: p  !! Fluid pressure (\(kg. m. s^{-1}\))
       PetscReal, intent(out):: t  !! Fluid temperature (\(^\circ C\))
       PetscInt, intent(out) :: err !! Error code
     end subroutine saturation_temperature

     subroutine saturation_pressure(self, t, p, err)
       !! Calculates saturation pressure as a function of temperature.
       import :: saturation_type, dp
       class(saturation_type), intent(in) :: self
       PetscReal, intent(in) :: t  !! Fluid temperature (\(^\circ C\))
       PetscReal, intent(out):: p  !! Fluid pressure (\(kg. m. s^{-1}\))
       PetscInt, intent(out) :: err  !! Error code
     end subroutine saturation_pressure

     subroutine region_init(self, thermo, extrapolate)
       !! Initializes region. The extrapolate parameter allows the
       !! region's methods to be called slightly out of their usual
       !! operating range if needed.
       import :: region_type, thermodynamics_type
       class(region_type), intent(in out) :: self
       class(thermodynamics_type), intent(in), target :: thermo
       PetscBool, intent(in), optional :: extrapolate
     end subroutine region_init

     subroutine region_destroy(self)
       !! Destroys region.
       import :: region_type
       class(region_type), intent(in out) :: self
     end subroutine region_destroy

     subroutine region_properties(self, param, props, err)
       !! Calculates fluid properties for region.
       import :: region_type, dp
       class(region_type), intent(in out) :: self
       PetscReal, intent(in) :: param(:)
       PetscReal, intent(out) :: props(:)
       PetscInt, intent(out) :: err
     end subroutine region_properties

     subroutine region_viscosity(self, temperature, pressure, density, viscosity)
       !! Calculates viscosity for region.
       import :: region_type, dp
       class(region_type), intent(in out) :: self
       PetscReal, intent(in) :: temperature, pressure, density
       PetscReal, intent(out) :: viscosity
     end subroutine region_viscosity

     subroutine thermodynamics_init_procedure(self, extrapolate)
       !! Initializes thermodynamics.
       import :: thermodynamics_type
       class(thermodynamics_type), intent(in out) :: self
       PetscBool, intent(in), optional :: extrapolate
     end subroutine thermodynamics_init_procedure

     subroutine thermodynamics_destroy_procedure(self)
       !! Destroys thermodynamics.
       import :: thermodynamics_type
       class(thermodynamics_type), intent(in out) :: self
     end subroutine thermodynamics_destroy_procedure

     PetscInt function thermodynamics_phase_composition_procedure(self, &
          region, pressure, temperature)
       !! Returns phase composition for given region, pressure and
       !! temperature. This is an integer with each binary bit
       !! representing the presence or absence of a particular phase.
       !! A negative value is returned if the routine fails.
       import :: thermodynamics_type
       class(thermodynamics_type), intent(in) :: self
       PetscInt, intent(in) :: region
       PetscReal, intent(in) :: pressure, temperature
     end function thermodynamics_phase_composition_procedure

  end interface

!------------------------------------------------------------------------

contains

!------------------------------------------------------------------------
! Saturation type
!------------------------------------------------------------------------

  subroutine saturation_init(self, thermo)
    !! Initialise saturation type.
    class(saturation_type), intent(in out) :: self
    class(thermodynamics_type), intent(in), target :: thermo

    self%thermo => thermo

  end subroutine saturation_init

!------------------------------------------------------------------------
! Region pointers
!------------------------------------------------------------------------

  subroutine pregion_set(self, tgt)
    !! Sets a region pointer. This is just a workaround to give
    !! tgt the 'target' attribute, which can't always be done as part of
    !! its declaration, e.g. if it's a component of a derived type.

    class(pregion_type), intent(in out) :: self
    class(region_type), target, intent(in) :: tgt

    self%ptr => tgt

  end subroutine pregion_set

!------------------------------------------------------------------------

end module thermodynamics_module
