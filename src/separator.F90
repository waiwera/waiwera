!   Copyright 2021 University of Auckland.

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

module separator_module
  !! Module for separators in source network, which separate two-phase
  !! fluid into water and steam.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use thermodynamics_module

  implicit none
  private

  PetscReal, parameter, public :: default_separator_pressure = 0.55e6_dp

  PetscInt, parameter, public :: num_separator_variables = 6
  PetscInt, parameter, public :: max_separator_variable_name_length = 24
  character(max_separator_variable_name_length), parameter, public :: &
       separator_variable_names(num_separator_variables) = [ &
       "separator_pressure  ", "steam_fraction      ", &
       "water_rate          ", "water_enthalpy      ", &
       "steam_rate          ", "steam_enthalpy      "]

  type, public :: separator_type
     !! Separator in source network.
     private
     PetscReal :: ref_water_enthalpy !! Reference enthalpy of water at separator pressure
     PetscReal :: ref_steam_enthalpy !! Reference enthalpy of steam at separator pressure
     PetscReal, pointer, public :: pressure !! Separator pressure
     PetscReal, pointer, public :: steam_fraction !! Steam fraction
     PetscReal, pointer, public :: water_rate !! Output separated water mass flow rate
     PetscReal, pointer, public :: water_enthalpy !! Output separated water enthalpy
     PetscReal, pointer, public :: steam_rate !! Output separated steam mass flow rate
     PetscReal, pointer, public :: steam_enthalpy !! Output separated steam enthalpy
     PetscBool, public :: on !! Whether separator is active
   contains
     private
     procedure, public :: init => separator_init
     procedure, public :: assign => separator_assign
     procedure, public :: separate => separator_separate
     procedure, public :: zero => separator_zero
     procedure, public :: destroy => separator_destroy
     procedure :: get_steam_fraction => separator_get_steam_fraction
     procedure :: get_separated_rates => separator_get_separated_rates
     procedure :: get_separated_enthalpies => separator_get_separated_enthalpies
  end type separator_type

contains

!------------------------------------------------------------------------

  subroutine separator_init(self, pressure, thermo)
    !! Initialise separator. (Separator assign() method must be called
    !! before use.)

    class(separator_type), intent(in out) :: self
    PetscReal, intent(in) :: pressure !! Separator pressure
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics
    ! Locals:
    PetscReal :: saturation_temperature
    PetscReal :: params(2), water_props(2), steam_props(2)
    PetscErrorCode :: err

    self%pressure = pressure
    self%on = (self%pressure > 0._dp)

    if (self%on) then

       call thermo%saturation%temperature(self%pressure, &
            saturation_temperature, err)
       params = [self%pressure, saturation_temperature]
       call thermo%water%properties(params, water_props, err)
       call thermo%steam%properties(params, steam_props, err)

       associate(water_density => water_props(1), &
            water_internal_energy => water_props(2), &
            steam_density => steam_props(1), &
            steam_internal_energy => steam_props(2))
         self%ref_water_enthalpy = water_internal_energy + &
              self%pressure / water_density
         self%ref_steam_enthalpy = steam_internal_energy + &
              self%pressure / steam_density
       end associate

    else
       call self%zero()
    end if

  end subroutine separator_init
  
!------------------------------------------------------------------------

  subroutine separator_assign(self, data, offset)
    !! Assigns pointers in separator object to elements in the data
    !! array, starting from the specified offset.

    class(separator_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! source data array
    PetscInt, intent(in) :: offset  !! source array offset

    self%pressure => data(offset)
    self%steam_fraction => data(offset + 1)
    self%water_rate => data(offset + 2)
    self%water_enthalpy => data(offset + 3)
    self%steam_rate => data(offset + 4)
    self%steam_enthalpy => data(offset + 5)

  end subroutine separator_assign

!------------------------------------------------------------------------

  subroutine separator_separate(self, rate, enthalpy)
    !! Calculates steam fraction and separated water and steam
    !! properties, from given input mass flow rate and enthalpy.

    class(separator_type), intent(in out) :: self
    PetscReal, intent(in) :: rate !! Input mass flow rate
    PetscReal, intent(in) :: enthalpy !! Input enthalpy

    call self%get_steam_fraction(enthalpy)
    call self%get_separated_rates(rate)
    call self%get_separated_enthalpies(enthalpy)
    
  end subroutine separator_separate

!------------------------------------------------------------------------

  subroutine separator_zero(self)
    !! Zeroes output quantities (e.g. when separator is not on).

    class(separator_type), intent(in out) :: self

    self%steam_fraction = 0._dp
    self%water_rate = 0._dp
    self%water_enthalpy = 0._dp
    self%steam_rate = 0._dp
    self%steam_enthalpy = 0._dp

  end subroutine separator_zero

!------------------------------------------------------------------------

  subroutine separator_get_steam_fraction(self, enthalpy)
    !! Calculates steam fraction for specified input enthalpy.

    class(separator_type), intent(in out) :: self
    PetscReal, intent(in) :: enthalpy !! Input enthalpy

    if (enthalpy <= self%ref_water_enthalpy) then
       self%steam_fraction = 0._dp
    else if (enthalpy <= self%ref_steam_enthalpy) then
       self%steam_fraction = (enthalpy - self%ref_water_enthalpy) / &
            (self%ref_steam_enthalpy - self%ref_water_enthalpy)
    else
       self%steam_fraction = 1._dp
    end if

  end subroutine separator_get_steam_fraction

!------------------------------------------------------------------------

  subroutine separator_get_separated_rates(self, rate)
    !! Calculates separated water and steam mass flow rates from the
    !! given input mass flow rate.

    class(separator_type), intent(in out) :: self
    PetscReal, intent(in) :: rate

    self%water_rate = (1._dp - self%steam_fraction) * rate
    self%steam_rate = self%steam_fraction * rate

  end subroutine separator_get_separated_rates

!------------------------------------------------------------------------

  subroutine separator_get_separated_enthalpies(self, enthalpy)
    !! Calculates separated water and steam enthalpies from the given
    !! input enthalpy.

    class(separator_type), intent(in out) :: self
    PetscReal, intent(in) :: enthalpy

    if (enthalpy <= self%ref_water_enthalpy) then
       self%water_enthalpy = enthalpy
       self%steam_enthalpy = 0._dp
    else if (enthalpy <= self%ref_steam_enthalpy) then
       self%water_enthalpy = self%ref_water_enthalpy
       self%steam_enthalpy = self%ref_steam_enthalpy
    else
       self%water_enthalpy = 0._dp
       self%steam_enthalpy = enthalpy
    end if

  end subroutine separator_get_separated_enthalpies

!------------------------------------------------------------------------

  subroutine separator_destroy(self)
    !! Destroys separator object.

    class(separator_type), intent(in out) :: self

    self%pressure => null()
    self%steam_fraction => null()
    self%water_rate => null()
    self%water_enthalpy => null()
    self%steam_rate => null()
    self%steam_enthalpy => null()

  end subroutine separator_destroy

!------------------------------------------------------------------------

end module separator_module
