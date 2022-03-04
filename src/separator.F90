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

  PetscInt, parameter, public :: num_separator_variables = 1
  PetscInt, parameter, public :: max_separator_variable_name_length = 24
  character(max_separator_variable_name_length), parameter, public :: &
       separator_variable_names(num_separator_variables) = [ &
       "steam_fraction      "]

  type :: separator_stage_type
     !! Separator stage.
     private
     PetscReal :: pressure !! Separator pressure for stage
     PetscReal :: ref_water_enthalpy !! Reference enthalpy of water at separator pressure
     PetscReal :: ref_steam_enthalpy !! Reference enthalpy of steam at separator pressure
     PetscReal, public :: steam_fraction !! Steam fraction for stage
     PetscReal, public :: water_rate !! Output separated water mass flow rate for stage
     PetscReal, public :: water_enthalpy !! Output separated water enthalpy for stage
     PetscReal, public :: steam_rate !! Output separated steam mass flow rate for stage
     PetscReal, public :: steam_enthalpy !! Output separated steam enthalpy for stage
   contains
     private
     procedure, public :: init => separator_stage_init
     procedure, public :: separate => separator_stage_separate
  end type separator_stage_type

  type, public :: separator_type
     !! Separator in source network.
     private
     PetscReal, pointer, public :: steam_fraction !! Steam fraction
     PetscBool, public :: on !! Whether separator is active
     type(separator_stage_type), allocatable :: stage(:) !! Separator stages
     PetscInt, public :: num_stages !! Number of separator stages
   contains
     private
     procedure, public :: init => separator_init
     procedure, public :: assign => separator_assign
     procedure, public :: separate => separator_separate
     procedure, public :: zero => separator_zero
     procedure, public :: destroy => separator_destroy
  end type separator_type

contains

!------------------------------------------------------------------------
! separator_stage_type routines
!------------------------------------------------------------------------

  subroutine separator_stage_init(self, pressure, thermo)
    !! Initialise separator stage.

    class(separator_stage_type), intent(in out) :: self
    PetscReal, intent(in) :: pressure !! Stage separator pressure
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics
    ! Locals:
    PetscReal :: saturation_temperature
    PetscReal :: params(2), water_props(2), steam_props(2)
    PetscErrorCode :: err

    self%pressure = pressure
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

  end subroutine separator_stage_init

!------------------------------------------------------------------------

  subroutine separator_stage_separate(self, rate, enthalpy)
    !! Calculates steam fraction and separated water and steam
    !! properties, from given input mass flow rate and enthalpy.

    class(separator_stage_type), intent(in out) :: self
    PetscReal, intent(in) :: rate !! Input mass flow rate
    PetscReal, intent(in) :: enthalpy !! Input enthalpy

    if (enthalpy <= self%ref_water_enthalpy) then
       self%steam_fraction = 0._dp
       self%water_enthalpy = enthalpy
       self%steam_enthalpy = 0._dp
    else if (enthalpy <= self%ref_steam_enthalpy) then
       self%steam_fraction = (enthalpy - self%ref_water_enthalpy) / &
            (self%ref_steam_enthalpy - self%ref_water_enthalpy)
       self%water_enthalpy = self%ref_water_enthalpy
       self%steam_enthalpy = self%ref_steam_enthalpy
    else
       self%steam_fraction = 1._dp
       self%water_enthalpy = 0._dp
       self%steam_enthalpy = enthalpy
    end if

    self%water_rate = (1._dp - self%steam_fraction) * rate
    self%steam_rate = self%steam_fraction * rate

  end subroutine separator_stage_separate

!------------------------------------------------------------------------
! separator_type routines
!------------------------------------------------------------------------

  subroutine separator_init(self, pressure, thermo)
    !! Initialise separator. (Separator assign() method must be called
    !! before use.)

    class(separator_type), intent(in out) :: self
    PetscReal, intent(in) :: pressure(:) !! Stage separator pressures
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics
    ! Locals:
    PetscInt :: i

    if (all(pressure > 0._dp)) then
       self%on = PETSC_TRUE
       self%num_stages = size(pressure)
       allocate(self%stage(self%num_stages))
       do i = 1, self%num_stages
          call self%stage(i)%init(pressure(i), thermo)
       end do
    else
       self%on = PETSC_FALSE
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

    self%steam_fraction => data(offset)

  end subroutine separator_assign

!------------------------------------------------------------------------

  subroutine separator_separate(self, rate, enthalpy, water_rate, &
       water_enthalpy, steam_rate, steam_enthalpy)
    !! Calculates steam fraction and separated water and steam
    !! properties, from given input mass flow rate and enthalpy.

    class(separator_type), intent(in out) :: self
    PetscReal, intent(in) :: rate !! Input mass flow rate
    PetscReal, intent(in) :: enthalpy !! Input enthalpy
    PetscReal, intent(out) :: water_rate !! Output water mass flow rate
    PetscReal, intent(out) :: water_enthalpy !! Output water enthalpy
    PetscReal, intent(out) :: steam_rate !! Output steam mass flow rate
    PetscReal, intent(out) :: steam_enthalpy !! Output steam enthalpy
    ! Locals:
    PetscInt :: i
    PetscReal :: total_steam_mass_rate, total_steam_energy_rate
    PetscReal :: q, h
    PetscReal, parameter :: tol = 1.e-9_dp

    q = rate
    h = enthalpy
    total_steam_mass_rate = 0._dp
    total_steam_energy_rate = 0._dp

    do i = 1, self%num_stages
       associate(stage => self%stage(i))
         call stage%separate(q, h)
         total_steam_mass_rate = total_steam_mass_rate + stage%steam_rate
         total_steam_energy_rate = total_steam_energy_rate + stage%steam_rate * &
              stage%steam_enthalpy
         q = stage%water_rate
         h = stage%water_enthalpy
       end associate
    end do

    water_rate = q
    water_enthalpy = h
    steam_rate = total_steam_mass_rate
    if (abs(total_steam_mass_rate) > tol) then
       steam_enthalpy = total_steam_energy_rate / total_steam_mass_rate
    else
       steam_enthalpy = 0._dp
    end if
    if (abs(rate) > tol) then
       self%steam_fraction = steam_rate / rate
    else
       self%steam_fraction = 0._dp
    end if

  end subroutine separator_separate

!------------------------------------------------------------------------

  subroutine separator_zero(self)
    !! Zeroes output quantities (e.g. when separator is not on).

    class(separator_type), intent(in out) :: self

    self%steam_fraction = 0._dp

  end subroutine separator_zero

!------------------------------------------------------------------------

  subroutine separator_destroy(self)
    !! Destroys separator object.

    class(separator_type), intent(in out) :: self

    self%steam_fraction => null()
    deallocate(self%stage)

  end subroutine separator_destroy

!------------------------------------------------------------------------

end module separator_module
