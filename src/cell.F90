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

module cell_module
  !! Defines types for accessing local quantities defined on a cell- geometry, rock
  !! and fluid properties.
  !! The components of these types all point to values in arrays obtained from
  !! parallel vectors.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use rock_module
  use fluid_module

  implicit none
  private

  type cell_type
     !! Type for accessing local cell geometry, rock and
     !! fluid properties.
     private
     PetscReal, pointer, public :: volume       !! cell volume
     PetscReal, pointer, contiguous, public :: centroid(:)  !! cell centroid
     type(rock_type), public :: rock   !! rock properties
     type(fluid_type), public :: fluid !! fluid properties
     PetscInt, public :: dof !! Number of degrees of freedom
   contains
     private
     procedure, public :: init => cell_init
     procedure, public :: assign_geometry => cell_assign_geometry
     procedure, public :: destroy => cell_destroy
     procedure, public :: balance => cell_balance
  end type cell_type

  PetscInt, parameter, public :: num_cell_variables = 2
  PetscInt, parameter, public :: &
       cell_variable_num_components(num_cell_variables) = &
       [3, 1] !! Number of components in each cell variable
  PetscInt, parameter, public :: max_cell_variable_name_length = 24
  character(max_cell_variable_name_length), parameter, public :: &
       cell_variable_names(num_cell_variables) = &
       [character(max_cell_variable_name_length):: "centroid", "volume"]

  public :: cell_type

contains

!------------------------------------------------------------------------

  subroutine cell_init(self, num_components, num_phases)
    !! Initialises a cell.

    class(cell_type), intent(in out) :: self
    PetscInt, intent(in) :: num_components !! Number of fluid components
    PetscInt, intent(in) :: num_phases  !! Number of fluid phases

    call self%fluid%init(num_components, num_phases)
    call self%rock%init()

    self%dof = sum(cell_variable_num_components)

  end subroutine cell_init

!------------------------------------------------------------------------

  subroutine cell_assign_geometry(self, data, offset)
    !! Assigns cell geometry pointers to values from specified data
    !! array, starting from the given offset.

    class(cell_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! array with geometry data
    PetscInt, intent(in)  :: offset  !! geometry array offset for this cell

    self%centroid => data(offset: offset + 2)
    self%volume => data(offset + 3)

  end subroutine cell_assign_geometry

!------------------------------------------------------------------------

  subroutine cell_destroy(self)
    !! Destroys a cell.

    class(cell_type), intent(in out) :: self

    nullify(self%volume)
    nullify(self%centroid)
    call self%fluid%destroy()
    call self%rock%destroy()
    
  end subroutine cell_destroy

!------------------------------------------------------------------------

  function cell_balance(self, num_primary) result(balance)
    !! Returns array containing mass balance (per unit volume) for each
    !! mass component in the cell, and energy balance for non-isothermal
    !! simulations.

    class(cell_type), intent(in) :: self
    PetscInt, intent(in) :: num_primary
    PetscReal :: balance(num_primary)
    ! Locals:
    PetscInt :: nc
    PetscReal :: er, ef
    PetscBool :: isothermal

    nc = self%fluid%num_components
    isothermal = (num_primary == nc)

    ! Mass balances:
    balance(1: nc) = self%rock%porosity * &
         self%fluid%component_density()

    if (.not. isothermal) then
       ! Energy balance:
       er = self%rock%energy(self%fluid%temperature)
       ef = self%fluid%energy()
       balance(num_primary) = self%rock%porosity * ef + &
            (1._dp - self%rock%porosity) * er
    end if

  end function cell_balance

!------------------------------------------------------------------------
  
end module cell_module
