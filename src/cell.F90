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
     procedure, public :: tracer_balance_coefs => cell_tracer_balance_coefs
     procedure, public :: diffusion_factor => cell_diffusion_factor
     procedure, public :: tortuosity => cell_tortuosity
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

  function cell_tracer_balance_coefs(self, tracer_phase_indices) result(coefs)
    !! Returns tracer balance coefficients for the cell. These are the
    !! coefficients which when multiplied by the tracer mass fractions
    !! give the tracer mass balances (per unit volume) in the cell.

    class(cell_type), intent(in) :: self
    PetscInt, intent(in) :: tracer_phase_indices(:) !! phase index for each tracer
    PetscReal :: coefs(size(tracer_phase_indices))
    ! Locals:
    PetscInt :: i

    coefs = 0._dp
    do i = 1, size(tracer_phase_indices)
       associate(phase => self%fluid%phase(tracer_phase_indices(i)))
         coefs(i) = self%rock%porosity * phase%saturation * phase%density
       end associate
    end do

  end function cell_tracer_balance_coefs

!------------------------------------------------------------------------

  PetscReal function cell_diffusion_factor(self, p) result(factor)
    !! Returns the diffusion factor for phase p, which multiplies the
    !! diffusion coefficient in tracer diffusion. The diffusion factor
    !! is given by the product of three quantities: rock porosity,
    !! fluid phase density (for phase p) and tortuosity.

    class(cell_type), intent(in) :: self
    PetscInt, intent(in) :: p !! phase index

    factor = self%rock%porosity * self%fluid%phase(p)%density * &
              self%tortuosity(p)

  end function cell_diffusion_factor

!------------------------------------------------------------------------

  PetscReal function cell_tortuosity(self, p) result(tortuosity)
    !! Returns the effective tortuosity in the cell, for the specified
    !! phase. In general, this is the product of two factors, a rock
    !! tortuosity and a fluid tortuosity. Here, the rock tortuosity is
    !! taken to be identically 1 and the fluid tortuosity is taken to
    !! be equal to the phase saturation, giving a constant-diffusivity
    !! formulation.

    class(cell_type), intent(in) :: self
    PetscInt, intent(in) :: p !! phase index
    ! Locals:
    PetscReal :: fluid_tortuosity
    PetscReal, parameter :: rock_tortuosity = 1._dp

    fluid_tortuosity = self%fluid%phase(p)%saturation
    tortuosity = rock_tortuosity * fluid_tortuosity

  end function cell_tortuosity

!------------------------------------------------------------------------

end module cell_module
