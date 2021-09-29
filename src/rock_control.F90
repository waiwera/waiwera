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

module rock_control_module
  !! Module for rock controls- for controlling rock parameters (e.g. permeability) over time.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use list_module
  use rock_module
  use interpolation_module

  implicit none
  private

  type, public, abstract :: rock_control_type
     !! Abstract type for rock control, controlling rock parameters
     !! over time, in one or more cells.
     private
   contains
     procedure(rock_control_destroy_procedure), public, deferred :: destroy
     procedure(rock_control_update_procedure), public, deferred :: update
  end type rock_control_type

  abstract interface

     subroutine rock_control_destroy_procedure(self)
       !! Destroys a rock control object.
       import :: rock_control_type
       class(rock_control_type), intent(in out) :: self
     end subroutine rock_control_destroy_procedure

     subroutine rock_control_update_procedure(self, t, &
          rock_data, rock_section, rock_range_start)
       use petscis
       !! Updates rock properties at the specified time.
       import :: rock_control_type
       class(rock_control_type), intent(in out) :: self
       PetscReal, intent(in) :: t
       PetscReal, pointer, contiguous, intent(in) :: rock_data(:)
       PetscSection, intent(in) :: rock_section
       PetscInt, intent(in) :: rock_range_start
     end subroutine rock_control_update_procedure
     
  end interface

end module rock_control_module
