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
  use control_module

  implicit none
  private

  type, public, extends(table_vector_control_type) :: permeability_table_rock_control_type
     !! Controls rock permeability via a table of values vs. time.
   contains
     procedure, public :: update => permeability_table_rock_control_update
  end type permeability_table_rock_control_type

  type, public, extends(table_vector_control_type) :: porosity_table_rock_control_type
     !! Controls rock porosity via a table of values vs. time.
   contains
     procedure, public :: update => porosity_table_rock_control_update
  end type porosity_table_rock_control_type

contains
  
!------------------------------------------------------------------------

  subroutine permeability_table_rock_control_update(self, time, &
       vector_array, section, range_start)
    !! Update routine for permeability rock control.

    use dm_utils_module, only: global_section_offset

    class(permeability_table_rock_control_type), intent(in out) :: self
    PetscReal, intent(in) :: time
    PetscReal, pointer, contiguous, intent(in) :: vector_array(:)
    PetscSection, intent(in) :: section
    PetscInt, intent(in) :: range_start
    ! Locals:
    PetscReal :: k(self%table%dim), permeability(3)
    type(rock_type) :: rock
    PetscInt :: i, c, rock_offset

    k = self%table%interpolate(time)
    select case (self%table%dim)
    case (1)
       permeability = k(1) ! scalar permeabilities
    case default
       permeability(1: self%table%dim) = k
    end select

    call rock%init()

    do i = 1, size(self%indices)
       c = self%indices(i)
       rock_offset = global_section_offset(section, c, range_start)
       call rock%assign(vector_array, rock_offset)
       rock%permeability = permeability
    end do

    call rock%destroy()

  end subroutine permeability_table_rock_control_update

!------------------------------------------------------------------------

  subroutine porosity_table_rock_control_update(self, time, &
       vector_array, section, range_start)
    !! Update routine for porosity rock control.

    use dm_utils_module, only: global_section_offset

    class(porosity_table_rock_control_type), intent(in out) :: self
    PetscReal, intent(in) :: time
    PetscReal, pointer, contiguous, intent(in) :: vector_array(:)
    PetscSection, intent(in) :: section
    PetscInt, intent(in) :: range_start
    ! Locals:
    PetscReal :: porosity
    type(rock_type) :: rock
    PetscInt :: i, c, rock_offset

    porosity = self%table%interpolate(time, 1)
    call rock%init()

    do i = 1, size(self%indices)
       c = self%indices(i)
       rock_offset = global_section_offset(section, c, range_start)
       call rock%assign(vector_array, rock_offset)
       rock%porosity = porosity
    end do

    call rock%destroy()

  end subroutine porosity_table_rock_control_update

!------------------------------------------------------------------------

end module rock_control_module
