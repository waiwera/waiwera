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

  type, abstract, extends(rock_control_type), public :: rock_control_table_type
     !! Controls a rock parameter (e.g. permeability) via a table of
     !! values vs. time.
     private
     PetscInt, allocatable, public :: cell_indices(:)
     type(interpolation_table_type), public :: table !! Table of values vs. time
   contains
     private
     procedure, public :: init => rock_control_table_init
     procedure, public :: destroy => rock_control_table_destroy
  end type rock_control_table_type

  type, extends(rock_control_table_type), public :: rock_control_permeability_table_type
     !! Controls rock permeability via a table of values vs. time.
   contains
     private
     procedure, public :: update => rock_control_permeability_table_update
  end type rock_control_permeability_table_type

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

contains

!------------------------------------------------------------------------
! Table rock control:
!------------------------------------------------------------------------

  subroutine rock_control_table_init(self, data, interpolation_type, &
       averaging_type, cell_indices, err)
    !! Initialises rock_control_table type.

    class(rock_control_table_type), intent(in out) :: self
    PetscReal, intent(in) :: data(:,:)
    PetscInt, intent(in) :: interpolation_type
    PetscInt, intent(in) :: averaging_type
    PetscInt, intent(in) :: cell_indices(:)
    PetscErrorCode, intent(out) :: err

    call self%table%init(data, interpolation_type, averaging_type, err)
    if (err == 0) then
       self%cell_indices = cell_indices
    end if

  end subroutine rock_control_table_init

!------------------------------------------------------------------------

  subroutine rock_control_table_destroy(self)
    !! Destroys rock_control_table_type object.

    class(rock_control_table_type), intent(in out) :: self

    call self%table%destroy()
    if (allocated(self%cell_indices)) deallocate(self%cell_indices)

  end subroutine rock_control_table_destroy

!------------------------------------------------------------------------

  subroutine rock_control_permeability_table_update(self, t, &
       rock_data, rock_section, rock_range_start)
    !! Update permeability for rock_control_permeabiltiy_table_type.

    use dm_utils_module, only: global_section_offset

    class(rock_control_permeability_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    PetscReal, pointer, contiguous, intent(in) :: rock_data(:)
    PetscSection, intent(in) :: rock_section
    PetscInt, intent(in) :: rock_range_start
    ! Locals:
    PetscReal :: k(self%table%dim), permeability(3)
    type(rock_type) :: rock
    PetscInt :: i, c, rock_offset

    k = self%table%interpolate(t)
    select case (self%table%dim)
    case (1)
       permeability = k(1) ! scalar permeabilities
    case default
       permeability(1: self%table%dim) = k
    end select

    call rock%init()

    do i = 1, size(self%cell_indices)
       c = self%cell_indices(i)
       rock_offset = global_section_offset(rock_section, c, &
            rock_range_start)
       call rock%assign(rock_data, rock_offset)
       rock%permeability = permeability
    end do

    call rock%destroy()

  end subroutine rock_control_permeability_table_update

!------------------------------------------------------------------------

end module rock_control_module
