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
  use eos_module

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

  type, public, abstract, extends(vector_vector_control_type) :: fluid_rock_control_type
     !! Controls rock properties according to fluid properties.
     private
     class(eos_type), pointer :: eos !! Equation of state
   contains
     procedure :: local_update => fluid_rock_control_local_update
     procedure, public :: update => fluid_rock_control_update
  end type fluid_rock_control_type

contains
  
!------------------------------------------------------------------------
! Permeability table rock control
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
! Porosity table rock control
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
! Fluid rock control
!------------------------------------------------------------------------

  subroutine fluid_rock_control_local_update(self, fluid, rock)
    !! Updates rock object according to fluid object. Derived types
    !! override this procedure.

    use fluid_module, only: fluid_type

    class(fluid_rock_control_type), intent(in out) :: self
    type(fluid_type), intent(in) :: fluid
    type(rock_type), intent(in out) :: rock

    continue

  end subroutine fluid_rock_control_local_update

!------------------------------------------------------------------------

  subroutine fluid_rock_control_update(self, &
       subject_array, subject_section, subject_range_start, &
       object_array, object_section, object_range_start)
    !! Updates rock vector based on the fluid vector.

    use dm_utils_module, only: global_section_offset
    use fluid_module, only: fluid_type

    class(fluid_rock_control_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: subject_array(:) !! Array on fluid vector
    PetscSection, intent(in) :: subject_section !! Global section for fluid vector
    PetscInt, intent(in) :: subject_range_start !! Range start for fluid vector
    PetscReal, pointer, contiguous, intent(in out) :: object_array(:) !! Array on rock vector
    PetscSection, intent(in) :: object_section !! Global section for rock vector
    PetscInt, intent(in) :: object_range_start !! Range start for rock vector
    ! Locals:
    PetscInt :: i, c, fluid_offset, rock_offset
    type(rock_type) :: rock
    type(fluid_type) :: fluid

    call fluid%init(self%eos%num_components, self%eos%num_phases)
    call rock%init()

    do i = 1, size(self%indices)
       c = self%indices(i)
       fluid_offset = global_section_offset(subject_section, c, subject_range_start)
       rock_offset = global_section_offset(object_section, c, object_range_start)
       call fluid%assign(subject_array, fluid_offset)
       call rock%assign(object_array, rock_offset)
       call self%local_update(fluid, rock)
    end do

    call fluid%destroy()
    call rock%destroy()

  end subroutine fluid_rock_control_update

!------------------------------------------------------------------------

end module rock_control_module
