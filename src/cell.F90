module cell_module
  !! Defines types for accessing local quantities defined on a cell- geometry, rock
  !! and fluid properties.
  !! The components of these types all point to values in arrays obtained from
  !! global parallel vectors.

  use rock_module
  use fluid_module

#include <petsc-finclude/petscdef.h>

  implicit none
  private

  type cell_type
     !! Type for accessing local cell geometry, rock and
     !! fluid properties.
     private
     PetscReal, pointer, public :: volume       !! cell volume
     PetscReal, pointer, public :: centroid(:)  !! cell centroid
     type(rock_type), public :: rock   !! rock properties
     type(fluid_type), public :: fluid !! fluid properties
   contains
     private
     procedure, public :: init => cell_init
     procedure, public :: assign => cell_assign
     procedure, public :: destroy => cell_destroy
     procedure, public :: dof => cell_dof
  end type cell_type

  public :: cell_type

contains

!------------------------------------------------------------------------

  subroutine cell_init(self, num_components, num_phases)
    !! Initialises a cell.

    class(cell_type), intent(in out) :: self
    PetscInt, intent(in) :: num_components !! Number of fluid components
    PetscInt, intent(in) :: num_phases     !! Number of fluid phases

    call self%fluid%init(num_components, num_phases)

  end subroutine cell_init

!------------------------------------------------------------------------

  subroutine cell_assign(self, geom_data, geom_offset, &
       rock_data, rock_offset, fluid_data, fluid_offset)
    !! Assigns pointers in a cell to values from the specified data arrays,
    !! starting from the given offset.

    class(cell_type), intent(in out) :: self
    PetscReal, target, intent(in) :: geom_data(:)  !! array with geometry data
    PetscInt, intent(in)  :: geom_offset  !! geometry array offset for this cell
    PetscReal, target, intent(in), optional :: rock_data(:)  !! array with rock data
    PetscInt, intent(in), optional  :: rock_offset  !! rock array offset for this cell
    PetscReal, target, intent(in), optional :: fluid_data(:)  !! array with fluid data
    PetscInt, intent(in), optional  :: fluid_offset  !! fluid array offset for this cell

    self%volume => geom_data(geom_offset + 1)
    self%centroid => geom_data(geom_offset + 2: geom_offset + 4)

    if ((present(rock_data)) .and. ((present(rock_offset)))) then
       call self%rock%assign(rock_data, rock_offset)
    end if
    if ((present(fluid_data)) .and. ((present(fluid_offset)))) then
       call self%fluid%assign(fluid_data, fluid_offset)
    end if

  end subroutine cell_assign

!------------------------------------------------------------------------

  subroutine cell_destroy(self)
    !! Destroys a cell.

    class(cell_type), intent(in out) :: self
    
    call self%fluid%destroy()
    
  end subroutine cell_destroy

!------------------------------------------------------------------------

  PetscInt function cell_dof(self)
    !! Returns number of degrees of freedom in a cell object.

    class(cell_type), intent(in) :: self
    ! Locals:
    PetscInt, parameter :: fixed_dof = 4

    cell_dof = fixed_dof

  end function cell_dof

!------------------------------------------------------------------------
  
end module cell_module
