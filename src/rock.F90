module rock_module
  !! Defines type for accessing local rock properties on cells and faces.

#include <petsc-finclude/petscdef.h>

  implicit none
  private

  type rock_type
     !! Local rock properties.
     PetscReal, pointer :: permeability(:)   !! permeability
     PetscReal, pointer :: heat_conductivity !! heat conductivity
     PetscReal, pointer :: porosity          !! porosity
     PetscReal, pointer :: density           !! grain density
     PetscReal, pointer :: specific_heat     !! specific heat
   contains
     private
     procedure, public :: assign => rock_assign
     procedure, public :: dof => rock_dof
  end type rock_type

  PetscInt, parameter :: num_rock_variables = 5
  PetscInt, parameter :: max_rock_variable_name_length = 32
  character(max_rock_variable_name_length), parameter :: &
       rock_variable_names(num_rock_variables) = &
       [character(max_rock_variable_name_length):: &
       "Permeability", "Heat conductivity", "Porosity", &
       "Density", "Specific heat"]
  PetscInt, parameter :: &
       rock_variable_num_components(num_rock_variables) = &
       [3, 1, 1, 1, 1]
  PetscInt, parameter :: &
       rock_variable_dim(num_rock_variables) = &
       [3, 3, 3, 3, 3]

  public :: rock_type, rock_variable_num_components, &
         rock_variable_dim, rock_variable_names

contains

!------------------------------------------------------------------------

  subroutine rock_assign(self, data, offset)
    !! Assigns pointers in a rock object to elements of the specified
    !! data array, starting from the given offset.

    class(rock_type), intent(in out) :: self
    PetscReal, target, intent(in) :: data(:)  !! rock data array
    PetscInt, intent(in) :: offset !! rock array offset

    self%permeability => data(offset + 1: offset + 3)
    self%heat_conductivity => data(offset + 4)
    self%porosity => data(offset + 5)
    self%density => data(offset + 6)
    self%specific_heat => data(offset + 7)

  end subroutine rock_assign
    
!------------------------------------------------------------------------

  PetscInt function rock_dof(self)
    !! Returns degrees of freedom in a rock object.

    class(rock_type), intent(in) :: self
    ! Locals:
    PetscInt, parameter :: fixed_dof = 7

    rock_dof = fixed_dof

  end function rock_dof

!------------------------------------------------------------------------

end module rock_module
