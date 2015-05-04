module rock_module
  !! Defines type for accessing local rock properties on cells and faces.

  use kinds_module

  implicit none
  private

#include <petsc-finclude/petscdef.h>

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
     procedure, public :: destroy => rock_destroy
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
  PetscInt, parameter :: max_rockname_length = 24

  ! Default rock properties:
  PetscReal, parameter, public :: default_porosity = 0.1_dp, default_density = 2200.0_dp
  PetscReal, parameter, public :: default_specific_heat = 1000._dp
  PetscReal, parameter, public :: default_heat_conductivity = 2.5_dp
  PetscReal, parameter, public :: default_permeability(3) = [1.e-13_dp, 1.e-13_dp, 1.e-13_dp]

  public :: rock_type, rock_variable_num_components, &
         rock_variable_dim, rock_variable_names, max_rockname_length

contains

!------------------------------------------------------------------------

  subroutine rock_assign(self, data, offset)
    !! Assigns pointers in a rock object to elements of the specified
    !! data array, starting from the given offset.

    class(rock_type), intent(in out) :: self
    PetscReal, target, intent(in) :: data(:)  !! rock data array
    PetscInt, intent(in) :: offset !! rock array offset

    self%permeability => data(offset: offset + 2)
    self%heat_conductivity => data(offset + 3)
    self%porosity => data(offset + 4)
    self%density => data(offset + 5)
    self%specific_heat => data(offset + 6)

  end subroutine rock_assign
    
!------------------------------------------------------------------------

  PetscInt function rock_dof(self)
    !! Returns degrees of freedom in a rock object.

    class(rock_type), intent(in) :: self

    rock_dof = sum(rock_variable_num_components)

  end function rock_dof

!------------------------------------------------------------------------

  subroutine rock_destroy(self)
    !! Destroys a rock object (nullifies all pointer components).

    class(rock_type), intent(in out) :: self

    nullify(self%permeability)
    nullify(self%heat_conductivity)
    nullify(self%porosity)
    nullify(self%density)
    nullify(self%specific_heat)

  end subroutine rock_destroy

!------------------------------------------------------------------------

end module rock_module
