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
     procedure, public :: populate => rock_populate
  end type rock_type

  public :: rock_type

contains

!------------------------------------------------------------------------

  subroutine rock_populate(self, data, offset)
    !! Populates a rock object with values from the specified data array,
    !! starting from the given offset.

    class(rock_type), intent(in out) :: self
    PetscReal, target, intent(in) :: data(:)  !! rock data array
    PetscInt, intent(in) :: offset !! rock array offset

    self%permeability => data(offset + 1: offset + 3)
    self%heat_conductivity => data(offset + 4)
    self%porosity => data(offset + 5)
    self%density => data(offset + 6)
    self%specific_heat => data(offset + 7)

  end subroutine rock_populate
    
!------------------------------------------------------------------------

end module rock_module
