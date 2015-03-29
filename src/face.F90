module face_module
  !! Defines type for accessing local quantities defined on a mesh face.

#include <petsc-finclude/petscdef.h>

  use cell_module

  implicit none
  private

  type face_type
     !! Type for accessing local face properties.
     private
     PetscReal, pointer, public :: area !! face area
     PetscReal, pointer, public :: distance(:) !! cell centroid distances on either side of the face
     type(cell_type), public :: cell(2)
   contains
     private
     procedure, public :: assign => face_assign
     procedure, public :: dof => face_dof
  end type face_type

  public :: face_type

contains

!------------------------------------------------------------------------

  subroutine face_assign(self, geom_data, geom_offset)
    !! Assigns pointers in a face to elements of the specified data array,
    !! starting from the given offset.

    class(face_type), intent(in out) :: self
    PetscReal, target, intent(in) :: geom_data(:)  !! array with face geometry data
    PetscInt, intent(in)  :: geom_offset  !! face geometry array offset for this cell
    
    self%area => geom_data(geom_offset + 1)
    self%distance => geom_data(geom_offset + 2: geom_offset + 3)

    ! TODO: assign self%cell(:) here

  end subroutine face_assign

!------------------------------------------------------------------------

  PetscInt function face_dof(self)
    !! Returns number of degrees of freedom in a face object.

    class(face_type), intent(in) :: self
    ! Locals:
    PetscInt, parameter :: fixed_dof = 3

    face_dof = fixed_dof

  end function face_dof

!------------------------------------------------------------------------

end module face_module
