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
     PetscReal, pointer, public :: normal(:) !! normal vector to face
     PetscReal, pointer, public :: area_normal(:) !! area-weighted normal vector (normal / area)
     type(cell_type), public :: cell(2)
   contains
     private
     procedure, public :: assign => face_assign
     procedure, public :: dof => face_dof
  end type face_type

  public :: face_type

contains

!------------------------------------------------------------------------

  subroutine face_assign(self, face_geom_data, face_geom_offset, cell_geom_data)
    !! Assigns pointers in a face to elements of the specified data array,
    !! starting from the given offset.

    class(face_type), intent(in out) :: self
    PetscReal, target, intent(in) :: face_geom_data(:)  !! array with face geometry data
    PetscInt, intent(in)  :: face_geom_offset  !! face geometry array offset for this cell
    PetscReal, target, intent(in) :: cell_geom_data(:)  !! array with cell geometry data
    
    self%area_normal => face_geom_data(face_geom_offset: face_geom_offset + 2)
    self%normal => face_geom_data(face_geom_offset + 3: face_geom_offset + 5)
    self%area => face_geom_data(face_geom_offset + 6)
    self%distance => face_geom_data(face_geom_offset + 7: face_geom_offset + 8)

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
