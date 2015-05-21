module boundary_module
  !! Module for implementation of simulation boundary conditions.

  implicit none
  private

#include <petsc-finclude/petsc.h90>

  character(len = 16), public :: open_boundary_label_name = "open boundary"

  public :: setup_boundary_labels

contains

!------------------------------------------------------------------------

  subroutine setup_boundary_labels(dm)
    !! Sets up labels identifying boundaries of the mesh (for e.g.
    !! applying boundary conditions.

    DM, intent(in) :: dm
    ! Locals:
    PetscErrorCode :: ierr
    PetscBool :: has_label

    call DMPlexHasLabel(dm, open_boundary_label_name, has_label, &
         ierr); CHKERRQ(ierr)
    if (.not.(has_label)) then
       call DMPlexCreateLabel(dm, open_boundary_label_name, &
            ierr); CHKERRQ(ierr)
       ! could read boundary faces from input here if needed- i.e. if labels
       ! not present in mesh file
    end if

  end subroutine setup_boundary_labels

!------------------------------------------------------------------------

end module boundary_module
