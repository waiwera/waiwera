module setup_module

  ! Setup and teardown routines for unit tests

  use IAPWS_module
  use IFC67_module

  implicit none

#include <finclude/petscsys.h>
#include <finclude/petscdef.h>

contains

!------------------------------------------------------------------------

  subroutine setup
    
    implicit none
    PetscErrorCode :: ierr
    
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    
    call IAPWS%init()
    call IFC67%init()

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown

    implicit none
    PetscErrorCode :: ierr

    call IAPWS%destroy()
    call IFC67%destroy()

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

end module setup_module
