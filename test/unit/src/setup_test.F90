module setup_module

  ! Setup and teardown routines for unit tests

#include <petsc/finclude/petscsys.h>

  use petscsys

  implicit none

contains

!------------------------------------------------------------------------

  subroutine setup

    ! Global setup routine

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr
    
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()
    
  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown

    ! Global teardown routine

    PetscErrorCode :: ierr

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

end module setup_module
