module setup_module

  ! Setup and teardown routines for unit tests

  use mpi_module

  implicit none

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdef.h>

contains

!------------------------------------------------------------------------

  subroutine setup

    ! Global setup routine

    PetscErrorCode :: ierr
    
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call mpi%init(PETSC_COMM_WORLD)
    
  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown

    ! Global teardown routine

    PetscErrorCode :: ierr

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

end module setup_module
