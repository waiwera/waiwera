module setup_module

  ! Setup and teardown routines for unit tests

  use mpi_module
  use IAPWS_module
  use IFC67_module

  implicit none

#include <petsc-finclude/petscsys.h>
#include <petsc-finclude/petscdef.h>

contains

!------------------------------------------------------------------------

  subroutine setup
    
    PetscErrorCode :: ierr
    
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)

    comm = PETSC_COMM_WORLD
    call MPI_COMM_RANK(comm, rank, ierr)
    
    call IAPWS%init()
    call IFC67%init()

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown

    PetscErrorCode :: ierr

    call IAPWS%destroy()
    call IFC67%destroy()

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

end module setup_module
