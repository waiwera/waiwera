module profiling_module
  !! Module for profiling code via PETSc log events.

  implicit none
  private

#include <petsc/finclude/petscsys.h>

  PetscClassId, public ::  log_class
  PetscLogDouble, public :: flops
  PetscLogEvent, public :: lhs_fn_event
  PetscLogEvent, public :: rhs_fn_event

  public :: init_profiling

contains

!------------------------------------------------------------------------

  subroutine init_profiling()
    !! Initialize code profiling.

    PetscErrorCode :: ierr

    call PetscClassIdRegister("supermodel", log_class, ierr); CHKERRQ(ierr)

    ! Register log events:
    call PetscLogEventRegister("lhs_function", log_class, lhs_fn_event, ierr)
    call PetscLogEventRegister("rhs_function", log_class, rhs_fn_event, ierr)

  end subroutine init_profiling

!------------------------------------------------------------------------

end module profiling_module
