module profiling_module
  !! Module for profiling code via PETSc log events.

  implicit none
  private

#include <petsc/finclude/petscsys.h>

  PetscClassId, public ::  log_class
  PetscLogDouble, public :: flops
  PetscLogEvent, public :: simulation_init_event
  PetscLogEvent, public :: fluid_init_event
  PetscLogEvent, public :: fluid_properties_event, fluid_transitions_event
  PetscLogEvent, public :: lhs_fn_event, rhs_fn_event
  PetscLogEvent, public :: output_event

  public :: init_profiling

contains

!------------------------------------------------------------------------

  subroutine init_profiling()
    !! Initialize code profiling.

    ! Locals:
    PetscErrorCode :: ierr

    call PetscClassIdRegister("supermodel", log_class, ierr); CHKERRQ(ierr)

    ! Register log events:
    call PetscLogEventRegister("sim_init", log_class, simulation_init_event, ierr)
    call PetscLogEventRegister("fluid_init", log_class, fluid_init_event, ierr)
    call PetscLogEventRegister("fluid_properties", log_class, fluid_properties_event, ierr)
    call PetscLogEventRegister("fluid_transitions", log_class, fluid_transitions_event, ierr)
    call PetscLogEventRegister("lhs_function", log_class, lhs_fn_event, ierr)
    call PetscLogEventRegister("rhs_function", log_class, rhs_fn_event, ierr)
    call PetscLogEventRegister("output", log_class, output_event, ierr)

  end subroutine init_profiling

!------------------------------------------------------------------------

end module profiling_module
