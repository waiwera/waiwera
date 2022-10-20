module salt_thermodynamics_test

  ! Tests for salt thermodynamics module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use salt_thermodynamics_module
  use zofu

  implicit none
  private

  public :: setup, teardown
  public :: test_halite_solubility

contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown()

    PetscErrorCode :: ierr

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

  subroutine test_halite_solubility(test)
    ! Halite solubility

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscReal :: temperature, expected, s
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       temperature = 20._dp
       expected = 0.264044_dp
       s = halite_solubility(temperature)
       call test%assert(expected, s, " 20 deg C")

       temperature = 200._dp
       expected = 0.31898_dp
       s = halite_solubility(temperature)
       call test%assert(expected, s, " 200 deg C")

       temperature = 300._dp
       expected = 0.37918_dp
       s = halite_solubility(temperature)
       call test%assert(expected, s, " 300 deg C")

       temperature = 350._dp
       expected = 0.41723_dp
       s = halite_solubility(temperature)
       call test%assert(expected, s, " 350 deg C")

    end if

  end subroutine test_halite_solubility

!------------------------------------------------------------------------

end module salt_thermodynamics_test
