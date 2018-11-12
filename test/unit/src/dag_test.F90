module dag_test

  ! Tests for dag module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use dag_module
  use zofu

  implicit none
  private

  public :: setup, teardown, test_dag

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

  subroutine test_dag(test)
    ! DAG

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(dag_type) :: dag
    PetscInt, allocatable :: order(:)
    PetscErrorCode :: err
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    call dag%init(count = 5)

    call dag%set_edges(1, [0])
    call dag%set_edges(2, [4, 0])
    call dag%set_edges(3, [4])

    call dag%sort(order, err)

    if (rank == 0) then
       call test%assert(0, err, "sort error")
       call test%assert([0, 1, 4, 2, 3], order, "sort order")
    end if

    call dag%set_edges(4, [2])
    call dag%sort(order, err)
    if (rank == 0) then
       call test%assert(1, err, "circular dependency")
    end if

  end subroutine test_dag

!------------------------------------------------------------------------

end module dag_test
