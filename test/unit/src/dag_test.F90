module dag_test

  ! Tests for dag module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use dag_module
  use fruit

  implicit none
  private

  public :: test_dag

contains

!------------------------------------------------------------------------

  subroutine test_dag
    ! DAG

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
       call assert_equals(0, err, "sort error")
       call assert_equals([0, 1, 4, 2, 3], order, 5, "sort order")
    end if

    call dag%set_edges(4, [2])
    call dag%sort(order, err)
    if (rank == 0) then
       call assert_equals(1, err, "circular dependency")
    end if

  end subroutine test_dag

!------------------------------------------------------------------------

end module dag_test
