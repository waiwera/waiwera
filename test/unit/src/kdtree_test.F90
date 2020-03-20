module kdtree_test

  ! Tests for k-d tree module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use zofu
  use kdtree_module

  implicit none
  private

  public :: setup, teardown
  public :: test_kdtree_linear

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

  subroutine test_kdtree_linear(test)

    !! Test k-d tree on linear mesh

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call linear_test_case(2)
       call linear_test_case(3)

    end if

  contains

    subroutine linear_test_case(dim)

      !! test case for dimension dim

      PetscInt, intent(in) :: dim
      ! Locals:
      type(kdtree_type) :: kdt
      PetscReal, allocatable :: data(:,:), x(:)
      PetscInt :: n, i
      PetscInt, allocatable :: index(:), expected_index(:)
      PetscErrorCode, allocatable :: err(:)
      character(12) :: msg

      write(msg, '(a, i0)') 'linear dim ', dim
      n = 1000
      allocate(data(dim, n), index(n), expected_index(n), err(n))
      do i = 1, n
         data(1, i) = i - 0.5_dp
         data(2:, i) = 0.5_dp
      end do

      call kdt%init(data)

      do i = 1, n
         x = data(:, i)
         call kdt%search(x, index(i), err = err(i))
         expected_index(i) = i
      end do

      call test%assert(index, expected_index, msg // ' indices')
      call test%assert(all(err == 0), msg // ' error')

      call kdt%destroy()
      deallocate(data, x, index, expected_index, err)

    end subroutine linear_test_case

  end subroutine test_kdtree_linear

!------------------------------------------------------------------------

end module kdtree_test
