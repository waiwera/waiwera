module powertable_test

  ! Tests for powertable module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use zofu
  use powertable_module

  implicit none
  private 

  public :: setup, teardown, setup_test
  public :: test_powertable_positive_powers, &
       test_powertable_negative_powers, test_powertable_mixed_powers, &
       test_powertable_multiple_configure

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

  subroutine setup_test(test)

    class(unit_test_type), intent(in out) :: test

    test%tolerance = 1.e-9

  end subroutine setup_test

!------------------------------------------------------------------------

  subroutine test_powertable_positive_powers(test)

    ! Powertable positive powers test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(powertable_type) :: p
    PetscInt, parameter :: powers(4) = [1, 2, 3, 4]
    PetscReal :: x
    PetscReal, allocatable :: xp(:)
    PetscInt :: lower, upper
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       lower = min(minval(powers), 0)
       upper = max(maxval(powers), 1)

       allocate(xp(lower:upper))

       call p%configure(powers)

       call test%assert(lower, p%lower, "lower")
       call test%assert(upper, p%upper, "upper")

       x = 1._dp
       xp = [1._dp, 1._dp, 1._dp, 1._dp]
       call p%compute(x)
       call test%assert(xp, p%power(powers), "1.0")

       x = 2._dp
       xp = [2._dp, 4._dp, 8._dp, 16._dp]
       call p%compute(x)
       call test%assert(xp, p%power(powers), "2.0")

       deallocate(xp)
       call p%destroy()

    end if

  end subroutine test_powertable_positive_powers

!------------------------------------------------------------------------

  subroutine test_powertable_negative_powers(test)

    ! Powertable negative powers test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(powertable_type) :: p
    PetscInt, parameter :: powers(3) = [-6, -3, -2]
    PetscReal :: x
    PetscReal, allocatable :: xp(:)
    PetscInt :: lower, upper
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       lower = min(minval(powers), 0)
       upper = max(maxval(powers), 1)

       allocate(xp(lower:upper))

       call p%configure(powers)

       call test%assert(lower, p%lower, "lower")
       call test%assert(upper, p%upper, "upper")

       x = 1._dp
       xp = [1._dp, 1._dp, 1._dp]
       call p%compute(x)
       call test%assert(xp, p%power(powers), "1.0")

       x = 2._dp
       xp = [0.015625_dp, 0.125_dp, 0.25_dp]
       call p%compute(x)
       call test%assert(xp, p%power(powers), "2.0")

       deallocate(xp)
       call p%destroy()

    end if

  end subroutine test_powertable_negative_powers

!------------------------------------------------------------------------

  subroutine test_powertable_mixed_powers(test)

    ! Powertable mixed powers test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(powertable_type) :: p
    PetscInt, parameter :: powers(4) = [-3, -1, 1, 4]
    PetscReal :: x
    PetscReal, allocatable :: xp(:)
    PetscInt :: lower, upper
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       lower = min(minval(powers), 0)
       upper = max(maxval(powers), 1)

       allocate(xp(lower:upper))

       call p%configure(powers)

       call test%assert(lower, p%lower, "lower")
       call test%assert(upper, p%upper, "upper")

       x = 1._dp
       xp = [1._dp, 1._dp, 1._dp, 1._dp]
       call p%compute(x)
       call test%assert(xp, p%power(powers), "1.0")

       x = 2._dp
       xp = [0.125_dp, 0.5_dp, 2._dp, 16._dp]
       call p%compute(x)
       call test%assert(xp, p%power(powers), "2.0")

       deallocate(xp)
       call p%destroy()

    end if

  end subroutine test_powertable_mixed_powers

!------------------------------------------------------------------------

  subroutine test_powertable_multiple_configure(test)

    ! Powertable multiple configuration test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(powertable_type) :: p
    PetscInt, parameter :: powers1(3) = [-4, -3, -1], powers2(2) = [3, 5]
    PetscInt, allocatable :: allpowers(:)
    PetscReal :: x
    PetscReal, allocatable :: xp(:)
    PetscInt :: lower, upper
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       allpowers = [powers1, powers2]
       lower = min(minval(allpowers), 0)
       upper = max(maxval(allpowers), 1)

       allocate(xp(lower:upper))

       call p%configure(powers1)
       call p%configure(powers2)

       call test%assert(lower, p%lower, "lower")
       call test%assert(upper, p%upper, "upper")

       x = 1._dp
       xp = [1._dp, 1._dp, 1._dp, 1._dp, 1._dp]
       call p%compute(x)
       call test%assert(xp, p%power(allpowers), "1.0")

       x = 2._dp
       xp = [0.0625_dp, 0.125_dp, 0.5_dp, 8._dp, 32._dp]
       call p%compute(x)
       call test%assert(xp, p%power(allpowers), "2.0")

       deallocate(xp, allpowers)
       call p%destroy()

    end if

  end subroutine test_powertable_multiple_configure

!------------------------------------------------------------------------

end module powertable_test
