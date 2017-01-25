module powertable_test

  ! Tests for powertable module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use fruit
  use powertable_module

  implicit none
  private 

  PetscReal, parameter :: tol = 1.e-9

  public :: test_powertable_positive_powers, &
       test_powertable_negative_powers, test_powertable_mixed_powers, &
       test_powertable_multiple_configure

  contains

!------------------------------------------------------------------------

    subroutine test_powertable_positive_powers

      ! Powertable positive powers test

      type(powertable_type) :: p
      PetscInt, parameter :: powers(4) = [1, 2, 3, 4]
      PetscReal :: x
      PetscReal, allocatable :: xp(:)
      PetscInt :: lower, upper

      if (mpi%rank == mpi%output_rank) then

         lower = min(minval(powers), 0)
         upper = max(maxval(powers), 1)

         allocate(xp(lower:upper))

         call p%configure(powers)

         call assert_equals(lower, p%lower, "lower")
         call assert_equals(upper, p%upper, "upper")

         x = 1._dp
         xp = [1._dp, 1._dp, 1._dp, 1._dp]
         call p%compute(x)
         call assert_equals(xp, p%power(powers), size(xp), tol, "1.0")

         x = 2._dp
         xp = [2._dp, 4._dp, 8._dp, 16._dp]
         call p%compute(x)
         call assert_equals(xp, p%power(powers), size(xp), tol, "2.0")

         deallocate(xp)
         call p%destroy()

      end if

    end subroutine test_powertable_positive_powers

!------------------------------------------------------------------------

    subroutine test_powertable_negative_powers

      ! Powertable negative powers test

      type(powertable_type) :: p
      PetscInt, parameter :: powers(3) = [-6, -3, -2]
      PetscReal :: x
      PetscReal, allocatable :: xp(:)
      PetscInt :: lower, upper

      if (mpi%rank == mpi%output_rank) then

         lower = min(minval(powers), 0)
         upper = max(maxval(powers), 1)

         allocate(xp(lower:upper))

         call p%configure(powers)

         call assert_equals(lower, p%lower, "lower")
         call assert_equals(upper, p%upper, "upper")

         x = 1._dp
         xp = [1._dp, 1._dp, 1._dp]
         call p%compute(x)
         call assert_equals(xp, p%power(powers), size(xp), tol, "1.0")

         x = 2._dp
         xp = [0.015625_dp, 0.125_dp, 0.25_dp]
         call p%compute(x)
         call assert_equals(xp, p%power(powers), size(xp), tol, "2.0")

         deallocate(xp)
         call p%destroy()

      end if

    end subroutine test_powertable_negative_powers

!------------------------------------------------------------------------

    subroutine test_powertable_mixed_powers

      ! Powertable mixed powers test

      type(powertable_type) :: p
      PetscInt, parameter :: powers(4) = [-3, -1, 1, 4]
      PetscReal :: x
      PetscReal, allocatable :: xp(:)
      PetscInt :: lower, upper

      if (mpi%rank == mpi%output_rank) then
         
         lower = min(minval(powers), 0)
         upper = max(maxval(powers), 1)

         allocate(xp(lower:upper))

         call p%configure(powers)

         call assert_equals(lower, p%lower, "lower")
         call assert_equals(upper, p%upper, "upper")

         x = 1._dp
         xp = [1._dp, 1._dp, 1._dp, 1._dp]
         call p%compute(x)
         call assert_equals(xp, p%power(powers), size(xp), tol, "1.0")

         x = 2._dp
         xp = [0.125_dp, 0.5_dp, 2._dp, 16._dp]
         call p%compute(x)
         call assert_equals(xp, p%power(powers), size(xp), tol, "2.0")

         deallocate(xp)
         call p%destroy()

      end if

    end subroutine test_powertable_mixed_powers

!------------------------------------------------------------------------

    subroutine test_powertable_multiple_configure

      ! Powertable multiple configuration test

      type(powertable_type) :: p
      PetscInt, parameter :: powers1(3) = [-4, -3, -1], powers2(2) = [3, 5]
      PetscInt, allocatable :: allpowers(:)
      PetscReal :: x
      PetscReal, allocatable :: xp(:)
      PetscInt :: lower, upper

      if (mpi%rank == mpi%output_rank) then
         
         allpowers = [powers1, powers2]
         lower = min(minval(allpowers), 0)
         upper = max(maxval(allpowers), 1)

         allocate(xp(lower:upper))

         call p%configure(powers1)
         call p%configure(powers2)

         call assert_equals(lower, p%lower, "lower")
         call assert_equals(upper, p%upper, "upper")

         x = 1._dp
         xp = [1._dp, 1._dp, 1._dp, 1._dp, 1._dp]
         call p%compute(x)
         call assert_equals(xp, p%power(allpowers), size(xp), tol, "1.0")

         x = 2._dp
         xp = [0.0625_dp, 0.125_dp, 0.5_dp, 8._dp, 32._dp]
         call p%compute(x)
         call assert_equals(xp, p%power(allpowers), size(xp), tol, "2.0")

         deallocate(xp, allpowers)
         call p%destroy()

      end if

    end subroutine test_powertable_multiple_configure

!------------------------------------------------------------------------

end module powertable_test
