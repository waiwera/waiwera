module root_finder_test

  ! Tests for root finder

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use fruit
  use root_finder_module

  implicit none
  private

  public :: test_root_finder_linear, test_root_finder_quadratic

contains

!------------------------------------------------------------------------

  subroutine test_root_finder_linear
    ! Test root finder on linear equation.

    type(root_finder_type) :: finder
    PetscMPIInt :: rank
    PetscInt :: ierr
    procedure(root_finder_function), pointer :: f
    PetscReal, parameter :: expected_root = 0.5_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    f => linear
    call finder%init(f)
    call finder%find()
    if (rank == 0) then
       call assert_equals(0, finder%err, "Linear function error")
       call assert_equals(expected_root, finder%root, finder%root_tolerance, &
            "Linear function root")
       call assert_true(finder%iterations <= 2, "Linear function iterations")
    end if
    call finder%destroy()

    ! Non-bracketing interval:
    call finder%init(f, [0.75_dp, 1._dp])
    call finder%find()
    if (rank == 0) then
       call assert_equals(ROOT_FINDER_INTERVAL_NOT_BRACKETED, &
            finder%err, "Linear function non-bracketing error")
    end if
    call finder%destroy()

  contains

    PetscReal function linear(x, context) result(y)
      PetscReal, intent(in) :: x
      class(*), pointer, intent(in out) :: context
      y = 0.5_dp - x
    end function linear

  end subroutine test_root_finder_linear

!------------------------------------------------------------------------

  subroutine test_root_finder_quadratic
    ! Test root finder on quadratic equation.

    type(root_finder_type) :: finder
    PetscMPIInt :: rank
    PetscInt :: ierr
    procedure(root_finder_function), pointer :: f
    PetscReal, parameter :: expected_root = 0.75_dp - sqrt(0.5_dp)

    f => quadratic
    call finder%init(f)
    call finder%find()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call assert_equals(0, finder%err, "Quadratic function error")
       call assert_equals(expected_root, finder%root, &
            finder%root_tolerance, "Quadratic function root")
       call assert_true(finder%iterations <= 7, "Quadratic function iterations")
    end if

    call finder%destroy()

  contains

    PetscReal function quadratic(x, context) result(y)
      PetscReal, intent(in) :: x
      class(*), pointer, intent(in out) :: context
      y = (x - 0.75_dp) ** 2 - 0.5_dp
    end function quadratic

  end subroutine test_root_finder_quadratic

!------------------------------------------------------------------------

end module root_finder_test
