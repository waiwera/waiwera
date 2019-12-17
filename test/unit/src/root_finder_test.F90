module root_finder_test

  ! Tests for root finder

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use zofu
  use root_finder_module

  implicit none
  private

  public :: setup, teardown
  public :: test_root_finder_linear, test_root_finder_quadratic, &
       test_root_finder_Zhang, test_root_finder_inverse_quadratic, &
       test_root_finder_saturation

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

  subroutine test_root_finder_linear(test)
    ! Linear equation

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(root_finder_type) :: finder
    PetscMPIInt :: rank
    PetscInt :: ierr
    procedure(root_finder_function), pointer :: f
    PetscReal, parameter :: expected_root = 0.5_dp
    PetscInt, parameter :: expected_iterations = 2

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    f => linear
    call finder%init(f)
    call finder%find()
    if (rank == 0) then
       call test%assert(0, finder%err, "Linear function error")
       call test%assert(expected_root, finder%root, &
            "Linear function root", finder%root_tolerance / finder%root)
       call test%assert(finder%iterations <= expected_iterations, &
            "Linear function iterations")
    end if
    call finder%destroy()

    ! Non-bracketing interval:
    call finder%init(f, [0.75_dp, 1._dp])
    call finder%find()
    if (rank == 0) then
       call test%assert(ROOT_FINDER_INTERVAL_NOT_BRACKETED, &
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

  subroutine test_root_finder_quadratic(test)
    ! Quadratic equation

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(root_finder_type) :: finder
    PetscMPIInt :: rank
    PetscInt :: ierr
    procedure(root_finder_function), pointer :: f
    PetscReal, parameter :: expected_root = 0.75_dp - sqrt(0.5_dp)
    PetscInt, parameter :: expected_iterations = 7

    f => quadratic
    call finder%init(f)
    call finder%find()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert(0, finder%err, "Quadratic function error")
       call test%assert(expected_root, finder%root, &
            "Quadratic function root", finder%root_tolerance / finder%root)
       call test%assert(finder%iterations <= expected_iterations, &
            "Quadratic function iterations")
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

  subroutine test_root_finder_Zhang(test)
    ! Zhang function

    ! Zhang (2011), "An Improvement to the Brent's Method", IJEA,
    ! vol. 2, pp. 21-26.

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(root_finder_type) :: finder
    PetscMPIInt :: rank
    PetscInt :: ierr
    procedure(root_finder_function), pointer :: f
    PetscReal, parameter :: expected_root = 0.8654740331015734_dp
    PetscInt, parameter :: expected_iterations = 12

    f => zhang
    call finder%init(f, [0._dp, 4._dp])
    call finder%find()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert(0, finder%err, "Zhang function error")
       call test%assert(expected_root, finder%root, &
            "Zhang function root", finder%root_tolerance / finder%root)
       call test%assert(finder%iterations <= expected_iterations, &
            "Zhang function iterations")
    end if

    call finder%destroy()

  contains

    PetscReal function zhang(x, context) result(y)
      PetscReal, intent(in) :: x
      class(*), pointer, intent(in out) :: context
      y = cos(x) - x ** 3
    end function zhang

  end subroutine test_root_finder_Zhang

!------------------------------------------------------------------------

  subroutine test_root_finder_inverse_quadratic(test)
    ! Inverse quadratic function

    ! From Stage (2013), "Comments on An Improvement to Brent's Method", IJEA,
    ! vol. 4(1).

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(root_finder_type) :: finder
    PetscMPIInt :: rank
    PetscInt :: ierr
    procedure(root_finder_function), pointer :: f
    PetscReal, parameter :: expected_root = 2._dp / 3._dp
    PetscInt, parameter :: expected_iterations = 18

    f => invquad
    call finder%init(f, [-10._dp, 10._dp])
    call finder%find()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert(0, finder%err, "Inverse quadratic function error")
       call test%assert(expected_root, finder%root, &
            "Inverse quadratic function root", finder%root_tolerance / finder%root)
       call test%assert(finder%iterations <= expected_iterations, &
            "Inverse quadratic function iterations")
    end if

    call finder%destroy()

  contains

    PetscReal function invquad(x, context) result(y)
      PetscReal, intent(in) :: x
      class(*), pointer, intent(in out) :: context
      associate(xs => x - 2._dp / 3._dp)
        y = sqrt(abs(xs))
        if (xs > 0._dp) then
           y = -y
        end if
      end associate
    end function invquad

  end subroutine test_root_finder_inverse_quadratic

!------------------------------------------------------------------------

  subroutine test_root_finder_saturation(test)
    ! Saturation line intersection

    use IAPWS_module
    use interpolation_module, only: interpolation_table_type

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(root_finder_type) :: finder
    type(interpolation_table_type), target :: inc
    class(*), pointer :: pinc
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscErrorCode :: err
    procedure(root_finder_function), pointer :: f
    type(IAPWS_type) :: thermo
    PetscInt, parameter :: num_vars = 2
    PetscReal :: var(num_vars), data(2, 1 + num_vars)
    PetscReal, parameter :: expected_temperature = 218.61315743282924_dp
    PetscInt, parameter :: expected_iterations = 6

    call thermo%init()
    data(:, 1) = [0._dp, 1._dp]
    data(:, 2) = [20.e5_dp, 23.e5_dp]
    data(:, 3) = [210._dp, 220._dp]
    call inc%init(data, err = err)
    call test%assert(0, err, "error")

    pinc => inc
    f => saturation_difference
    call finder%init(f, context = pinc)
    call finder%find()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert(0, finder%err, "Saturation line error")
       var = inc%interpolate(finder%root)
       associate(T => var(2))
         call test%assert(expected_temperature, T, &
              "Saturation line temperature", finder%root_tolerance / finder%root)
       end associate
       call test%assert(finder%iterations <= expected_iterations, &
            "Saturation line iterations")
    end if

    call finder%destroy()
    call inc%destroy()
    call thermo%destroy()

  contains

    PetscReal function saturation_difference(x, context) result(dp)
      ! Returns pressure difference between point 0 <= x <= 1 along
      ! line between (P0, T0) and (P1, T1) and saturation curve.
      PetscReal, intent(in) :: x
      class(*), pointer, intent(in out) :: context
      ! Locals:
      PetscReal :: var(num_vars), Ps
      PetscInt :: err
      associate(P => var(1), T => var(2))
        select type (context)
        type is (interpolation_table_type)
           var = context%interpolate_at_index(x)
           call thermo%saturation%pressure(T, Ps, err)
        end select
        dp = Ps - P
      end associate
    end function saturation_difference

  end subroutine test_root_finder_saturation

!------------------------------------------------------------------------

end module root_finder_test
