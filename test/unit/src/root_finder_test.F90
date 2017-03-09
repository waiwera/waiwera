module root_finder_test

  ! Tests for root finder

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use fruit
  use root_finder_module

  implicit none
  private

  type :: saturation_context_type
     private
     PetscReal, public :: T0, T1, P0, P1
   contains
     private
     procedure, public :: PT => saturation_context_PT
  end type saturation_context_type

  public :: test_root_finder_linear, test_root_finder_quadratic, &
       test_root_finder_Zhang, test_root_finder_inverse_quadratic, &
       test_root_finder_saturation

contains

!------------------------------------------------------------------------

  subroutine saturation_context_PT(self, x, P, T)
    ! Pressure and temperature as functions of x.
    class(saturation_context_type), intent(in) :: self
    PetscReal, intent(in) :: x
    PetscReal, intent(out) :: P, T
    P = (1._dp - x) * self%P0 + x * self%P1
    T = (1._dp - x) * self%T0 + x * self%T1
  end subroutine saturation_context_PT

!------------------------------------------------------------------------

  subroutine test_root_finder_linear
    ! Linear equation.

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
       call assert_equals(0, finder%err, "Linear function error")
       call assert_equals(expected_root, finder%root, finder%root_tolerance, &
            "Linear function root")
       call assert_true(finder%iterations <= expected_iterations, &
            "Linear function iterations")
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
    ! Quadratic equation.

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

       call assert_equals(0, finder%err, "Quadratic function error")
       call assert_equals(expected_root, finder%root, &
            finder%root_tolerance, "Quadratic function root")
       call assert_true(finder%iterations <= expected_iterations, &
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

  subroutine test_root_finder_Zhang
    ! Zhang's function

    ! Zhang (2011), "An Improvement to the Brent’s Method", IJEA,
    ! vol. 2, pp. 21-26.

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

       call assert_equals(0, finder%err, "Zhang function error")
       call assert_equals(expected_root, finder%root, &
            finder%root_tolerance, "Zhang function root")
       call assert_true(finder%iterations <= expected_iterations, &
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

  subroutine test_root_finder_inverse_quadratic
    ! Inverse quadratic function

    ! From Stage (2013), "Comments on An Improvement to Brent's Method", IJEA,
    ! vol. 4(1).

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

       call assert_equals(0, finder%err, "Inverse quadratic function error")
       call assert_equals(expected_root, finder%root, &
            finder%root_tolerance, "Inverse quadratic function root")
       call assert_true(finder%iterations <= expected_iterations, &
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

  subroutine test_root_finder_saturation
    ! Saturation line intersection

    use IAPWS_module

    type(root_finder_type) :: finder
    type(saturation_context_type), target :: context
    class(*), pointer :: pcontext
    PetscMPIInt :: rank
    PetscInt :: ierr
    procedure(root_finder_function), pointer :: f
    type(IAPWS_type) :: thermo
    PetscReal :: P, T
    PetscReal, parameter :: expected_temperature = 218.61315743282924_dp
    PetscInt, parameter :: expected_iterations = 6

    call thermo%init()
    context%P0 = 20.e5
    context%T0 = 210._dp
    context%P1 = 23.e5
    context%T1 = 220._dp
    pcontext => context
    f => saturation_difference
    call finder%init(f, context = pcontext)
    call finder%find()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call assert_equals(0, finder%err, "Saturation line error")
       call context%PT(finder%root, P, T)
       call assert_equals(expected_temperature, T, &
            finder%root_tolerance, "Saturation line temperature")
       call assert_true(finder%iterations <= expected_iterations, &
            "Saturation line iterations")
    end if

    call finder%destroy()
    call thermo%destroy()

  contains

    PetscReal function saturation_difference(x, context) result(dp)
      ! Returns pressure difference between point 0 <= x <= 1 along
      ! line between (P0, T0) and (P1, T1) and saturation curve.
      PetscReal, intent(in) :: x
      class(*), pointer, intent(in out) :: context
      ! Locals:
      PetscReal :: P, T, Ps
      PetscInt :: err
      select type (context)
      type is (saturation_context_type)
         call context%PT(x, P, T)
         call thermo%saturation%pressure(T, Ps, err)
      end select
      dp = Ps - P
    end function saturation_difference

  end subroutine test_root_finder_saturation

!------------------------------------------------------------------------

end module root_finder_test
