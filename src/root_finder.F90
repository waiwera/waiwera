!   Copyright 2017 University of Auckland.

!   This file is part of Waiwera.

!   Waiwera is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.

!   Waiwera is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.

!   You should have received a copy of the GNU Lesser General Public License
!   along with Waiwera.  If not, see <http://www.gnu.org/licenses/>.

module root_finder_module
  !! Finding roots of a function of one variable.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module

  implicit none
  private

  PetscInt, parameter, public :: ROOT_FINDER_INTERVAL_NOT_BRACKETED = 1
  PetscInt, parameter, public :: ROOT_FINDER_ITERATIONS_EXCEEDED = 2

  type, public :: root_finder_type
     !! Root finder type.
     private
     procedure(root_finder_function), pointer, nopass, public :: f !! Function to find root of
     class(*), pointer, public :: context  !! Context data for the function f
     PetscReal, public :: interval(2) !! Bracketing interval for the root
     PetscReal, public :: root !! Root value found
     PetscReal, public :: root_tolerance !! Convergence tolerance for root value
     PetscReal, public :: function_tolerance !! Convergence tolerance for function value
     PetscInt, public :: iterations !! Number of iterations taken
     PetscInt, public :: max_iterations !! Maximum allowable number of iterations
     PetscErrorCode, public :: err !! Error flag
   contains
     private
     procedure, public :: init => root_finder_init
     procedure, public :: destroy => root_finder_destroy
     procedure, public :: find => root_finder_find
  end type root_finder_type

  interface
     PetscReal function root_finder_function(x, context)
       PetscReal, intent(in) :: x
       class(*), pointer, intent(in out) :: context
     end function root_finder_function
  end interface

  public :: root_finder_function

contains

!------------------------------------------------------------------------

  subroutine root_finder_init(self, f, interval, root_tolerance, &
       function_tolerance, max_iterations, context)
    !! Initialise root finder.
    class(root_finder_type), intent(in out) :: self

    procedure(root_finder_function), pointer, intent(in) :: f
    PetscReal, intent(in), optional :: interval(2)
    PetscReal, intent(in), optional :: root_tolerance
    PetscReal, intent(in), optional :: function_tolerance
    PetscInt, intent(in), optional :: max_iterations
    class(*), pointer, intent(in), optional :: context
    ! Locals:
    PetscReal, parameter :: default_interval(2) = [0._dp, 1._dp]
    PetscReal, parameter :: default_root_tolerance = 1.e-8_dp
    PetscReal, parameter :: default_function_tolerance = 1.e-8_dp
    PetscInt, parameter :: default_max_iterations = 100

    self%f => f
    if (present(interval)) then
       self%interval = interval
    else
       self%interval = default_interval
    end if
    if (present(root_tolerance)) then
       self%root_tolerance = root_tolerance
    else
       self%root_tolerance = default_root_tolerance
    end if
    if (present(function_tolerance)) then
       self%function_tolerance = function_tolerance
    else
       self%function_tolerance = default_function_tolerance
    end if
    if (present(max_iterations)) then
       self%max_iterations = max_iterations
    else
       self%max_iterations = default_max_iterations
    end if
    if (present(context)) then
       self%context => context
    else
       self%context => null()
    end if

    self%iterations = 0
    self%root = 0._dp
    self%err = 0

  end subroutine root_finder_init

!------------------------------------------------------------------------

  subroutine root_finder_destroy(self)
    !! Destroy root finder.
    class(root_finder_type), intent(in out) :: self

    self%f => null()
    self%context => null()

  end subroutine root_finder_destroy

!------------------------------------------------------------------------

  subroutine root_finder_find(self)
    !! Find root, using Brent's method. Algorithm based on zbrent() in
    !! Press, Teukolsky, Vetterling and Flannery, "Numerical Recipes
    !! in Fortran", 2nd ed., 1992.

    class(root_finder_type), intent(in out) :: self
    ! Locals:
    PetscReal :: a, b, c, d, e
    PetscReal :: fa, fb, fc
    PetscReal :: dx, p, pc, q, r, s
    PetscInt :: iter
    PetscBool :: found
    PetscReal, parameter :: small = 1.e-16_dp

    self%iterations = 0
    self%err = 0
    self%root = 0._dp
    found = PETSC_FALSE

    a = self%interval(1)
    b = self%interval(2)
    fa = self%f(a, self%context)
    fb = self%f(b, self%context)

    if (fa * fb > 0._dp) then
       self%err = ROOT_FINDER_INTERVAL_NOT_BRACKETED
    else

       c = b
       fc = fb

       do iter = 1, self%max_iterations

          if (fb * fc > 0._dp) then
             c = a
             fc = fa
             d = b - a
             e = d
          end if

          if (abs(fc) < abs(fb)) then
             a = b
             b = c
             c = a
             fa = fb
             fb = fc
             fc = fa
          end if

          dx = 0.5_dp * (c - b)

          if ((abs(dx) <= self%root_tolerance) .or. &
               (abs(fb) <= self%function_tolerance)) then

             ! Root found:
             found = PETSC_TRUE
             exit

          else

             if ((abs(e) >= self%root_tolerance) .and. &
                  (abs(fa) > abs(fb))) then

                ! Inverse quadratic interpolation:
                s = fb / fa
                if (abs(a - c) <= small) then
                   p = 2._dp * dx * s
                   q = 1._dp - s
                else
                   q = fa / fc
                   r = fb / fc
                   p = s * (2._dp * dx * q * (q - r) - (b - a) * (r - 1._dp))
                   q = (q - 1._dp) * (r - 1._dp) * (s - 1._dp)
                end if

                ! Check bounds:
                if (p > 0._dp) then
                   q = -q
                else
                   p = -p
                end if

                pc = min(3._dp * dx * q - abs(self%root_tolerance * q), &
                     abs(e * q))
                if (2._dp * p < pc) then
                   ! Accept interpolation:
                   e = d
                   d = p / q
                else
                   ! Interpolation failed- use bisection:
                   d = dx
                   e = d
                end if

             else
                ! Bounds decreasing too slowly- use bisection:
                   d = dx
                   e = d
             end if

             a = b
             fa = fb
             if (abs(d) > self%root_tolerance) then
                b = b + d
             else
                b = b + sign(self%root_tolerance, dx)
             end if
             fb = self%f(b, self%context)

          end if

       end do

       self%root = b
       self%iterations = iter
       if (.not. found) then
          self%err = ROOT_FINDER_ITERATIONS_EXCEEDED
       end if

    end if

  end subroutine root_finder_find

!------------------------------------------------------------------------

end module root_finder_module
