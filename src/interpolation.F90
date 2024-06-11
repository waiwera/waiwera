!   Copyright 2016 University of Auckland.

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

module interpolation_module
  !! Module for interpolation tables.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module

  implicit none
  private

  PetscInt, parameter, public :: INTERP_LINEAR = 0, INTERP_STEP = 1, INTERP_PCHIP = 2
  PetscInt, parameter, public :: INTERP_AVERAGING_ENDPOINT = 0, &
       INTERP_AVERAGING_INTEGRATE = 1
  PetscInt, parameter, public :: max_interpolation_str_length = 16
  character(max_interpolation_str_length), parameter, public :: &
       default_interpolation_str = "linear"
  PetscInt, parameter, public :: max_averaging_str_length = 16
  character(max_averaging_str_length), parameter, public :: &
       default_averaging_str = "integrate"

  type, public :: interpolation_coordinate_type
     !! Coordinate axis for interpolation.
     private
     PetscReal, allocatable, public :: val(:) !! Coordinate values
     PetscInt, public :: size !! Number of coordinate values
     PetscInt, public :: index !! Current index along the coordinate axis
   contains
     procedure, public :: init => interpolation_coordinate_init
     procedure, public :: destroy => interpolation_coordinate_destroy
     procedure, public :: find => interpolation_coordinate_find
  end type interpolation_coordinate_type

  type, public :: interpolation_table_type
     !! Interpolation table type, with default linear interpolation.
     private
     type(interpolation_coordinate_type), public :: coord !! Data coordinates
     PetscReal, allocatable, public :: val(:,:) !! Data values
     PetscInt, public :: dim !! Dimension of arrays being interpolated
     PetscInt, public :: interpolation_type !! Interpolation type
     PetscInt, public :: averaging_type !! Averaging type
     PetscBool, public :: continuous !! Whether interpolant is continuous or not
     procedure(averaging_function), pointer, public :: average_array_internal
   contains
     private
     procedure :: interpolation_table_init
     procedure :: interpolation_table_init_default_averaging
     generic, public :: init => interpolation_table_init, &
          interpolation_table_init_default_averaging
     procedure, public :: interpolant => interpolation_table_linear_interpolant
     procedure, public :: inverse_interpolant => interpolation_table_linear_inverse_interpolant
     procedure, public :: set_averaging_type => interpolation_table_set_averaging_type
     procedure, public :: destroy => interpolation_table_destroy
     procedure, public :: interpolate_at_index => interpolation_table_interpolate_at_index
     procedure, public :: interpolate_component_at_index => &
          interpolation_table_interpolate_component_at_index
     procedure, public :: interpolation_table_interpolate
     procedure, public :: interpolation_table_interpolate_component
     generic, public :: interpolate => interpolation_table_interpolate, &
          interpolation_table_interpolate_component
     procedure, public :: find_component_at_index => &
          interpolation_table_find_component_at_index
     procedure :: interpolation_table_average_endpoint
     procedure :: interpolation_table_average_integrate
     procedure :: average_array => interpolation_table_average_array
     procedure :: average_component => interpolation_table_average_component
     generic, public :: average => average_array, average_component
  end type interpolation_table_type

  type, public :: pinterpolation_table_type
     !! Pointer to interpolation table.
     class(interpolation_table_type), pointer, public :: ptr
  end type pinterpolation_table_type

  type, public, extends(interpolation_table_type) :: interpolation_table_step_type
     !! Step (piecewise constant) interpolation table
   contains
     procedure, public :: interpolation_table_init => interpolation_table_step_init
     procedure, public :: interpolant => interpolation_table_step_interpolant
     procedure, public :: inverse_interpolant => interpolation_table_step_inverse_interpolant
  end type interpolation_table_step_type

  type, public, extends(interpolation_table_type) :: interpolation_table_pchip_type
     !! PCHIP (Piecewise Cubic Hermite Interpolation Polynomial) interpolation table
     private
     PetscReal, allocatable :: deriv(:,:) !! Derivatives array for each dimension
   contains
     procedure :: get_derivatives => interpolation_table_pchip_get_derivatives
     procedure, public :: interpolation_table_init => interpolation_table_pchip_init
     procedure, public :: interpolant => interpolation_table_pchip_interpolant
     procedure, public :: inverse_interpolant => interpolation_table_pchip_inverse_interpolant
     procedure, public :: destroy => interpolation_table_pchip_destroy
  end type interpolation_table_pchip_type

  interface

     function averaging_function(self, interval)
       !! Function for averaging data value array over a given x
       !! interval.
       import :: interpolation_table_type
       class(interpolation_table_type), intent(in out) :: self
       PetscReal, intent(in) :: interval(2) !! x interval
       PetscReal :: averaging_function(self%dim)
     end function averaging_function

  end interface

  public :: interpolation_type_from_str, averaging_type_from_str

contains

!------------------------------------------------------------------------

  PetscInt function interpolation_type_from_str(str) result(interpolation_type)
    !! Returns interpolation type from string specification.

    use utils_module, only: str_to_lower

    character(len = *), intent(in) :: str
    select case (str_to_lower(str))
    case ("linear")
       interpolation_type = INTERP_LINEAR
    case ("step")
       interpolation_type = INTERP_STEP
    case ("pchip")
       interpolation_type = INTERP_PCHIP
    end select

  end function interpolation_type_from_str

!------------------------------------------------------------------------

  PetscInt function averaging_type_from_str(str) result(averaging_type)
    !! Returns averaging type from string specification.

    use utils_module, only: str_to_lower

    character(len = *), intent(in) :: str
    select case (str_to_lower(str))
    case ("endpoint")
       averaging_type = INTERP_AVERAGING_ENDPOINT
    case ("integrate")
       averaging_type = INTERP_AVERAGING_INTEGRATE
    end select

  end function averaging_type_from_str

!------------------------------------------------------------------------
! interpolation_coordinate_type:
!------------------------------------------------------------------------

  subroutine interpolation_coordinate_init(self, values)
    !! Initialises coordinate axis.

    class(interpolation_coordinate_type), intent(in out) :: self
    PetscReal, intent(in) :: values(:)

    self%val = values
    self%size = size(values, 1)
    self%index = 1

  end subroutine interpolation_coordinate_init

!------------------------------------------------------------------------

  subroutine interpolation_coordinate_destroy(self)
    !! Destroys coordinate axis.
    class(interpolation_coordinate_type), intent(in out) :: self

    if (allocated(self%val)) deallocate(self%val)

  end subroutine interpolation_coordinate_destroy

!------------------------------------------------------------------------

  subroutine interpolation_coordinate_find(self, x)
    !! Updates index property so that val(index) <= x < val(index + 1).
    !! If x is below the lower coordinate limit, index = 0.
    !! If x is above the upper coordinate limit, index = size.
    !! Search algorithm adapted from Press, Teukolsky, Vetterling and Flannery,
    !! "Numerical Recipes in Fortran", 2nd ed., 1992.

    class(interpolation_coordinate_type), intent(in out) :: self
    PetscReal, intent(in) :: x !! x value to interpolate at
    ! Locals:
    PetscInt :: i1, i2

    if (x <= self%val(1)) then
       self%index = 0
    else if (x >= self%val(self%size)) then
       self%index = self%size
    else
       call bracket(x, i1, i2)
       call bisect(x, i1, i2)
       self%index = i1
    end if

  contains

!........................................................................

    subroutine bracket(x, i1, i2)
      !! Find indices i1, i2 which bracket value x, using self%index
      !! as an initial estimate of i1 if possible.

      PetscReal, intent(in) :: x
      PetscInt, intent(out) :: i1, i2
      ! Locals:
      PetscInt :: inc
      PetscBool :: found


      if ((self%index >= 1) .and. (self%index <= self%size)) then

         i1 = self%index
         inc = 1
         found = PETSC_FALSE
         if (x >= self%val(i1)) then

            ! Search up:
            do while (.not. found)
               i2 = i1 + inc
               if (i2 > self%size) then
                  i2 = self%size + 1
                  found = PETSC_TRUE
               else if (x >= self%val(i2)) then
                  i1 = i2
                  inc = inc + inc
               else
                  found = PETSC_TRUE
               end if
            end do
         else

            ! Search down:
            i2 = i1
            do while (.not. found)
               i1 = i2 - inc
               if (i1 < 1) then
                  i1 = 0
                  found = PETSC_TRUE
               else if (x < self%val(i1)) then
                  i2 = i1
                  inc = inc + inc
               else
                  found = PETSC_TRUE
               end if
            end do
         end if

      else
         i1 = 0
         i2 = self%size + 1
      end if

    end subroutine bracket

!........................................................................

    subroutine bisect(x, i1, i2)
      !! Do bisection to narrow initial bracketing interval i1, i2 so
      !! that i2 = i1 + 1.

      PetscReal, intent(in) :: x
      PetscInt, intent(in out) :: i1, i2
      ! Locals:
      PetscInt :: im

      do while (i2 - i1 > 1)
         im = (i1 + i2) / 2
         if (x >= self%val(im)) then
            i1 = im
         else
            i2 = im
         end if
      end do

    end subroutine bisect

  end subroutine interpolation_coordinate_find

!------------------------------------------------------------------------
! interpolation_table_type:  
!------------------------------------------------------------------------
  
  subroutine interpolation_table_init(self, array, averaging_type)
    !! Initialises interpolation table.  The values of array(:, 1) are
    !! the x coordinate values to interpolate between, while
    !! array(:, 2:) are the corresponding data values. The array is
    !! assumed to have size (at least) 2 in the second dimension.

    use utils_module, only: array_sorted

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: array(:,:)
    PetscInt, intent(in) :: averaging_type
    ! Locals:
    PetscInt :: n(2), i
    PetscInt, allocatable :: isort(:)
    PetscReal, allocatable :: sorted_array(:, :)
    PetscErrorCode :: ierr

    n = shape(array)
    self%dim = n(2) - 1

    ! Sort data if necessary:
    if (.not. array_sorted(array(:, 1))) then
       isort = [(i - 1, i = 1, n(1))]
       call PetscSortRealWithPermutation(n(1), array(:, 1), &
            isort, ierr); CHKERRQ(ierr)
       allocate(sorted_array, mold = array)
       do i = 1, n(1)
          sorted_array(i, :) = array(isort(i) + 1, :)
       end do
       deallocate(isort)
    else
       sorted_array = array
    end if

    call self%coord%init(sorted_array(:, 1))
    allocate(self%val(self%dim, n(1)))
    self%val = transpose(sorted_array(:, 2:))
    deallocate(sorted_array)
    self%interpolation_type = INTERP_LINEAR
    call self%set_averaging_type(averaging_type)
    self%continuous = PETSC_TRUE

  end subroutine interpolation_table_init

  subroutine interpolation_table_init_default_averaging(self, array)
    !! Initialises with default endpoint averaging.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: array(:,:)

    call self%init(array, INTERP_AVERAGING_ENDPOINT)

  end subroutine interpolation_table_init_default_averaging

!------------------------------------------------------------------------

  subroutine interpolation_table_set_averaging_type(self, averaging_type)
    !! Sets averaging function for the table.

    class(interpolation_table_type), intent(in out) :: self
    PetscInt, intent(in) :: averaging_type

    self%averaging_type = averaging_type

    select case (averaging_type)
    case (INTERP_AVERAGING_ENDPOINT)
       self%average_array_internal => interpolation_table_average_endpoint
    case (INTERP_AVERAGING_INTEGRATE)
       self%average_array_internal => interpolation_table_average_integrate
    case default
       self%average_array_internal => interpolation_table_average_endpoint
    end select

  end subroutine interpolation_table_set_averaging_type

!------------------------------------------------------------------------

  function interpolation_table_linear_interpolant(self, x) result(y)
    !! Piecewise linear interpolant.

    class(interpolation_table_type), intent(in) :: self
    PetscReal, intent(in) :: x
    PetscReal :: y(self%dim)
    ! Locals:
    PetscReal :: xi

    associate(xv => self%coord%val(self%coord%index: self%coord%index + 1), &
         v => self%val(:, self%coord%index: self%coord%index + 1))
      xi = (x - xv(1)) / (xv(2) - xv(1))
      y = (1._dp - xi) * v(:, 1) + xi * v(:, 2)
    end associate

  end function interpolation_table_linear_interpolant

!------------------------------------------------------------------------

  subroutine interpolation_table_linear_inverse_interpolant(self, y, &
       component, x, err)
    !! Piecewise linear inverse interpolant.

    class(interpolation_table_type), intent(in) :: self
    PetscReal, intent(in) :: y
    PetscInt, intent(in) :: component
    PetscReal, intent(out) :: x
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: xi, vmax, vs(2), ys
    PetscReal, parameter :: tol = 1.e-8_dp

    err = 0

    associate(xv => self%coord%val(self%coord%index: self%coord%index + 1), &
         v => self%val(component, self%coord%index : self%coord%index + 1))
      vmax = max(abs(v(1)), abs(v(2)))
      if (abs(v(2) - v(1)) >= tol * vmax) then
         vs = v / vmax
         ys = y / vmax
         xi = (ys - vs(1)) / (vs(2) - vs(1))
         x = (1._dp - xi) * xv(1) + xi * xv(2)
      else
         err = 1
      end if
    end associate

  end subroutine interpolation_table_linear_inverse_interpolant

!------------------------------------------------------------------------

  subroutine interpolation_table_destroy(self)
    !! Destroys interpolation table.

    class(interpolation_table_type), intent(in out) :: self

    call self%coord%destroy()
    if (allocated(self%val)) deallocate(self%val)

  end subroutine interpolation_table_destroy

!------------------------------------------------------------------------

  function interpolation_table_interpolate_at_index(self, x) result(y)
    !! Returns interpolated y value for the given x, using the current
    !! coordinate index.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: x !! x value to interpolate at
    PetscReal :: y(self%dim)

    if (self%coord%index <= 0) then
       y = self%val(:, 1)
    else if (self%coord%index >= self%coord%size) then
       y = self%val(:, self%coord%size)
    else
       y = self%interpolant(x)
    end if

  end function interpolation_table_interpolate_at_index

!------------------------------------------------------------------------

  function interpolation_table_interpolate_component_at_index( &
       self, x, component) result(yi)
    !! Returns interpolated y value component for the given x, using
    !! the current coordinate index.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: x !! x value to interpolate at
    PetscInt, intent(in) :: component !! Component to interpolate
    PetscReal :: yi
    ! Locals:
    PetscReal :: y(self%dim)

    y = self%interpolate_at_index(x)
    yi = y(component)

  end function interpolation_table_interpolate_component_at_index

!------------------------------------------------------------------------

  function interpolation_table_interpolate(self, x) result(y)
    !! Returns interpolated y value for the given x.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: x !! x value to interpolate at
    PetscReal :: y(self%dim)

    call self%coord%find(x)
    y = self%interpolate_at_index(x)

  end function interpolation_table_interpolate

!------------------------------------------------------------------------

  function interpolation_table_interpolate_component(self, x, &
       component) result(yi)
    !! Returns interpolated y value for the given x.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: x !! x value to interpolate at
    PetscInt, intent(in) :: component !! Component to interpolate
    PetscReal :: yi

    call self%coord%find(x)
    yi = self%interpolate_component_at_index(x, component)

  end function interpolation_table_interpolate_component

!------------------------------------------------------------------------

  subroutine interpolation_table_find_component_at_index(self, &
       yi, component, x, err)
    !! At current index, finds x value corresponding to the given
    !! value yi of the specified component of y.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: yi
    PetscInt, intent(in) :: component
    PetscReal, intent(out) :: x
    PetscErrorCode, intent(out) :: err

    if ((self%coord%index <= 0) .or. &
         (self%coord%index >= self%coord%size)) then
       err = 1
    else
       call self%inverse_interpolant(yi, component, x, err)
    end if

  end subroutine interpolation_table_find_component_at_index

!------------------------------------------------------------------------

  function interpolation_table_average_endpoint(self, interval) &
       result(y)
    !! Returns y value averaged over the specified x interval. Values
    !! at the end points are interpolated first, then these two values are
    !! averaged to give the result.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: interval(2) !! x interval
    PetscReal :: y(self%dim)
    ! Locals:
    PetscReal :: yend1(self%dim), yend2(self%dim)

    yend1 = self%interpolate(interval(1))
    yend2 = self%interpolate(interval(2))

    y = 0.5_dp * (yend1 + yend2)

  end function interpolation_table_average_endpoint

!------------------------------------------------------------------------

  function interpolation_table_average_integrate(self, interval) &
       result(y)
    !! Returns y value averaged over the specified x interval,
    !! calculated by integration. Values of y are assumed to vary
    !! between data points according to the table's interpolation type.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: interval(2) !! x interval
    PetscReal :: y(self%dim)
    ! Locals:
    PetscReal :: integral(self%dim), dx, xav
    PetscReal :: x1, x2, y1(self%dim), y2(self%dim)
    PetscInt :: istart, i
    PetscBool :: finished
    PetscReal, parameter :: tol = 1.e-15_dp

    dx = interval(2) - interval(1)
    if (dx >= tol) then

       integral = 0._dp
       x1 = interval(1)
       y1 = self%interpolate(x1)
       istart = self%coord%index
       finished = PETSC_FALSE

       do i = istart, self%coord%size - 1
          x2 = self%coord%val(i + 1)
          if (x2 > interval(2)) then
             x2 = interval(2)
             finished = PETSC_TRUE
          end if
          self%coord%index = i
          y2 = self%interpolate_at_index(x2)
          call update_integral(x1, x2, y1, y2, integral)
          x1 = x2
          if (self%continuous) then
             y1 = y2
          else
             self%coord%index = i + 1
             y1 = self%interpolate_at_index(x1)
          end if
          if (finished) exit
       end do

       if (x1 < interval(2)) then
          ! data ran out before end of interval:
          x2 = interval(2)
          y2 = self%interpolate(x2)
          call update_integral(x1, x2, y1, y2, integral)
       end if

       y = integral / dx

    else
       ! Interval close to zero- just interpolate at midpoint:
       xav = 0.5_dp *(interval(1) + interval(2))
       y = self%interpolate(xav)
    end if

  contains

    subroutine update_integral(x1, x2, y1, y2, integral)
      !! Add contribution from last data interval to integral.
      PetscReal, intent(in) :: x1, x2, y1(:), y2(:)
      PetscReal, intent(in out) :: integral(:)
      integral = integral + (x2 - x1) * 0.5_dp * (y1 + y2)
    end subroutine update_integral

  end function interpolation_table_average_integrate

!------------------------------------------------------------------------

  function interpolation_table_average_array(self, interval) &
       result(y)
    !! Returns y value averaged over the specified x interval.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: interval(2) !! x interval
    PetscReal :: y(self%dim)

    y = self%average_array_internal(interval)

  end function interpolation_table_average_array

!------------------------------------------------------------------------

  function interpolation_table_average_component(self, interval, index) &
       result(yi)
    !! Returns specified component of y value averaged over the
    !! specified x interval.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: interval(2) !! x interval
    PetscInt, intent(in) :: index !! Component to average
    PetscReal :: yi
    ! Locals:
    PetscReal :: y(self%dim)

    y = self%average_array(interval)
    yi = y(index)

  end function interpolation_table_average_component

!------------------------------------------------------------------------
! Step interpolation
!------------------------------------------------------------------------

  subroutine interpolation_table_step_init(self, array, averaging_type)
    !! Initialises piecewise constant interpolation table.

    class(interpolation_table_step_type), intent(in out) :: self
    PetscReal, intent(in) :: array(:,:)
    PetscInt, intent(in) :: averaging_type

    call self%interpolation_table_type%init(array, averaging_type)
    self%interpolation_type = INTERP_STEP
    self%continuous = PETSC_FALSE

  end subroutine interpolation_table_step_init

!------------------------------------------------------------------------

  function interpolation_table_step_interpolant(self, x) result(y)
    !! Piecewise constant interpolant.

    class(interpolation_table_step_type), intent(in) :: self
    PetscReal, intent(in) :: x
    PetscReal :: y(self%dim)

    y = self%val(:, self%coord%index)

  end function interpolation_table_step_interpolant

!------------------------------------------------------------------------

  subroutine interpolation_table_step_inverse_interpolant(self, y, &
       component, x, err)
    !! Piecewise constant inverse interpolant (undefined).

    class(interpolation_table_step_type), intent(in) :: self
    PetscReal, intent(in) :: y
    PetscInt, intent(in) :: component
    PetscReal, intent(out) :: x
    PetscErrorCode, intent(out) :: err

    x = 0._dp
    err = 1

  end subroutine interpolation_table_step_inverse_interpolant

!------------------------------------------------------------------------
! PCHIP (Piecewise Cubic Hermite Interpolation Polynomial) interpolation
!------------------------------------------------------------------------

  subroutine interpolation_table_pchip_get_derivatives(self)
    !! Gets PCHIP derivatives arrays for each data dimension.

    ! Adapted from routines originally in the public-domain SLATEC
    ! library (https://www.netlib.org/slatec/pchip/) and translated
    ! into modern Fortran in:

    ! PCHIP: Piecewise Cubic Hermite Interpolation Package
    ! https://github.com/jacobwilliams/PCHIP

    ! Copyright (c) 2019-2024, Jacob Williams
    ! All rights reserved.

    ! Redistribution and use in source and binary forms, with or without modification,
    ! are permitted provided that the following conditions are met:

    ! * Redistributions of source code must retain the above copyright notice, this
    !   list of conditions and the following disclaimer.

    ! * Redistributions in binary form must reproduce the above copyright notice, this
    !   list of conditions and the following disclaimer in the documentation and/or
    !   other materials provided with the distribution.

    ! * The names of its contributors may not be used to endorse or promote products
    !   derived from this software without specific prior written permission.

    ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    ! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    ! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
    ! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    ! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    ! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
    ! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    ! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    ! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    class(interpolation_table_pchip_type), intent(in out) :: self
    ! Locals:
    PetscInt :: c

    allocate(self%deriv(self%dim, self%coord%size))

    do c = 1, self%dim
       call pchip_deriv(self%coord%val, self%val(c, :), self%deriv(c, :))
    end do

  contains

!........................................................................

    subroutine pchip_deriv(x, f, d)

      use utils_module, only: sign_test

      PetscReal, intent(in) :: x(:), f(:)
      PetscReal, intent(out) :: d(:)
      ! Locals:
      integer :: n, i, s
      PetscReal :: h1, h2, del1, del2, drat1, drat2
      PetscReal :: hsum, w1, w2, dmin, dmax

      n = size(x)

      if (n == 1) then
         d = 0._dp ! constant
      else
         h1 = x(2) - x(1)
         del1 = (f(2) - f(1)) / h1

         if (n == 2) then
            d = del1 ! linear
         else
            ! n >= 3:
            h2 = x(3) - x(2)
            del2 = (f(3) - f(2)) / h2

            hsum = h1 + h2
            w1 = (h1 + hsum) / hsum
            w2 = -h1 / hsum
            d(1) = w1 * del1 + w2 * del2
            if (sign_test(d(1), del1) <= 0) then
               d(1) = 0._dp
            else if (sign_test(del1, del2) < 0) then
               dmax = 3._dp * del1
               if (abs(d(1)) > abs(dmax))  d(1) = dmax
            end if

            do i = 2, n - 1

               if (i > 2) then
                  h1 = h2
                  h2 = x(i + 1) - x(i)
                  hsum = h1 + h2
                  del1 = del2
                  del2 = (f(i+1) - f(i)) / h2
               end if

               s = sign_test(del1, del2)
               if (s > 0) then
                  w1 = (hsum + h1) / (3._dp * hsum)
                  w2 = (hsum + h2) / (3._dp * hsum)
                  dmax = max(abs(del1), abs(del2))
                  dmin = min(abs(del1), abs(del2))
                  drat1 = del1 / dmax
                  drat2 = del2 / dmax
                  d(i) = dmin / (w1 * drat1 + w2 * drat2)
               else
                  d(i) = 0._dp
               end if

            end do

            w1 = -h2 / hsum
            w2 = (h2 + hsum) / hsum
            d(n) = w1 * del1 + w2 * del2
            if (sign_test(d(n), del2) <= 0) then
               d(n) = 0._dp
            else if (sign_test(del1, del2) < 0) then
               dmax = 3._dp * del2
               if (abs(d(n)) > abs(dmax)) d(n) = dmax
            end if

         end if
      end if

    end subroutine pchip_deriv

  end subroutine interpolation_table_pchip_get_derivatives

!------------------------------------------------------------------------

  subroutine interpolation_table_pchip_init(self, array, averaging_type)
    !! Initialises PCHIP interpolation table.

    class(interpolation_table_pchip_type), intent(in out) :: self
    PetscReal, intent(in) :: array(:,:)
    PetscInt, intent(in) :: averaging_type

    call self%interpolation_table_type%init(array, averaging_type)
    self%interpolation_type = INTERP_PCHIP
    call self%get_derivatives()

  end subroutine interpolation_table_pchip_init

!------------------------------------------------------------------------

  function interpolation_table_pchip_interpolant(self, x) result(y)
    !! PCHIP interpolant.

    class(interpolation_table_pchip_type), intent(in) :: self
    PetscReal, intent(in) :: x
    PetscReal :: y(self%dim)
    ! Locals:
    PetscInt :: c
    PetscReal :: h, dx, delta, del1, del2, c2, c3

    associate(xv => self%coord%val(self%coord%index: self%coord%index + 1))

      h = xv(2) - xv(1)
      dx = x - xv(1)

      do c = 1, self%dim
         associate(v => self%val(c, self%coord%index: self%coord%index + 1), &
              d => self%deriv(c, self%coord%index: self%coord%index + 1))
           delta = (v(2) - v(1)) / h
           del1 = (d(1) - delta) / h
           del2 = (d(2) - delta)  / h
           c2 = -(2._dp * del1 + del2)
           c3 = (del1 + del2) / h
           y(c) = v(1) + dx * (d(1) + dx * (c2 + dx * c3))
         end associate
      end do

    end associate

  end function interpolation_table_pchip_interpolant

!------------------------------------------------------------------------

  subroutine interpolation_table_pchip_inverse_interpolant(self, y, &
       component, x, err)
    !! PCHIP inverse interpolant (undefined).

    class(interpolation_table_pchip_type), intent(in) :: self
    PetscReal, intent(in) :: y
    PetscInt, intent(in) :: component
    PetscReal, intent(out) :: x
    PetscErrorCode, intent(out) :: err

    x = 0._dp
    err = 1

  end subroutine interpolation_table_pchip_inverse_interpolant

!------------------------------------------------------------------------

  subroutine interpolation_table_pchip_destroy(self)
    !! Destroys PCHIP interpolation table.

    class(interpolation_table_pchip_type), intent(in out) :: self

    call self%interpolation_table_type%destroy()
    deallocate(self%deriv)

  end subroutine interpolation_table_pchip_destroy

!------------------------------------------------------------------------
  
end module interpolation_module
