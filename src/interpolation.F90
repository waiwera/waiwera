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

  PetscInt, parameter, public :: INTERP_LINEAR = 0, INTERP_STEP = 1
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
     !! Interpolation table type.
     private
     type(interpolation_coordinate_type), public :: coord !! Data coordinates
     PetscReal, allocatable, public :: val(:,:) !! Data values
     PetscInt, public :: dim !! Dimension of arrays being interpolated
     PetscInt, public :: interpolation_type !! Interpolation type
     PetscInt, public :: averaging_type !! Averaging type
     procedure(interpolation_function), pointer, nopass, public :: interpolant
     procedure(inverse_interpolation_function), pointer, nopass, public :: inverse_interpolant
     procedure(averaging_function), pointer, public :: average_array_internal
   contains
     private
     procedure :: interpolation_table_init
     procedure :: interpolation_table_init_default_averaging
     procedure :: interpolation_table_init_default_all
     generic, public :: init => interpolation_table_init, &
          interpolation_table_init_default_averaging, &
          interpolation_table_init_default_all
     procedure, public :: set_interpolation_type => interpolation_table_set_interpolation_type
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

  interface

     function interpolation_function(coord, val, x, dim, index)
       !! Function for interpolating value array at a given coordinate
       !! value x between index and index + 1.
       PetscReal, intent(in) :: coord(:), val(:,:)
       PetscReal, intent(in) :: x
       PetscInt, intent(in) :: dim, index
       PetscReal :: interpolation_function(dim)
     end function interpolation_function

     subroutine inverse_interpolation_function(coord, val, y, index, &
          component, x, err)
       !! Routine for inverting interpolation function.
       PetscReal, intent(in) :: coord(:), val(:,:)
       PetscReal, intent(in) :: y
       PetscInt, intent(in) :: index, component
       PetscReal, intent(out) :: x
       PetscErrorCode, intent(out) :: err
     end subroutine inverse_interpolation_function

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
! Interpolation functions:
!------------------------------------------------------------------------

  function interpolant_linear(coord, val, x, dim, index) result(y)
    !! Linear interpolation function.

    PetscReal, intent(in) :: coord(:), val(:,:)
    PetscReal, intent(in) :: x
    PetscInt, intent(in) :: dim, index
    PetscReal :: y(dim)
    ! Locals:
    PetscReal :: xi

    xi = (x - coord(index)) / (coord(index + 1) - coord(index))
    y = (1._dp - xi) * val(:, index) + xi * val(:, index + 1)

  end function interpolant_linear

!------------------------------------------------------------------------

  subroutine inverse_interpolant_linear(coord, val, y, index, &
       component, x, err)
    !! Routine for inverting linear interpolation function.
    PetscReal, intent(in) :: coord(:), val(:,:)
    PetscReal, intent(in) :: y
    PetscInt, intent(in) :: index, component
    PetscReal, intent(out) :: x
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: d, xi
    PetscReal, parameter :: tol = 1.e-8_dp

    err = 0

    d = val(component, index + 1) - val(component, index)

    if (abs(d) >= tol) then
       xi = (y - val(component, index)) / d
       x = (1._dp - xi) * coord(index) + xi * coord(index + 1)
    else
       err = 1
    end if

  end subroutine inverse_interpolant_linear

!------------------------------------------------------------------------

  function interpolant_step(coord, val, x, dim, index) result(y)
    !! Piecewise constant step interpolation function, with constant
    !! value set to the data value at index.

    PetscReal, intent(in) :: coord(:), val(:,:)
    PetscReal, intent(in) :: x
    PetscInt, intent(in) :: dim, index
    PetscReal :: y(dim)

    y = val(:, index)

  end function interpolant_step

!------------------------------------------------------------------------
! interpolation_coordinate_type:
!------------------------------------------------------------------------

  subroutine interpolation_coordinate_init(self, values, err)
    !! Initialises coordinate axis. Returns err > 0 if the specified
    !! values array is not sorted.

    use utils_module, only: array_sorted

    class(interpolation_coordinate_type), intent(in out) :: self
    PetscReal, intent(in) :: values(:)
    PetscErrorCode, intent(out) :: err

    if (array_sorted(values)) then
       self%val = values
       self%size = size(values, 1)
       self%index = 1
       err = 0
    else
       err = 1
    end if

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
  
  subroutine interpolation_table_init(self, array, interpolation_type, &
       averaging_type, err)
    !! Initialises interpolation table.  The values of array(:, 1) are
    !! the sorted x coordinate values to interpolate between, while
    !! array(:, 2:) are the corresponding data values. The array is
    !! assumed to have size (at least) 2 in the second dimension. The
    !! error code err returns > 0 if the table could not be
    !! initialised.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: array(:,:)
    PetscInt, intent(in) :: interpolation_type
    PetscInt, intent(in) :: averaging_type
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: n(2)

    call self%coord%init(array(:, 1), err)
    if (err == 0) then
       n = shape(array)
       self%dim = n(2) - 1
       allocate(self%val(self%dim, n(1)))
       self%val = transpose(array(:, 2:))
       call self%set_interpolation_type(interpolation_type)
       call self%set_averaging_type(averaging_type)
    end if

  end subroutine interpolation_table_init

  subroutine interpolation_table_init_default_averaging(self, array, &
       interpolation_type, err)
    !! Initialises with default endpoint averaging.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: array(:,:)
    PetscInt, intent(in) :: interpolation_type
    PetscErrorCode, intent(out) :: err

    call self%init(array, interpolation_type, INTERP_AVERAGING_ENDPOINT, err)

  end subroutine interpolation_table_init_default_averaging

  subroutine interpolation_table_init_default_all(self, array, err)
    !! Initialises with default linear interpolation and endpoint averaging.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: array(:,:)
    PetscErrorCode, intent(out) :: err

    call self%init(array, INTERP_LINEAR, INTERP_AVERAGING_ENDPOINT, err)

  end subroutine interpolation_table_init_default_all

!------------------------------------------------------------------------

  subroutine interpolation_table_set_interpolation_type(self, interpolation_type)
    !! Sets interpolation function for the table.

    class(interpolation_table_type), intent(in out) :: self
    PetscInt, intent(in) :: interpolation_type

    self%interpolation_type = interpolation_type

    select case (interpolation_type)
    case (INTERP_LINEAR)
       self%interpolant => interpolant_linear
       self%inverse_interpolant => inverse_interpolant_linear
    case (INTERP_STEP)
       self%interpolant => interpolant_step
       self%inverse_interpolant => null()
    case default
       self%interpolant => interpolant_linear
       self%inverse_interpolant => inverse_interpolant_linear
    end select

  end subroutine interpolation_table_set_interpolation_type

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
       y = self%interpolant(self%coord%val, self%val, x, self%dim, &
            self%coord%index)
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
       call self%inverse_interpolant(self%coord%val, self%val, &
            yi, self%coord%index, component, x, err)
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
          y1 = y2
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
  
end module interpolation_module
