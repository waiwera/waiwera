module interpolation_module
  !! Module for interpolation tables.

  use kinds_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>
#include <petsc/finclude/petscsys.h>

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
     PetscReal, allocatable, public :: val(:) !! Data values
     procedure(interpolation_function), pointer, nopass, public :: interpolant
     procedure(averaging_function), pointer, public :: average
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
     procedure :: interpolate_at_index => interpolation_table_interpolate_at_index
     procedure, public :: interpolate => interpolation_table_interpolate
     procedure :: interpolation_table_average_endpoint
     procedure :: interpolation_table_average_integrate
  end type interpolation_table_type

  interface

     PetscReal function interpolation_function(coord, val, x, index)
       !! Function for interpolating value array at a given coordinate
       !! value x between index and index + 1.
       PetscReal, intent(in) :: coord(:), val(:)
       PetscReal, intent(in) :: x
       PetscInt, intent(in) :: index
     end function interpolation_function

     PetscReal function averaging_function(self, interval)
       !! Function for averaging data value array over a given x
       !! interval.
       import :: interpolation_table_type
       class(interpolation_table_type), intent(in out) :: self
       PetscReal, intent(in) :: interval(2) !! x interval
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

  PetscReal function interpolant_linear(coord, val, x, index) result(y)
    !! Linear interpolation function.

    PetscReal, intent(in) :: coord(:), val(:)
    PetscReal, intent(in) :: x
    PetscInt, intent(in) :: index
    ! Locals:
    PetscReal :: theta

    theta = (x - coord(index)) / (coord(index + 1) - coord(index))
    y = (1._dp - theta) * val(index) + theta * val(index + 1)

  end function interpolant_linear

!------------------------------------------------------------------------

  PetscReal function interpolant_step(coord, val, x, index) result(y)
    !! Piecewise constant step interpolation function, with constant
    !! value set to the data value at index.

    PetscReal, intent(in) :: coord(:), val(:)
    PetscReal, intent(in) :: x
    PetscInt, intent(in) :: index

    y = val(index)

  end function interpolant_step

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

    deallocate(self%val)

  end subroutine interpolation_coordinate_destroy

!------------------------------------------------------------------------

  subroutine interpolation_coordinate_find(self, x)
    !! Updates index property so that val(index) <= x < val(index + 1).
    !! If x is below the lower coordinate limit, index = 0.
    !! If x is above the upper coordinate limit, index = size.

    class(interpolation_coordinate_type), intent(in out) :: self
    PetscReal, intent(in) :: x !! x value to interpolate at
    ! Locals:
    PetscInt :: i, end_index, direction

    if (x <= self%val(1)) then
       ! Below lower coordinate limit:
       self%index = 0
    else if (x >= self%val(self%size)) then
       ! Above upper table limit:
       self%index = self%size
    else

       ! Check starting index:
       if (self%index <= 0) then
          self%index = 1
       else if (self%index >= self%size) then
          self%index = self%size - 1
       end if

       ! Determine search direction:
       if (x < self%val(self%index)) then
          end_index = 1
          direction = -1
       else
          end_index = self%size - 1
          direction = 1
       end if

       do i = self%index, end_index, direction
          if ((self%val(i) <= x) .and. (x < self%val(i + 1))) then
             self%index = i
             exit
          end if
       end do

    end if

  end subroutine interpolation_coordinate_find

!------------------------------------------------------------------------
! interpolation_table_type:  
!------------------------------------------------------------------------
  
  subroutine interpolation_table_init(self, array, interpolation_type, &
       averaging_type)
    !! Initialises interpolation table.  The values of array(:,1) are
    !! the sorted x coordinate values to interpolate between, while
    !! array(:,2) are the corresponding data values. The array is
    !! assumed to have size (at least) 2 in the second dimension.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: array(:,:)
    PetscInt, intent(in) :: interpolation_type
    PetscInt, intent(in) :: averaging_type

    call self%coord%init(array(:, 1))
    self%val = array(:, 2)
    call self%set_interpolation_type(interpolation_type)
    call self%set_averaging_type(averaging_type)

  end subroutine interpolation_table_init

  subroutine interpolation_table_init_default_averaging(self, array, &
       interpolation_type)
    !! Initialises with default endpoint averaging.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: array(:,:)
    PetscInt, intent(in) :: interpolation_type

    call self%init(array, interpolation_type, INTERP_AVERAGING_ENDPOINT)

  end subroutine interpolation_table_init_default_averaging

  subroutine interpolation_table_init_default_all(self, array)
    !! Initialises with default linear interpolation and endpoint averaging.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: array(:,:)

    call self%init(array, INTERP_LINEAR, INTERP_AVERAGING_ENDPOINT)

  end subroutine interpolation_table_init_default_all

!------------------------------------------------------------------------

  subroutine interpolation_table_set_interpolation_type(self, interpolation_type)
    !! Sets interpolation function for the table.

    class(interpolation_table_type), intent(in out) :: self
    PetscInt, intent(in) :: interpolation_type

    select case (interpolation_type)
    case (INTERP_LINEAR)
       self%interpolant => interpolant_linear
    case (INTERP_STEP)
       self%interpolant => interpolant_step
    case default
       self%interpolant => interpolant_linear
    end select

  end subroutine interpolation_table_set_interpolation_type

!------------------------------------------------------------------------

  subroutine interpolation_table_set_averaging_type(self, averaging_type)
    !! Sets averaging function for the table.

    class(interpolation_table_type), intent(in out) :: self
    PetscInt, intent(in) :: averaging_type

    select case (averaging_type)
    case (INTERP_AVERAGING_ENDPOINT)
       self%average => interpolation_table_average_endpoint
    case (INTERP_AVERAGING_INTEGRATE)
       self%average => interpolation_table_average_integrate
    case default
       self%average => interpolation_table_average_endpoint
    end select

  end subroutine interpolation_table_set_averaging_type

!------------------------------------------------------------------------

  subroutine interpolation_table_destroy(self)
    !! Destroys interpolation table.

    class(interpolation_table_type), intent(in out) :: self

    call self%coord%destroy()
    deallocate(self%val)

  end subroutine interpolation_table_destroy

!------------------------------------------------------------------------

  PetscReal function interpolation_table_interpolate_at_index(self, x) result(y)
    !! Returns interpolated y value for the given x, using the current
    !! coordinate index.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: x !! x value to interpolate at

    if (self%coord%index <= 0) then
       y = self%val(1)
    else if (self%coord%index >= self%coord%size) then
       y = self%val(self%coord%size)
    else
       y = self%interpolant(self%coord%val, self%val, x, self%coord%index)
    end if

  end function interpolation_table_interpolate_at_index

!------------------------------------------------------------------------

  PetscReal function interpolation_table_interpolate(self, x) result(y)
    !! Returns interpolated y value for the given x.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: x !! x value to interpolate at

    call self%coord%find(x)
    y = self%interpolate_at_index(x)

  end function interpolation_table_interpolate

!------------------------------------------------------------------------

  PetscReal function interpolation_table_average_endpoint(self, interval) &
       result(y)
    !! Returns y value averaged over the specified x interval. Values
    !! at the end points are interpolated first, then these two values are
    !! averaged to give the result.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: interval(2) !! x interval
    ! Locals:
    PetscReal :: yend(2)

    yend(1) = self%interpolate(interval(1))
    yend(2) = self%interpolate(interval(2))

    y = 0.5_dp * sum(yend)

  end function interpolation_table_average_endpoint

!------------------------------------------------------------------------

  PetscReal function interpolation_table_average_integrate(self, interval) &
       result(y)
    !! Returns y value averaged over the specified x interval,
    !! calculated by integration. Values of y are assumed to vary
    !! between data points according to the table's interpolation type.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: interval(2) !! x interval
    ! Locals:
    PetscReal :: integral, dx, xav
    PetscReal :: x1, x2, y1, y2
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
      PetscReal, intent(in) :: x1, x2, y1, y2
      PetscReal, intent(in out) :: integral
      integral = integral + (x2 - x1) * 0.5_dp * (y1 + y2)
    end subroutine update_integral

  end function interpolation_table_average_integrate

!------------------------------------------------------------------------

end module interpolation_module
