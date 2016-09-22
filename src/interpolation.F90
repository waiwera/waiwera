module interpolation_module
  !! Module for interpolation tables.

  use kinds_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>
#include <petsc/finclude/petscsys.h>

  PetscInt, parameter, public :: INTERP_LINEAR = 0, INTERP_STEP = 1, &
       INTERP_STEP_AVERAGE = 2

  type, public :: interpolation_table_type
     !! Interpolation table type.
     private
     PetscReal, allocatable, public :: data(:,:) !! Data array of (x,y) values
     PetscInt :: size !! Number of x values in the data array
     PetscInt, public :: index !! Current index in the table
     procedure(interpolation_function), pointer, nopass, public :: interpolant
   contains
     private
     procedure :: interpolation_table_init
     procedure :: interpolation_table_init_default
     generic, public :: init => interpolation_table_init, &
          interpolation_table_init_default
     procedure, public :: set_type => interpolation_table_set_type
     procedure, public :: destroy => interpolation_table_destroy
     procedure, public :: find => interpolation_table_find
     procedure, public :: interpolate => interpolation_table_interpolate
  end type interpolation_table_type

  interface

     PetscReal function interpolation_function(data, x, index)
       !! Function for interpolating data array at a given x value
       !! between given index and index + 1.
       PetscReal, intent(in) :: data(:,:)
       PetscReal, intent(in) :: x
       PetscInt, intent(in) :: index
     end function interpolation_function

  end interface

contains

!------------------------------------------------------------------------
! Interpolation functions:
!------------------------------------------------------------------------

  PetscReal function interpolant_linear(data, x, index) result(y)
    !! Linear interpolation function.

    PetscReal, intent(in) :: data(:,:)
    PetscReal, intent(in) :: x
    PetscInt, intent(in) :: index
    ! Locals:
    PetscReal :: theta

    theta = (x - data(index, 1)) / (data(index + 1, 1) - data(index, 1))
    y = (1._dp - theta) * data(index, 2) + theta * data(index + 1, 2)

  end function interpolant_linear

!------------------------------------------------------------------------

  PetscReal function interpolant_step(data, x, index) result(y)
    !! Piecewise constant step interpolation function, with constant
    !! value set to the data value at index.

    PetscReal, intent(in) :: data(:,:)
    PetscReal, intent(in) :: x
    PetscInt, intent(in) :: index

    y = data(index, 2)

  end function interpolant_step

!------------------------------------------------------------------------

  PetscReal function interpolant_step_average(data, x, index) result(y)
    !! Piecewise constant step interpolation function, with constant
    !! value set to the average of the data values at index and index + 1.

    PetscReal, intent(in) :: data(:,:)
    PetscReal, intent(in) :: x
    PetscInt, intent(in) :: index

    y = 0.5_dp * (data(index, 2) + data(index + 1, 2))

  end function interpolant_step_average

!------------------------------------------------------------------------
! interpolation_table_type:  
!------------------------------------------------------------------------
  
  subroutine interpolation_table_init(self, array, interpolation_type)
    !! Initialises interpolation table.  The values of array(:,1) are
    !! the sorted x values to interpolate between, while array(:,2)
    !! are the corresponding y values. The array is assumed to have
    !! size (at least) 2 in the second dimension.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: array(:,:)
    PetscInt, intent(in) :: interpolation_type

    self%data = array
    self%size = size(self%data, 1)
    self%index = 1
    call self%set_type(interpolation_type)

  end subroutine interpolation_table_init

  subroutine interpolation_table_init_default(self, array)
    !! Initialises with default linear interpolation.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: array(:,:)

    call self%init(array, INTERP_LINEAR)

  end subroutine interpolation_table_init_default

!------------------------------------------------------------------------

  subroutine interpolation_table_set_type(self, interpolation_type)
    !! Sets interpolation function for the table.

    class(interpolation_table_type), intent(in out) :: self
    PetscInt, intent(in) :: interpolation_type

    select case (interpolation_type)
    case (INTERP_LINEAR)
       self%interpolant => interpolant_linear
    case (INTERP_STEP)
       self%interpolant => interpolant_step
    case (INTERP_STEP_AVERAGE)
       self%interpolant => interpolant_step_average
    case default
       self%interpolant => interpolant_linear
    end select

  end subroutine interpolation_table_set_type

!------------------------------------------------------------------------

  subroutine interpolation_table_destroy(self)
    !! Destroys interpolation table.

    class(interpolation_table_type), intent(in out) :: self

    deallocate(self%data)

  end subroutine interpolation_table_destroy

!------------------------------------------------------------------------

  subroutine interpolation_table_find(self, x)
    !! Updates index property so that data(index, 1) <= x < data(index + 1, 1).
    !! If x is below the lower table limit, index = 0.
    !! If x is above the upper table limit, index = size.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: x !! x value to interpolate at
    ! Locals:
    PetscInt :: i, end_index, direction

    if (x <= self%data(1, 1)) then
       ! Below lower table limit- use y value at lower end:
       self%index = 0
    else if (x >= self%data(self%size, 1)) then
       ! Above upper table limit- use y value at upper end:
       self%index = self%size
    else

       ! Check starting index:
       if (self%index <= 0) then
          self%index = 1
       else if (self%index >= self%size) then
          self%index = self%size - 1
       end if

       ! Determine search direction:
       if (x < self%data(self%index, 1)) then
          end_index = 1
          direction = -1
       else
          end_index = self%size - 1
          direction = 1
       end if

       do i = self%index, end_index, direction
          if ((self%data(i, 1) <= x) .and. (x < self%data(i + 1, 1))) then
             self%index = i
             exit
          end if
       end do

    end if

  end subroutine interpolation_table_find

!------------------------------------------------------------------------

  PetscReal function interpolation_table_interpolate(self, x) result(y)
    !! Returns interpolated y value for the given x.

    class(interpolation_table_type), intent(in out) :: self
    PetscReal, intent(in) :: x !! x value to interpolate at

    call self%find(x)

    if (self%index <= 0) then
       y = self%data(1, 2)
    else if (self%index >= self%size) then
       y = self%data(self%size, 2)
    else
       y = self%interpolant(self%data, x, self%index)
    end if

  end function interpolation_table_interpolate

!------------------------------------------------------------------------

end module interpolation_module
