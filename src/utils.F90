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

module utils_module
  !! Utility functions for string handling, formatting, file names
  !! etc. and constants.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module

  implicit none
  private

  PetscReal, parameter, public :: pi = 4._dp * atan(1._dp)

  interface polynomial
     module procedure polynomial_single
     module procedure polynomial_multiple
  end interface polynomial

  interface array_cumulative_sum
     module procedure array_cumulative_sum_real
     module procedure array_cumulative_sum_integer
  end interface array_cumulative_sum

  interface array_sorted
     module procedure array_sorted_int
     module procedure array_sorted_real
  end interface array_sorted

  public :: str_to_lower, &
       int_str_len, str_array_index, &
       split_filename, change_filename_extension, &
       date_time_str, degrees_to_radians, rotation_matrix_2d, &
       polynomial, array_pair_sum, array_cumulative_sum, &
       array_exclusive_products, get_mpi_int_gather_array, &
       array_sorted, clock_elapsed_time
  
contains

!------------------------------------------------------------------------

  elemental function str_to_lower(a) result(b)
    !! Converts a string to all lower case.

    character(len = *), intent(in) :: a !! Input string
    character(len = len(a)) :: b !! Output lowercase string
    integer :: i,j

    b = a
    do i = 1, len(b)
       j = iachar(b(i:i))
       if (j >= iachar("A") .and. j <= iachar("Z") ) then
          b(i:i) = achar(iachar(b(i:i)) + 32)
       end if
    end do

  end function str_to_lower

!------------------------------------------------------------------------

  recursive function int_str_len(i) result (w)
    !! Returns minimum length of string needed to represent a given
    !! integer i.

    PetscInt, intent(in) :: i !! Input integer
    PetscInt :: w !! Output string length

    if (i == 0) then
       w = 1
    else if (i > 0) then
       w = 1 + int(log10(real(i)))
    else
       w = 1 + int_str_len(-i)
    end if

  end function int_str_len

!------------------------------------------------------------------------

  PetscInt function str_array_index(str, str_array) result(index)
    !! Returns index of given string in an array of strings (or -1 if
    !! it isn't in there).

    character(len = *), intent(in) :: str
    character(len = *), intent(in) :: str_array(:)
    ! Locals:
    PetscInt :: i, n

    index = -1
    n = size(str_array)
    do i = 1, n
       if (trim(str) == trim(str_array(i))) then
          index = i
          exit
       end if
    end do

  end function str_array_index

!------------------------------------------------------------------------

  subroutine split_filename(filename, base, ext)
    !! Splits filename into base and extension.

    character(*), intent(in) :: filename !! File name
    character(:), allocatable, intent(out) :: base !! Base part of filename
    character(:), allocatable, intent(out) :: ext !! File extension
    ! Locals:
    PetscInt:: i, n, base_end, ext_start
  
    n = len(filename)
    i = scan(filename, '.', PETSC_TRUE)
    if ((i > 0) .and. (i < n)) then
       base_end = i - 1
       ext_start = i + 1
       base = filename(1: base_end)
       ext = filename(ext_start: n)
    else
       base = filename
       ext = ""
    end if

  end subroutine split_filename

!------------------------------------------------------------------------

  function change_filename_extension(filename, ext) result(new_filename)
    !! Changes filename extension.

    character(*), intent(in) :: filename !! File name
    character(*), intent(in) :: ext !! New file extension
    character(:), allocatable :: new_filename !! Output file name
    ! Locals:
    character(:), allocatable :: base, oldext

    call split_filename(filename, base, oldext)
    new_filename = trim(base) // '.' // trim(ext)

    deallocate(base, oldext)

  end function change_filename_extension

!------------------------------------------------------------------------

  character(25) function date_time_str()
    !! Returns string with current date and time.

    ! Locals:
    character(8) :: datestr
    character(10) :: timestr
    character(5) :: zonestr

    call date_and_time(datestr, timestr, zonestr)

    date_time_str = datestr // ' ' // timestr // ' ' // zonestr

  end function date_time_str

!------------------------------------------------------------------------

  PetscReal function degrees_to_radians(degrees) result(radians)
    !! Converts angle from degrees to radians.

    PetscReal, intent(in) :: degrees

    radians = degrees * pi / 180._dp

  end function degrees_to_radians

!------------------------------------------------------------------------

  function rotation_matrix_2d(angle) result(M)
    !! Returns a 2x2 rotation matrix corresponding to the given angle
    !! (anti-clockwise, in radians).

    PetscReal, intent(in) :: angle
    PetscReal :: M(2, 2)
    ! Locals:
    PetscReal :: c, s

    c = cos(angle)
    s = sin(angle)
    M = reshape([c, -s, s, c], [2, 2])

  end function rotation_matrix_2d

!------------------------------------------------------------------------

  function polynomial_single(a, x) result(p)
    !! Evaluate polynomial a1 + a2*x + a3 * x^2 + ..., using Horner's
    !! method.

    PetscReal, intent(in) :: a(:)
    PetscReal, intent(in) :: x
    PetscReal :: p
    ! Locals:
    PetscInt :: i

    associate(n => size(a))
      p = a(n)
      do i = n - 1, 1, -1
         p = a(i) + x * p
      end do
    end associate

  end function polynomial_single

!------------------------------------------------------------------------

  function polynomial_multiple(a, x) result(p)
    !! Evaluate polynomials a(:, 1) + a(:, 2) * x + a(:, 3) * x^2 +
    !! ..., using Horner's method.

    PetscReal, intent(in) :: a(:, :)
    PetscReal, intent(in) :: x
    PetscReal :: p(size(a, 1))
    ! Locals:
    PetscInt :: i

    associate(n => size(a, 2))
      p = a(:, n)
      do i = n - 1, 1, -1
         p = a(:, i) + x * p
      end do
    end associate

  end function polynomial_multiple

!------------------------------------------------------------------------

  PetscReal function array_pair_sum(a) result(s)
    !! Returns sum of products of consecutive pairs in an array
    !! (including the pair formed by the last and first elements).

    PetscReal, intent(in) :: a(:)
    ! Locals:
    PetscInt :: i, i1

    s = 0._dp
    associate(n => size(a))
      if (n == 1) then
         s = a(1)
      else
         do i = 1, n
            i1 = i + 1
            if (i1 > n) i1 = i1 - n
            s = s + a(i) * a(i1)
         end do
      end if
    end associate

  end function array_pair_sum

!------------------------------------------------------------------------

  function array_cumulative_sum_real(a) result(s)
    !! Cumulative sums of a real array.

    PetscReal, intent(in) :: a(:)
    PetscReal :: s(size(a))
    ! Locals:
    PetscInt :: i

    associate(n => size(a))
      if (n > 0) then
         s(1) = a(1)
         do i = 2, n
            s(i) = s(i - 1) + a(i)
         end do
      end if
    end associate

  end function array_cumulative_sum_real

!------------------------------------------------------------------------

  function array_cumulative_sum_integer(a) result(s)
    !! Cumulative sums of an integer array.

    PetscInt, intent(in) :: a(:)
    PetscInt :: s(size(a))
    ! Locals:
    PetscInt :: i

    associate(n => size(a))
      if (n > 0) then
         s(1) = a(1)
         do i = 2, n
            s(i) = s(i - 1) + a(i)
         end do
      end if
    end associate

  end function array_cumulative_sum_integer

!------------------------------------------------------------------------

  function array_exclusive_products(a) result(p)
    !! Returns products of array elements, excluding successive
    !! elements of the array. If the array has only one element, the
    !! returned result is 1.

    PetscReal, intent(in) :: a(:)
    PetscReal :: p(size(a))
    ! Locals:
    PetscInt :: i
    PetscReal :: x(size(a))

    do i = 1, size(a)
       x = a
       x(i) = 1._dp
       p(i) = product(x)
    end do

  end function array_exclusive_products

!------------------------------------------------------------------------

  function get_mpi_int_gather_array() result(array)
    !! Returns integer array for use in MPI gather call. This is of
    !! size equal to the number of processes on the root rank, and
    !! size 1 on other ranks (needs to be allocated even though it is
    !! not actually used.)

    PetscInt, allocatable :: array(:)
    ! Locals:
    PetscMPIInt :: rank, num_procs
    PetscInt :: size
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)

    if (rank == 0) then
       size = num_procs
    else ! have to allocate non-zero size, even if not actually used:
       size = 1
    end if
    allocate(array(size))

  end function get_mpi_int_gather_array

!------------------------------------------------------------------------

  PetscBool function array_sorted_int(a) result(sorted)
    !! Returns true if specified integer array a is monotonically
    !! increasing.

    PetscInt, intent(in) :: a(:)
    ! Locals:
    PetscInt :: i

    sorted = PETSC_TRUE
    associate(n => size(a))
      do i = 1, n - 1
         if (a(i + 1) < a(i)) then
            sorted = PETSC_FALSE
            exit
         end if
      end do
    end associate

  end function array_sorted_int

!------------------------------------------------------------------------

  PetscBool function array_sorted_real(a) result(sorted)
    !! Returns true if specified real array a is monotonically
    !! increasing.

    PetscReal, intent(in) :: a(:)
    ! Locals:
    PetscInt :: i

    sorted = PETSC_TRUE
    associate(n => size(a))
      do i = 1, n - 1
         if (a(i + 1) < a(i)) then
            sorted = PETSC_FALSE
            exit
         end if
      end do
    end associate

  end function array_sorted_real

!------------------------------------------------------------------------

  PetscReal function clock_elapsed_time(start)
    !! Returns elapsed time from start clock time, using
    !! the Fortran system_clock() function.

    integer(int32), intent(in) :: start
    ! Locals:
    integer(int32) :: end, rate

    call system_clock(end, rate)
    clock_elapsed_time = real(real(end - start, real32) / &
         real(rate, real32), dp)

  end function clock_elapsed_time

!------------------------------------------------------------------------

end module utils_module
