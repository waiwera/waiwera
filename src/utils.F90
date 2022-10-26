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

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module

  implicit none
  private

  PetscReal, parameter, public :: pi = 4._dp * atan(1._dp)

  interface polynomial
     module procedure polynomial_single
     module procedure polynomial_multiple
  end interface polynomial

  interface polynomial_derivative
     module procedure polynomial_derivative_single
     module procedure polynomial_derivative_multiple
  end interface polynomial_derivative

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
       polynomial, polynomial_derivative, &
       array_pair_sum, array_cumulative_sum, &
       array_exclusive_products, array_sorted, &
       array_indices_in_int_array, clock_elapsed_time, &
       array_is_permutation, invert_indices, &
       array_progressive_limit, newton1d

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

  function polynomial_derivative_single(a) result(da)
    !! Takes coefficient array of a polynomial and returns the
    !! coefficient array of its derivative.

    PetscReal, intent(in) :: a(:)
    PetscReal :: da(size(a) - 1)
    ! Locals:
    PetscInt :: i

    do i = 1, size(a) - 1
       da(i) = i * a(i + 1)
    end do

  end function polynomial_derivative_single

!------------------------------------------------------------------------

  function polynomial_derivative_multiple(a) result(da)
    !! Takes coefficient array for multiple polynomials and returns
    !! the coefficient array of their derivatives.

    PetscReal, intent(in) :: a(:, :)
    PetscReal :: da(size(a, 1), size(a, 2) - 1)
    ! Locals:
    PetscInt :: i

    do i = 1, size(a, 2) - 1
       da(:, i) = i * a(:, i + 1)
    end do

  end function polynomial_derivative_multiple

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

  function array_indices_in_int_array(a, b) result(indices)
    !! Returns (1-based) indices of the elements of integer array b in
    !! array a. It is assumed that a and b are permutations of each
    !! other.

    PetscInt, intent(in) :: a(:), b(:)
    PetscInt :: indices(size(a))
    ! Locals:
    PetscInt :: fa(size(a)), fb(size(a)), fbi(size(a))
    PetscInt :: i
    PetscErrorCode :: ierr

    associate (n => size(a))

      ! Find sort permutations for a and b:
      fa = [(i, i = 0, n - 1)]
      call PetscSortIntWithPermutation(n, a, fa, ierr); CHKERRQ(ierr)
      fa = fa + 1
      fb = [(i, i = 0, n - 1)]
      call PetscSortIntWithPermutation(n, b, fb, ierr); CHKERRQ(ierr)
      fb = fb + 1

      ! Invert permutation for b:
      fbi = -1
      do i = 1, n
         fbi(fb(i)) = i
      end do

      do i = 1, n
         indices(i) = fa(fbi(i))
      end do

    end associate

  end function array_indices_in_int_array

!------------------------------------------------------------------------

  PetscReal function clock_elapsed_time(start)
    !! Returns elapsed time from start clock time, using
    !! the Fortran system_clock() function.

    use iso_fortran_env, only: int32, real32

    integer(int32), intent(in) :: start
    ! Locals:
    integer(int32) :: end, rate

    call system_clock(end, rate)
    clock_elapsed_time = real(real(end - start, real32) / &
         real(rate, real32), dp)

  end function clock_elapsed_time

!------------------------------------------------------------------------

  PetscBool function array_is_permutation(a)
    !! Returns true if the integer array a is a permutation.

    PetscInt, intent(in) :: a(:)
    ! Locals:
    PetscInt, allocatable :: count(:)
    PetscInt :: i, amin, amax

    associate(n => size(a))
      amin = minval(a)
      amax = amin + n - 1
      allocate(count(amin: amax))
      count = 0
      array_is_permutation = PETSC_TRUE
      do i = 1, n
         if (a(i) <= amax) then
            count(a(i)) = count(a(i)) + 1
         else
            array_is_permutation = PETSC_FALSE
            exit
         end if
      end do
      array_is_permutation = (array_is_permutation .and. all(count == 1))
      deallocate(count)
    end associate

  end function array_is_permutation

!------------------------------------------------------------------------

  IS function invert_indices(indices, name)
    !! Returns index set containing the inverse of the mapping
    !! represented by the values contained in the indices array on all
    !! processes. The resulting index set contains values only on the
    !! root process and is given the specified object name.

    use mpi_utils_module, only: get_mpi_int_gather_array

    PetscInt, intent(in) :: indices(:)
    character(*), intent(in) :: name
    ! Locals:
    PetscInt :: i, local_size
    PetscMPIInt :: rank, num_procs, num_all, is_count
    PetscInt, allocatable :: counts(:), displacements(:)
    PetscInt, allocatable :: indices_all(:), global_indices(:)
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)
    local_size = size(indices)

    counts = get_mpi_int_gather_array()
    displacements = get_mpi_int_gather_array()
    call MPI_gather(local_size, 1, MPI_INTEGER, counts, 1, &
         MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       displacements = [[0], &
            array_cumulative_sum(counts(1: num_procs - 1))]
       num_all = sum(counts)
       is_count = num_all
    else
       num_all = 1
       is_count = 0
    end if

    allocate(indices_all(0: num_all - 1), global_indices(0: num_all - 1))
    global_indices = -1
    call MPI_gatherv(indices, local_size, MPI_INTEGER, &
         indices_all, counts, displacements, &
         MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       do i = 0, num_all - 1
          global_indices(indices_all(i)) = i
       end do
    end if
    deallocate(indices_all, counts, displacements)

    call ISCreateGeneral(PETSC_COMM_WORLD, is_count, &
         global_indices, PETSC_COPY_VALUES, invert_indices, ierr)
    CHKERRQ(ierr)
    call PetscObjectSetName(invert_indices, name, ierr)
    deallocate(global_indices)

  end function invert_indices

!------------------------------------------------------------------------

  function array_progressive_limit(a, total, order) result(limit)
    !! Returns array of limits required to apply to real array a
    !! progressively so that it sums to the specified total: each
    !! element of a is limited in order until the required total is
    !! reached. The optional order array specifies a permutation so
    !! that the limiting can be carried out in any order.

    PetscReal, intent(in) :: a(:) !! Quantities to limit
    PetscReal, intent(in) :: total !! Target total of a array
    PetscInt, intent(in), optional :: order(:) !! Permutation for limiting order
    PetscReal :: limit(size(a)) !! Output limits for a
    ! Locals:
    PetscInt :: limit_order(size(a))
    PetscInt :: i, j
    PetscReal :: sum_a, next_sum_a

    associate(n => size(a))

      if (present(order)) then
         limit_order = order
      else
         limit_order = [(i, i = 1, n)]
      end if

      limit = 0._dp
      sum_a = 0._dp
      do i = 1, n
         j = limit_order(i)
         next_sum_a = sum_a + a(j)
         if (next_sum_a > total) then
            limit(j) = total - sum_a
            exit
         else
            limit(j) = a(j)
            sum_a = next_sum_a
         end if
      end do

    end associate

  end function array_progressive_limit

!------------------------------------------------------------------------

  subroutine newton1d(f, x, x_increment, tolerance, max_iterations, err)
    !! 1-D Newton solve to find f(x) = 0, for the specified relative
    !! variable increment, function tolerance and maximum number of
    !! iterations. The error flag returns nonzero if there were any
    !! errors in function evaluation or the iteration limit was
    !! exceeded.

    interface
       PetscReal function f(x, err)
         PetscReal, intent(in) :: x
         PetscErrorCode, intent(out) :: err
       end function f
    end interface

    PetscReal, intent(in out) :: x !! Starting and final variable value
    PetscReal, intent(in) :: x_increment !! Relative variable increment
    PetscReal, intent(in) :: tolerance !! Function tolerance
    PetscInt, intent(in) :: max_iterations !! Maximum number of iterations
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: dx, fx, fxd, df
    PetscInt :: i
    PetscBool :: found

    dx = x_increment * x
    found = PETSC_FALSE

    do i = 1, max_iterations
       fx = f(x, err)
       if (err == 0) then
          if (abs(fx) <= tolerance) then
             found = PETSC_TRUE
             exit
          else
             fxd = f(x + dx, err)
             if (err == 0) then
                df = (fxd - fx) / dx
                x = x - fx / df
             else
                exit
             end if
          end if
       else
          exit
       end if
    end do

    if ((err == 0) .and. (.not.(found))) then
       err = 1
    end if

  end subroutine newton1d

!------------------------------------------------------------------------

end module utils_module
