module utils_test

  ! Tests for utils module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use zofu
  use utils_module
  use kinds_module

  implicit none
  private 

  public :: setup, teardown, setup_test
  public :: test_str_to_lower, &
       test_split_filename, test_change_filename_extension, &
       test_int_str_len, test_str_array_index, &
       test_degrees_to_radians, test_rotation_matrix_2d, &
       test_polynomial, test_multipolynomial, test_polynomial_derivative, &
       test_polynomial_integral, &
       test_array_pair_sum, test_array_cumulative_sum, &
       test_array_exclusive_products, test_array_sorted, &
       test_array_indices_in_int_array, test_is_permutation, &
       test_array_progressive_limit, test_newton1d, &
       test_is_permutation_of, test_array_unique, test_sign_test, &
       test_hermite_spline

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

  subroutine setup_test(test)

    class(unit_test_type), intent(in out) :: test

    test%tolerance = 1.e-9

  end subroutine setup_test

!------------------------------------------------------------------------

  subroutine test_str_to_lower(test)

    ! Test str_to_lower()
    
    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: strlen = 8
    character(len = strlen) :: str, lower_str
    character(len = strlen), parameter :: expected = "abcabc12"
    PetscMPIInt :: rank
    PetscInt :: ierr
    
    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       str = "abcABC12"
       lower_str = str_to_lower(str)
       call test%assert(expected, lower_str, 'str_to_lower')
    end if

  end subroutine test_str_to_lower

!------------------------------------------------------------------------

  subroutine test_split_filename(test)

    ! Test split_filename()

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    character(:), allocatable :: filename
    character(:), allocatable :: expected_base, expected_ext
    character(:), allocatable :: base, ext
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
    
       filename  = "model.h5"
       expected_base = "model"
       expected_ext = "h5"
       call split_filename(filename, base, ext)
       call test%assert(expected_base, base, 'base ' // filename)
       call test%assert(expected_ext, ext, 'extension ' // filename)

       filename = "/path/to/model.json"
       expected_base = "/path/to/model"
       expected_ext = "json"
       call split_filename(filename, base, ext)
       call test%assert(expected_base, base, 'base ' // filename)
       call test%assert(expected_ext, ext, 'extension ' // filename)

       filename = "/path/to/model"
       expected_base = "/path/to/model"
       expected_ext = " "
       call split_filename(filename, base, ext)
       call test%assert(expected_base, base, 'base ' // filename)
       call test%assert(expected_ext, ext, 'extension ' // filename)

    end if
    
  end subroutine test_split_filename

!------------------------------------------------------------------------

  subroutine test_change_filename_extension(test)

    ! Test change_filename_extension()

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    character(:), allocatable :: filename, ext
    character(:), allocatable :: new_filename, expected_filename
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
    
       filename = "model.h5"
       ext = "log"
       expected_filename = "model.log"
       new_filename = change_filename_extension(filename, ext)
       call test%assert(expected_filename, new_filename, 'filename')

       filename = "/path/to/model.json"
       ext = "log"
       expected_filename = "/path/to/model.log"
       new_filename = change_filename_extension(filename, ext)
       call test%assert(expected_filename, new_filename, 'filename')

       filename = "/path/to/file"
       ext = "json"
       expected_filename = "/path/to/file.json"
       new_filename = change_filename_extension(filename, ext)
       call test%assert(expected_filename, new_filename, 'filename')

    end if
    

  end subroutine test_change_filename_extension

!------------------------------------------------------------------------

  subroutine test_int_str_len(test)

    ! Test int_str_len()

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert(1, int_str_len(0), '0')
       call test%assert(1, int_str_len(5), '5')
       call test%assert(2, int_str_len(79), '79')
       call test%assert(4, int_str_len(2001), '2001')
       call test%assert(3, int_str_len(-47), '-47')

    end if

  end subroutine test_int_str_len

!------------------------------------------------------------------------

  subroutine test_str_array_index(test)

    ! Test str_array_index()

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: strlen = 3
    character(len = strlen) :: arr(3), str
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       arr = ["foo", "baz", "bar"]

       str = "foo"
       call test%assert(1, str_array_index(str, arr), str)
       str = "bar"
       call test%assert(3, str_array_index(str, arr), str)
       str = "pog"
       call test%assert(-1, str_array_index(str, arr), str)

    end if

  end subroutine test_str_array_index

!------------------------------------------------------------------------

  subroutine test_degrees_to_radians(test)

    ! Test degrees_to_radians()

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       call test%assert(pi / 2._dp, &
            degrees_to_radians(90._dp), '90 degrees')
       call test%assert(pi, &
            degrees_to_radians(180._dp), '180 degrees')
    end if

  end subroutine test_degrees_to_radians

!------------------------------------------------------------------------

  subroutine test_rotation_matrix_2d(test)

    ! Test rotation_matrix_2d

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscReal :: angle, rotation(4)
    PetscReal :: expected_rotation(4)
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       angle = 0._dp
       rotation = reshape(rotation_matrix_2d(angle), [4])
       expected_rotation = [1._dp, 0._dp, 0._dp, 1._dp]
       call test%assert(expected_rotation, &
            rotation, 'angle = 0')

       angle = 1._dp
       rotation = reshape(rotation_matrix_2d(angle), [4])
       expected_rotation = [&
            0.540302305868_dp, -0.841470984808_dp,&
            0.841470984808_dp, 0.540302305868_dp]
       call test%assert(expected_rotation, &
            rotation, 'angle = 1 rad')

    end if

  end subroutine test_rotation_matrix_2d

!------------------------------------------------------------------------

  subroutine test_polynomial(test)
    ! Test polynomial

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscReal :: x
    PetscReal, parameter :: a(5) = [1._dp, 1._dp, 0.5_dp, &
         1._dp / 6._dp, 1._dp / 24._dp]
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       x = 0._dp
       call test%assert(1._dp, polynomial(a, x), '0')

       x = 1._dp
       call test%assert(2.708333333333_dp, polynomial(a, x), '1')

       x = -1._dp
       call test%assert(0.375_dp, polynomial(a, x), '-1')

       x = 2.3_dp
       call test%assert(9.1388375_dp, polynomial(a, x), '2.3')

    end if

  end subroutine test_polynomial

!------------------------------------------------------------------------

  subroutine test_multipolynomial(test)
    ! Test multipolynomial

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscReal :: x
    PetscInt, parameter :: m = 2, n = 5
    PetscReal, parameter :: a(m, n) = reshape([&
         1._dp, 1._dp, &
         1._dp, 0._dp, &
         0.5_dp, -0.5_dp, &
         1._dp / 6._dp, 0._dp, &
         1._dp / 24._dp, 1._dp / 24._dp], &
         [m, n])
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       x = 0._dp
       call test%assert([1._dp, 1._dp], polynomial(a, x), '0')

       x = 1._dp
       call test%assert([2.708333333333_dp, 0.541666666667_dp], &
            polynomial(a, x), '1')

       x = -1._dp
       call test%assert([0.375_dp, 0.541666666667_dp], &
            polynomial(a, x), '-1')

       x = 2.3_dp
       call test%assert([9.1388375_dp, -0.47899583333_dp], &
            polynomial(a, x), '2.3')

    end if

  end subroutine test_multipolynomial

!------------------------------------------------------------------------

  subroutine test_polynomial_derivative(test)
    ! Test polynomial derivative

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscReal, parameter :: a(5) = [1._dp, 1._dp, 0.5_dp, &
         1._dp / 6._dp, 1._dp / 24._dp]
    PetscReal, allocatable :: da(:)
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       da = polynomial_derivative(a)
       call test%assert([1._dp, 1._dp, 0.5_dp, 1._dp / 6._dp], da, 'coefs')

    end if

  end subroutine test_polynomial_derivative

!------------------------------------------------------------------------

  subroutine test_polynomial_integral(test)
    ! Test polynomial integral

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscReal, parameter :: a(4) = [-1._dp, 2._dp, 0.5_dp, -0.25_dp]
    PetscReal, allocatable :: ai(:), a2(:)
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       ai = polynomial_integral(a)
       call test%assert([0._dp, -1._dp, 1._dp, 1._dp / 6._dp, &
            -1._dp / 16._dp], ai, 'coefs')
       a2 = polynomial_derivative(ai)
       call test%assert(a, a2, 'deriv of integral')

    end if

  end subroutine test_polynomial_integral

!------------------------------------------------------------------------

  subroutine test_array_pair_sum(test)
    ! Test array_pair_sum

    use kinds_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert(2._dp, array_pair_sum([2._dp]), '1-D')
       call test%assert(12._dp, array_pair_sum([2._dp, 3._dp]), '2-D')
       call test%assert(26._dp, array_pair_sum([2._dp, 3._dp, 4._dp]), '3-D')

    end if

  end subroutine test_array_pair_sum

!------------------------------------------------------------------------

  subroutine test_array_cumulative_sum(test)
    ! Test array_cumulative_sum

    use kinds_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert([PetscInt::], &
            array_cumulative_sum([PetscInt::]), '0-D')
       call test%assert([2._dp], &
            array_cumulative_sum([2._dp]), '1-D')
       call test%assert([2._dp, 5._dp], &
            array_cumulative_sum([2._dp, 3._dp]), '2-D')
       call test%assert([2._dp, 5._dp, 9._dp], &
            array_cumulative_sum([2._dp, 3._dp, 4._dp]), '3-D')
       call test%assert([2, 5, 9], &
            array_cumulative_sum([2, 3, 4]), 'integer')
    end if

  end subroutine test_array_cumulative_sum

!------------------------------------------------------------------------

  subroutine test_array_exclusive_products(test)
    ! Test array_exclusive_products

    use kinds_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert([1._dp], &
            array_exclusive_products([2._dp]), '1-D')
       call test%assert([3._dp, 2._dp], &
            array_exclusive_products([2._dp, 3._dp]), '2-D')
       call test%assert([12._dp, 8._dp, 6._dp], &
            array_exclusive_products([2._dp, 3._dp, 4._dp]), '3-D')

    end if

  end subroutine test_array_exclusive_products

!------------------------------------------------------------------------

    subroutine test_array_sorted(test)
    ! Test array sorted

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert(array_sorted([-1, 3, 4]), '[-1, 3, 4]')
       call test%assert(array_sorted([-1._dp, 3._dp, 4._dp]), '[-1., 3., 4.]')
       call test%assert(.not. array_sorted([1, 3, 2, 5]), '[1, 3, 2, 5]')
       call test%assert(.not. array_sorted([1._dp, 3._dp, 2._dp, 5._dp]), '[1., 3., 2., 5.]')
       call test%assert(array_sorted([1._dp, 3._dp, 3._dp]), '[1., 3., 3.]')

    end if

  end subroutine test_array_sorted

!------------------------------------------------------------------------

  subroutine test_array_indices_in_int_array(test)
    ! Test array_indices_in_int_array

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert([1], array_indices_in_int_array([3], [3]), 'case 1')
       call test%assert([1, 2], array_indices_in_int_array([0, 1], [0, 1]), 'case 2')
       call test%assert([2, 1], array_indices_in_int_array([1, 0], [0, 1]), 'case 3')
       call test%assert([2, 3, 1], array_indices_in_int_array([3, 4, 1], [4, 1, 3]), 'case 4')

    end if

  end subroutine test_array_indices_in_int_array

!------------------------------------------------------------------------

  subroutine test_is_permutation(test)
    ! Test array_is_permutation

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscInt, allocatable :: a(:)

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert(array_is_permutation([0]), 'case 1')
       call test%assert(array_is_permutation([1]), 'case 2')
       call test%assert(array_is_permutation([1, 2]), 'case 3')
       call test%assert(array_is_permutation([0, 1]), 'case 4')
       call test%assert(.not. array_is_permutation([1, 1]), 'case 5')
       call test%assert(array_is_permutation([2, 1]), 'case 6')

       allocate(a(0: 2))
       a(0) = 1; a(1) = 2; a(2) = 0
       call test%assert(array_is_permutation(a), 'case 7')
       deallocate(a)

       allocate(a(0: 2))
       a(0) = 2; a(1) = 3; a(2) = 1
       call test%assert(array_is_permutation(a), 'case 8')
       deallocate(a)

       call test%assert(.not. array_is_permutation([1, 2, 4]), 'case 9')
       call test%assert(.not. array_is_permutation([2, 1, 2]), 'case 10')

       call test%assert(array_is_permutation([3, 0, 5, 4, 1, 2]), 'case 11')

    end if

  end subroutine test_is_permutation

!------------------------------------------------------------------------

  subroutine test_is_permutation_of(test)
    ! Test array_is_permutation_of

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert(array_is_permutation_of([0], [0]), 'case 1')
       call test%assert(.not. array_is_permutation_of([0], [1]), 'case 2')
       call test%assert(.not. array_is_permutation_of([0], [0, 1]), 'case 3')
       call test%assert(array_is_permutation_of([1, 3, -2], [3, -2, 1]), 'case 4')
       call test%assert(array_is_permutation_of([1, 3, 3], [3, 1, 3]), 'case 5')

    end if

  end subroutine test_is_permutation_of

!------------------------------------------------------------------------

  subroutine test_array_unique(test)
    ! Test array_uniquie

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscInt, allocatable :: empty(:)

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       allocate(empty(0))
       call test%assert(empty, array_unique(empty), 'case 0')
       call test%assert([0], array_unique([0]), 'case 1')
       call test%assert([0, 1], array_unique([0, 1]), 'case 2')
       call test%assert([0, 1], array_unique([0, 1, 0]), 'case 3')
       call test%assert([0, 1], array_unique([0, 1, 0, 1]), 'case 4')
       call test%assert([-1], array_unique([-1, -1, -1]), 'case 5')
    end if

  end subroutine test_array_unique

!------------------------------------------------------------------------


  subroutine test_array_progressive_limit(test)
    ! Test array progressive limit

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert([1._dp, 2._dp, 3._dp], &
            array_progressive_limit([1._dp, 2._dp, 3._dp], 10._dp), 'case 1')
       call test%assert([1._dp, 2._dp, 1.5_dp], &
            array_progressive_limit([1._dp, 2._dp, 3._dp], 4.5_dp), 'case 2')
       call test%assert([1._dp, 2._dp, 0._dp], &
            array_progressive_limit([1._dp, 2._dp, 3._dp], 3._dp), 'case 3')
       call test%assert([1._dp, 1._dp, 0._dp], &
            array_progressive_limit([1._dp, 2._dp, 3._dp], 2._dp), 'case 4')
       call test%assert([0.25_dp, 0._dp, 0._dp], &
            array_progressive_limit([1._dp, 2._dp, 3._dp], 0.25_dp), 'case 5')
       call test%assert([0._dp, 0._dp, 0._dp], &
            array_progressive_limit([1._dp, 2._dp, 3._dp], 0._dp), 'case 6')
       call test%assert([1._dp, 2._dp, 3._dp], &
            array_progressive_limit([1._dp, 2._dp, 3._dp], 10._dp, &
            [3, 2, 1]), 'case 7')
       call test%assert([0.5_dp, 2._dp, 3._dp], &
            array_progressive_limit([1._dp, 2._dp, 3._dp], 5.5_dp, &
            [3, 2, 1]), 'case 8')
       call test%assert([0._dp, 2._dp, 3._dp], &
            array_progressive_limit([1._dp, 2._dp, 3._dp], 5._dp, &
            [2, 3, 1]), 'case 9')
       call test%assert([0._dp, 0.2_dp, 0._dp], &
            array_progressive_limit([1._dp, 2._dp, 3._dp], 0.2_dp, &
            [2, 3, 1]), 'case 10')
       call test%assert([0._dp, 0._dp, 0._dp], &
            array_progressive_limit([1._dp, 2._dp, 3._dp], 0._dp, &
            [3, 1, 2]), 'case 11')

    end if

  end subroutine test_array_progressive_limit

!------------------------------------------------------------------------

  subroutine test_newton1d(test)
    ! Test newton1d

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal :: x
    PetscReal, parameter :: poly(3) = [-15._dp, 2._dp, 1._dp]
    PetscReal, parameter :: inc = 1.e-8_dp
    PetscReal, parameter :: ftol = 1.e-10_dp, xtol = 1.e-10_dp
    PetscErrorCode :: err

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       x = 0._dp
       call newton1d(f1, x, ftol, xtol, 1, inc, err)
       call test%assert(0._dp, x, 'case 1 value')
       call test%assert(0, err, 'case 1 error')

       x = 1._dp
       call newton1d(f1, x, ftol, xtol, 3, inc, err)
       call test%assert(0._dp, x, 'case 2 value')
       call test%assert(0, err, 'case 2 error')

       x = 1._dp
       call newton1d(f2, x, ftol, xtol, 5, inc, err)
       call test%assert(sqrt(2._dp), x, 'case 3 value')
       call test%assert(0, err, 'case 3 error')

       x = 1._dp
       call newton1d(f2, df2, x, ftol, xtol, 5, err)
       call test%assert(sqrt(2._dp), x, 'case 4 value')
       call test%assert(0, err, 'case 4 error')

       x = 1._dp
       call newton1d(poly, x, ftol, xtol, 6, err)
       call test%assert(3._dp, x, 'case 5 value')
       call test%assert(0, err, 'case 5 error')

    end if

  contains

    PetscReal function f1(x, err)
      PetscReal, intent(in) :: x
      PetscErrorCode, intent(out) :: err
      f1 = x
      err = 0
    end function f1

    PetscReal function f2(x, err)
      PetscReal, intent(in) :: x
      PetscErrorCode, intent(out) :: err
      f2 = x * x - 2._dp
      err = 0
    end function f2

    PetscReal function df2(x, err)
      PetscReal, intent(in) :: x
      PetscErrorCode, intent(out) :: err
      df2 = 2._dp * x
      err = 0
    end function df2

  end subroutine test_newton1d

!------------------------------------------------------------------------

  subroutine test_sign_test(test)
    ! Test sign_test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert(1, sign_test(1._dp, 1._dp), '1,1')
       call test%assert(0, sign_test(0._dp, 0._dp), '0,0')
       call test%assert(0, sign_test(0._dp, 2._dp), '0,2')
       call test%assert(1, sign_test(-2._dp, -1._dp), '-2,-1')
       call test%assert(-1, sign_test(-1.5_dp, 2._dp), '-1.5,2')
       call test%assert(-1, sign_test(1.5_dp, -2._dp), '1.5,-2')

    end if

  end subroutine test_sign_test

!------------------------------------------------------------------------

  subroutine test_hermite_spline(test)
    ! Test hermite_spline

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert(1._dp, hermite_spline(0._dp), '0')
       call test%assert(0._dp, hermite_spline(1._dp), '1')
       call test%assert(0.5_dp, hermite_spline(0.5_dp), '0.5')
       call test%assert(0.896_dp, hermite_spline(0.2_dp), '0.2')
       call test%assert(0.216_dp, hermite_spline(0.7_dp), '0.7')

    end if

  end subroutine test_hermite_spline

!------------------------------------------------------------------------

end module utils_test
