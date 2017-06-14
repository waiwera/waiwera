module utils_test

  ! Tests for utils module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use fruit
  use utils_module

  implicit none
  private 

  public :: test_str_to_upper, test_str_to_lower, &
       test_split_filename, test_change_filename_extension, &
       test_int_str_len, test_str_array_index, &
       test_degrees_to_radians, test_rotation_matrix_2d, &
       test_polynomial, test_multipolynomial, test_array_pair_sum, &
       test_array_cumulative_sum, test_array_exclusive_products

contains

!------------------------------------------------------------------------

  subroutine test_str_to_upper

    ! Test str_to_upper()
    
    PetscInt, parameter :: strlen = 8
    character(len = strlen) :: str, upper_str
    character(len = strlen), parameter :: expected = "ABCABC12"
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       str = "abcABC12"
       upper_str = str_to_upper(str)
       call assert_equals(expected, upper_str, 'str_to_upper')
    end if

  end subroutine test_str_to_upper

!------------------------------------------------------------------------

  subroutine test_str_to_lower

    ! Test str_to_lower()
    
    PetscInt, parameter :: strlen = 8
    character(len = strlen) :: str, lower_str
    character(len = strlen), parameter :: expected = "abcabc12"
    PetscMPIInt :: rank
    PetscInt :: ierr
    
    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       str = "abcABC12"
       lower_str = str_to_lower(str)
       call assert_equals(expected, lower_str, 'str_to_lower')
    end if

  end subroutine test_str_to_lower

!------------------------------------------------------------------------

  subroutine test_split_filename

    ! Test split_filename()

    character(:), allocatable :: filename
    character(:), allocatable :: expected_base, expected_ext
    character(:), allocatable :: base, ext
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
    
       allocate(filename, source = "model.h5")
       allocate(expected_base, source = "model")
       allocate(expected_ext, source = "h5")
       call split_filename(filename, base, ext)
       call assert_equals(expected_base, base, 'base ' // filename)
       call assert_equals(expected_ext, ext, 'extension ' // filename)

       filename = "/path/to/model.json"
       expected_base = "/path/to/model"
       expected_ext = "json"
       call split_filename(filename, base, ext)
       call assert_equals(expected_base, base, 'base ' // filename)
       call assert_equals(expected_ext, ext, 'extension ' // filename)

       filename = "/path/to/model"
       expected_base = "/path/to/model"
       expected_ext = " "
       call split_filename(filename, base, ext)
       call assert_equals(expected_base, base, 'base ' // filename)
       call assert_equals(expected_ext, ext, 'extension ' // filename)

    end if
    
  end subroutine test_split_filename

!------------------------------------------------------------------------

  subroutine test_change_filename_extension

    ! Test change_filename_extension()

    character(:), allocatable :: filename, ext
    character(:), allocatable :: new_filename, expected_filename
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
    
       allocate(filename, source = "model.h5")
       allocate(ext, source = "log")
       allocate(expected_filename, source = "model.log")
       allocate(new_filename, source = change_filename_extension(filename, ext))
       call assert_equals(expected_filename, new_filename, 'filename')

       filename = "/path/to/model.json"
       ext = "log"
       expected_filename = "/path/to/model.log"
       new_filename = change_filename_extension(filename, ext)
       call assert_equals(expected_filename, new_filename, 'filename')

       filename = "/path/to/file"
       ext = "json"
       expected_filename = "/path/to/file.json"
       new_filename = change_filename_extension(filename, ext)
       call assert_equals(expected_filename, new_filename, 'filename')

    end if
    

  end subroutine test_change_filename_extension

!------------------------------------------------------------------------

  subroutine test_int_str_len

    ! Test int_str_len()

    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call assert_equals(1, int_str_len(0), '0')
       call assert_equals(1, int_str_len(5), '5')
       call assert_equals(2, int_str_len(79), '79')
       call assert_equals(4, int_str_len(2001), '2001')
       call assert_equals(3, int_str_len(-47), '-47')

    end if

  end subroutine test_int_str_len

!------------------------------------------------------------------------

  subroutine test_str_array_index

    ! Test str_array_index()

    PetscInt, parameter :: strlen = 3
    character(len = strlen) :: arr(3), str
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       arr = ["foo", "baz", "bar"]

       str = "foo"
       call assert_equals(1, str_array_index(str, arr), str)
       str = "bar"
       call assert_equals(3, str_array_index(str, arr), str)
       str = "pog"
       call assert_equals(-1, str_array_index(str, arr), str)

    end if

  end subroutine test_str_array_index

!------------------------------------------------------------------------

  subroutine test_degrees_to_radians

    ! Test degrees_to_radians()

    use kinds_module

    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: tol = 1.e-9_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       call assert_equals(pi / 2._dp, &
            degrees_to_radians(90._dp), tol, '90 degrees')
       call assert_equals(pi, &
            degrees_to_radians(180._dp), tol, '180 degrees')
    end if

  end subroutine test_degrees_to_radians

!------------------------------------------------------------------------

  subroutine test_rotation_matrix_2d

    ! Test rotation_matrix_2d

    use kinds_module

    PetscReal :: angle, rotation(4)
    PetscReal :: expected_rotation(4)
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: tol = 1.e-9_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       angle = 0._dp
       rotation = reshape(rotation_matrix_2d(angle), [4])
       expected_rotation = [1._dp, 0._dp, 0._dp, 1._dp]
       call assert_equals(expected_rotation, &
            rotation, 4, tol, 'angle = 0')

       angle = 1._dp
       rotation = reshape(rotation_matrix_2d(angle), [4])
       expected_rotation = [&
            0.540302305868_dp, -0.841470984808_dp,&
            0.841470984808_dp, 0.540302305868_dp]
       call assert_equals(expected_rotation, &
            rotation, 4, tol, 'angle = 1 rad')

    end if

  end subroutine test_rotation_matrix_2d

!------------------------------------------------------------------------

  subroutine test_polynomial
    ! Test polynomial

    use kinds_module

    PetscReal :: x
    PetscReal, parameter :: a(5) = [1._dp, 1._dp, 0.5_dp, &
         1._dp / 6._dp, 1._dp / 24._dp]
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: tol = 1.e-9_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       x = 0._dp
       call assert_equals(1._dp, polynomial(a, x), tol, '0')

       x = 1._dp
       call assert_equals(2.708333333333_dp, polynomial(a, x), tol, '1')

       x = -1._dp
       call assert_equals(0.375_dp, polynomial(a, x), tol, '-1')

       x = 2.3_dp
       call assert_equals(9.1388375_dp, polynomial(a, x), tol, '2.3')

    end if

  end subroutine test_polynomial

!------------------------------------------------------------------------

  subroutine test_multipolynomial
    ! Test multipolynomial

    use kinds_module

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
    PetscReal, parameter :: tol = 1.e-9_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       x = 0._dp
       call assert_equals([1._dp, 1._dp], polynomial(a, x), m, tol, '0')

       x = 1._dp
       call assert_equals([2.708333333333_dp, 0.541666666667_dp], &
            polynomial(a, x), m, tol, '1')

       x = -1._dp
       call assert_equals([0.375_dp, 0.541666666667_dp], &
            polynomial(a, x), m, tol, '-1')

       x = 2.3_dp
       call assert_equals([9.1388375_dp, -0.47899583333_dp], &
            polynomial(a, x), m, tol, '2.3')

    end if

  end subroutine test_multipolynomial

!------------------------------------------------------------------------

  subroutine test_array_pair_sum
    ! Test array_pair_sum

    use kinds_module

    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: tol = 1.e-9_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call assert_equals(2._dp, array_pair_sum([2._dp]), tol, '1-D')
       call assert_equals(12._dp, array_pair_sum([2._dp, 3._dp]), tol, '2-D')
       call assert_equals(26._dp, array_pair_sum([2._dp, 3._dp, 4._dp]), &
            tol, '3-D')

    end if

  end subroutine test_array_pair_sum

!------------------------------------------------------------------------

  subroutine test_array_cumulative_sum
    ! Test array_cumulative_sum

    use kinds_module

    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: tol = 1.e-9_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call assert_equals([2._dp], &
            array_cumulative_sum([2._dp]), 1, tol, '1-D')
       call assert_equals([2._dp, 5._dp], &
            array_cumulative_sum([2._dp, 3._dp]), 2, tol, '2-D')
       call assert_equals([2._dp, 5._dp, 9._dp], &
            array_cumulative_sum([2._dp, 3._dp, 4._dp]), 3, &
            tol, '3-D')

    end if

  end subroutine test_array_cumulative_sum

!------------------------------------------------------------------------

  subroutine test_array_exclusive_products
    ! Test array_exclusive_products

    use kinds_module

    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: tol = 1.e-9_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call assert_equals([1._dp], &
            array_exclusive_products([2._dp]), 1, tol, '1-D')
       call assert_equals([3._dp, 2._dp], &
            array_exclusive_products([2._dp, 3._dp]), 2, tol, '2-D')
       call assert_equals([12._dp, 8._dp, 6._dp], &
            array_exclusive_products([2._dp, 3._dp, 4._dp]), 3, &
            tol, '3-D')

    end if

  end subroutine test_array_exclusive_products

!------------------------------------------------------------------------

end module utils_test
