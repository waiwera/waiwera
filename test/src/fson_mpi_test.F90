module fson_mpi_test

  ! Tests for fson_mpi module

  use fruit
  use fson
  use kinds_module
  use mpi_module
  use fson_mpi_module

  implicit none
  private 

#include <petsc-finclude/petscsys.h>
#include <petsc-finclude/petscdef.h>

  character(len = 32), parameter :: filename = "data/fson_mpi/test_fson_mpi.json"
  public :: test_fson_mpi_int, test_fson_mpi_real, test_fson_mpi_double
  public :: test_fson_mpi_logical, test_fson_mpi_char

contains

!------------------------------------------------------------------------

  subroutine test_fson_mpi_int

    ! Test fson_get_mpi integer routines

    type(fson_value), pointer :: json
    PetscInt :: val
    PetscInt, parameter :: expected = 7
    PetscInt, allocatable :: arr(:)
    PetscInt, parameter :: expected_arr(6) = [3, -1, 4, -1, 5, -9]
    PetscInt, allocatable :: arr_2d(:,:)
    PetscInt, parameter :: expected_arr_2d(3,2) = &
         transpose(reshape([1, 2, 3, 4, 5, 6], [2,3]))

    if (mpi%rank == mpi%input_rank) then
       json => fson_parse(filename)
    end if

    call fson_get_mpi(json, "int", val=val)
    call assert_equals(expected, val, "int")
    call fson_get_mpi(json, "missing_int", expected, val)
    call assert_equals(expected, val, "int")
    call fson_get_mpi(json, "int_array", val=arr)
    call assert_equals(expected_arr, arr, size(expected_arr), "int array")
    deallocate(arr)
    call fson_get_mpi(json, "int_array_2d", val=arr_2d)
    call assert_equals(expected_arr_2d, arr_2d, size(expected_arr_2d,1), &
         size(expected_arr_2d,2), "2d int array")
    deallocate(arr_2d)

    if (mpi%rank == mpi%input_rank) then
       call fson_destroy(json)
    end if

  end subroutine test_fson_mpi_int

!------------------------------------------------------------------------

  subroutine test_fson_mpi_real

    ! Test fson_get_mpi real routines

    type(fson_value), pointer :: json
    real :: val
    real, parameter :: tol = 1.e-7
    real, parameter :: expected = 2.71818284
    real, allocatable :: arr(:)
    real, parameter :: expected_arr(5) = [200.0, -1.8, 5.22004, 78.6, -1000.5]
    real, allocatable :: arr_2d(:,:)
    real, parameter :: expected_arr_2d(2,3) = &
         transpose(reshape([1., 2., 3., 4., 5., 6.], [3,2]))

    if (mpi%rank == mpi%input_rank) then
       json => fson_parse(filename)
    end if

    call fson_get_mpi(json, "real", val=val)
    call assert_equals(expected, val, tol, "real")
    call fson_get_mpi(json, "missing_real", expected, val)
    call assert_equals(expected, val, tol, "real")
    call fson_get_mpi(json, "real_array", val=arr)
    call assert_equals(expected_arr, arr, size(expected_arr), "real array")
    deallocate(arr)
    call fson_get_mpi(json, "real_array_2d", val=arr_2d)
    call assert_equals(expected_arr_2d, arr_2d, size(expected_arr_2d,1), &
         size(expected_arr_2d,2), "2d real array")
    deallocate(arr_2d)

    if (mpi%rank == mpi%input_rank) then
       call fson_destroy(json)
    end if

  end subroutine test_fson_mpi_real

!------------------------------------------------------------------------

  subroutine test_fson_mpi_double

    ! Test fson_get_mpi double routines

    type(fson_value), pointer :: json
    real(dp) :: val
    real(dp), parameter :: tol = 1.e-10_dp
    real(dp), parameter :: expected = 1.61803398875_dp
    real(dp), allocatable :: arr(:)
    real(dp), parameter :: expected_arr(5) = [-200.0_dp, &
         1.8_dp, -5.22004_dp, -78.6_dp, 1000.5_dp]
    real(dp), allocatable :: arr_2d(:,:)
    real(dp), parameter :: expected_arr_2d(2,2) = &
         transpose(reshape([-1._dp, 1._dp, 2._dp, -3._dp], [2,2]))
    if (mpi%rank == mpi%input_rank) then
       json => fson_parse(filename)
    end if

    call fson_get_mpi(json, "double", val=val)
    call assert_equals(expected, val, tol, "double")
    call fson_get_mpi(json, "missing_double", expected, val)
    call assert_equals(expected, val, tol, "double")
    call fson_get_mpi(json, "double_array", val=arr)
    call assert_equals(expected_arr, arr, &
         size(expected_arr), "double array")
    deallocate(arr)
    call fson_get_mpi(json, "double_array_2d", val=arr_2d)
    call assert_equals(expected_arr_2d, arr_2d, size(expected_arr_2d,1), &
         size(expected_arr_2d,2), "2d double array")
    deallocate(arr_2d)

    if (mpi%rank == mpi%input_rank) then
       call fson_destroy(json)
    end if

  end subroutine test_fson_mpi_double

!------------------------------------------------------------------------

  subroutine test_fson_mpi_logical

    ! Test fson_get_mpi logical routines

    type(fson_value), pointer :: json
    PetscBool :: val
    PetscBool, parameter :: expected = .false.
    PetscBool, allocatable :: arr(:)
    PetscBool, parameter :: expected_arr(4) = [.false., .false., .true., .true.]
    PetscBool, allocatable :: arr_2d(:,:)
    PetscBool, parameter :: expected_arr_2d(3,2) = &
         transpose(reshape([.true., .false., .false., .false., .true., .true.],&
         [2,3]))

    if (mpi%rank == mpi%input_rank) then
       json => fson_parse(filename)
    end if

    call fson_get_mpi(json, "logical", val=val)
    call assert_equals(expected, val, "logical")
    call fson_get_mpi(json, "missing_logical", expected, val)
    call assert_equals(expected, val, "logical")
    call fson_get_mpi(json, "logical_array", val=arr)
    call assert_equals(expected_arr, arr, &
         size(expected_arr), "logical array")
    deallocate(arr)
    call fson_get_mpi(json, "logical_array_2d", val=arr_2d)
    call assert_equals(expected_arr_2d, arr_2d, size(expected_arr_2d,1), &
         size(expected_arr_2d,2), "2d logical array")
    deallocate(arr_2d)

    if (mpi%rank == mpi%input_rank) then
       call fson_destroy(json)
    end if

  end subroutine test_fson_mpi_logical

!------------------------------------------------------------------------

  subroutine test_fson_mpi_char

    ! Test fson_get_mpi character routines

    type(fson_value), pointer :: json
    PetscInt, parameter :: char_len = 12
    character(len = char_len) :: val
    character(len = char_len), parameter :: expected = "etaoinshrdlu"

    if (mpi%rank == mpi%input_rank) then
       json => fson_parse(filename)
    end if

    call fson_get_mpi(json, "character", val=val)
    call assert_equals(expected, val, "character")
    call fson_get_mpi(json, "missing_character", expected, val)
    call assert_equals(expected, val, "character")

    if (mpi%rank == mpi%input_rank) then
       call fson_destroy(json)
    end if

  end subroutine test_fson_mpi_char

!------------------------------------------------------------------------

end module fson_mpi_test
