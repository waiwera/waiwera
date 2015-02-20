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

  public :: test_fson_mpi

contains

!------------------------------------------------------------------------

  subroutine test_fson_mpi

    ! Test fson_get_mpi

    type(fson_value), pointer :: json
    character(len = 24), parameter :: filename = "data/test_fson_mpi.json"
    integer, parameter :: char_len = 12
    real, parameter :: real_tol = 1.e-7
    real(dp), parameter :: double_tol = 1.e-10_dp
    integer :: int_val
    real :: real_val
    real(dp) :: double_val
    logical :: logical_val
    character(len = char_len) :: char_val
    integer, parameter :: expected_int = 7
    real, parameter :: expected_real = 2.71818284
    real(dp), parameter :: expected_double = 1.61803398875_dp
    logical, parameter :: expected_logical = .false.
    character(len = char_len), parameter :: expected_char = "etaoinshrdlu"
    integer, parameter :: missing_int = -3
    integer, allocatable :: int_array(:)
    integer, parameter :: expected_int_array(6) = [3, -1, 4, -1, 5, -9]
    real, allocatable :: real_array(:)
    real, parameter :: expected_real_array(5) = [200.0, -1.8, 5.22004, 78.6, -1000.5]
    real(dp), allocatable :: double_array(:)
    real(dp), parameter :: expected_double_array(5) = [-200.0_dp, &
         1.8_dp, -5.22004_dp, -78.6_dp, 1000.5_dp]
    logical, allocatable :: logical_array(:)
    logical, parameter :: expected_logical_array(4) = [.false., .false., .true., .true.]

    if (mpi%rank == mpi%input_rank) then
       json => fson_parse(filename)
    end if

    call fson_get_mpi(json, "int", 0, int_val)
    call assert_equals(expected_int, int_val, "int")
    call fson_get_mpi(json, "missing_int", expected_int, int_val)
    call assert_equals(expected_int, int_val, "int")
    call fson_get_mpi(json, "int_array", [0], int_array)
    call assert_equals(expected_int_array, int_array, &
         size(expected_int_array), "int array")
    deallocate(int_array)

    call fson_get_mpi(json, "real", 0.0, real_val)
    call assert_equals(expected_real, real_val, real_tol, "real")
    call fson_get_mpi(json, "missing_real", expected_real, real_val)
    call assert_equals(expected_real, real_val, real_tol, "real")
    call fson_get_mpi(json, "real_array", [0.0], real_array)
    call assert_equals(expected_real_array, real_array, &
         size(expected_real_array), "real array")
    deallocate(real_array)

    call fson_get_mpi(json, "double", 0.0_dp, double_val)
    call assert_equals(expected_double, double_val, double_tol, "double")
    call fson_get_mpi(json, "missing_double", expected_double, double_val)
    call assert_equals(expected_double, double_val, double_tol, "double")
    call fson_get_mpi(json, "double_array", [0.0_dp], double_array)
    call assert_equals(expected_double_array, double_array, &
         size(expected_double_array), "double array")
    deallocate(double_array)

    call fson_get_mpi(json, "logical", .true., logical_val)
    call assert_equals(expected_logical, logical_val, "logical")
    call fson_get_mpi(json, "missing_logical", expected_logical, logical_val)
    call assert_equals(expected_logical, logical_val, "logical")
    call fson_get_mpi(json, "logical_array", [.true.], logical_array)
    call assert_equals(expected_logical_array, logical_array, &
         size(expected_logical_array), "logical array")
    deallocate(logical_array)

    call fson_get_mpi(json, "character", "", char_val)
    call assert_equals(expected_char, char_val, "character")
    call fson_get_mpi(json, "missing_character", expected_char, char_val)
    call assert_equals(expected_char, char_val, "character")

    if (mpi%rank == mpi%input_rank) then
       call fson_destroy(json)
    end if

  end subroutine test_fson_mpi

!------------------------------------------------------------------------

end module fson_mpi_test
