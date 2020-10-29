module fson_mpi_test

  ! Tests for fson_mpi module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use zofu
  use fson
  use kinds_module
  use fson_mpi_module

  implicit none
  private 

  character(len = 512) :: data_path, filename
  public :: setup, teardown, setup_test
  public :: test_fson_mpi_int, test_fson_mpi_real, test_fson_mpi_double
  public :: test_fson_mpi_logical, test_fson_mpi_char
  public :: test_fson_mpi_array_rank, test_fson_get_name_mpi
  public :: test_fson_value_next_mpi

contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: ios

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()

    call get_environment_variable('WAIWERA_TEST_DATA_PATH', &
         data_path, status = ios)
    if (ios /= 0) data_path = ''
    filename = trim(adjustl(data_path)) // "fson_mpi/test_fson_mpi.json"

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown()

    PetscErrorCode :: ierr

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

  subroutine setup_test(test)

    class(unit_test_type), intent(in out) :: test

    test%tolerance = 1.e-7

  end subroutine setup_test

!------------------------------------------------------------------------

  subroutine test_fson_mpi_int(test)

    ! Test fson_get_mpi integer routines

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    PetscInt :: val
    PetscInt, parameter :: expected = 7
    PetscInt, allocatable :: arr(:)
    PetscInt, parameter :: expected_arr(6) = [3, -1, 4, -1, 5, -9]
    PetscInt, allocatable :: arr_2d(:,:)
    PetscInt, parameter :: expected_arr_2d(3,2) = &
         transpose(reshape([1, 2, 3, 4, 5, 6], [2,3]))

    json => fson_parse_mpi(filename)

    call fson_get_mpi(json, "int", val=val)
    call test%assert(expected, val, "int")
    call fson_get_mpi(json, "missing_int", expected, val)
    call test%assert(expected, val, "int")
    call fson_get_mpi(json, "int_array", val=arr)
    call test%assert(expected_arr, arr, "int array")
    deallocate(arr)
    call fson_get_mpi(json, "int_array_2d", val=arr_2d)
    call test%assert(expected_arr_2d, arr_2d, "2d int array")
    deallocate(arr_2d)
    call fson_get_mpi(json, "empty_array", val=arr)
    call test%assert(PETSC_FALSE, allocated(arr), "empty array")

    call fson_destroy_mpi(json)

  end subroutine test_fson_mpi_int

!------------------------------------------------------------------------

  subroutine test_fson_mpi_real(test)

    ! Test fson_get_mpi real routines

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    real :: val
    real, parameter :: expected = 2.71818284
    real, allocatable :: arr(:)
    real, parameter :: expected_arr(5) = [200.0, -1.8, 5.22004, 78.6, -1000.5]
    real, allocatable :: arr_2d(:,:)
    real, parameter :: expected_arr_2d(2,3) = &
         transpose(reshape([1., 2., 3., 4., 5., 6.], [3,2]))

    json => fson_parse_mpi(filename)

    call fson_get_mpi(json, "real", val=val)
    call test%assert(expected, val, "real")
    call fson_get_mpi(json, "missing_real", expected, val)
    call test%assert(expected, val, "real")
    call fson_get_mpi(json, "real_array", val=arr)
    call test%assert(expected_arr, arr, "real array")
    deallocate(arr)
    call fson_get_mpi(json, "real_array_2d", val=arr_2d)
    call test%assert(expected_arr_2d, arr_2d, "2d real array")
    deallocate(arr_2d)
    call fson_get_mpi(json, "empty_array", val=arr)
    call test%assert(PETSC_FALSE, allocated(arr), "empty array")

    call fson_destroy_mpi(json)

  end subroutine test_fson_mpi_real

!------------------------------------------------------------------------

  subroutine test_fson_mpi_double(test)

    ! Test fson_get_mpi double routines

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    PetscReal :: val
    PetscReal, parameter :: tol = 1.e-10_dp
    PetscReal, parameter :: expected = 1.61803398875_dp
    PetscReal, allocatable :: arr(:)
    PetscReal, parameter :: expected_arr(5) = [-200.0_dp, &
         1.8_dp, -5.22004_dp, -78.6_dp, 1000.5_dp]
    PetscReal, allocatable :: arr_2d(:,:)
    PetscReal, parameter :: expected_arr_2d(2,2) = &
         transpose(reshape([-1._dp, 1._dp, 2._dp, -3._dp], [2,2]))

    json => fson_parse_mpi(filename)

    call fson_get_mpi(json, "double", val=val)
    call test%assert(expected, val, "double", tol = tol)
    call fson_get_mpi(json, "missing_double", expected, val)
    call test%assert(expected, val, "double", tol = tol)
    call fson_get_mpi(json, "double_array", val=arr)
    call test%assert(expected_arr, arr, "double array")
    deallocate(arr)
    call fson_get_mpi(json, "double_array_2d", val=arr_2d)
    call test%assert(expected_arr_2d, arr_2d, "2d double array")
    deallocate(arr_2d)
    call fson_get_mpi(json, "empty_array", val=arr)
    call test%assert(PETSC_FALSE, allocated(arr), "empty array")

    call fson_destroy_mpi(json)

  end subroutine test_fson_mpi_double

!------------------------------------------------------------------------

  subroutine test_fson_mpi_logical(test)

    ! Test fson_get_mpi logical routines

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    PetscBool :: val
    PetscBool, parameter :: expected = PETSC_FALSE
    PetscBool, allocatable :: arr(:)
    PetscBool, parameter :: expected_arr(4) = [PETSC_FALSE, PETSC_FALSE, &
         PETSC_TRUE, PETSC_TRUE]
    PetscBool, allocatable :: arr_2d(:,:)
    PetscBool, parameter :: expected_arr_2d(3,2) = &
         transpose(reshape([PETSC_TRUE, PETSC_FALSE, PETSC_FALSE, &
         PETSC_FALSE, PETSC_TRUE, PETSC_TRUE],&
         [2,3]))

    json => fson_parse_mpi(filename)

    call fson_get_mpi(json, "logical", val=val)
    call test%assert(expected, val, "logical")
    call fson_get_mpi(json, "missing_logical", expected, val)
    call test%assert(expected, val, "logical")
    call fson_get_mpi(json, "logical_array", val=arr)
    call test%assert(expected_arr, arr, "logical array")
    deallocate(arr)
    call fson_get_mpi(json, "logical_array_2d", val=arr_2d)
    call test%assert(expected_arr_2d, arr_2d, "2d logical array")
    deallocate(arr_2d)
    call fson_get_mpi(json, "empty_array", val=arr)
    call test%assert(PETSC_FALSE, allocated(arr), "empty array")

    call fson_destroy_mpi(json)

  end subroutine test_fson_mpi_logical

!------------------------------------------------------------------------

  subroutine test_fson_mpi_char(test)

    ! Test fson_get_mpi character routines

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    PetscInt, parameter :: char_len = 12
    character(len = char_len) :: val
    character(len = char_len), parameter :: expected = "etaoinshrdlu"
    character(len = char_len), allocatable :: arr(:)
    character(len = char_len), parameter :: default_arr(2) = ["a", "b"]
    character(len = char_len), parameter :: expected_arr(3) = &
         ["foo        ", "eric       ", "egeszegedre"]

    json => fson_parse_mpi(filename)

    call fson_get_mpi(json, "character", val=val)
    call test%assert(expected, val, "character")
    call fson_get_mpi(json, "missing_character", expected, val)
    call test%assert(expected, val, "character")

    call fson_get_mpi(json, "char_array", default_arr, char_len, arr)
    call test%assert(expected_arr, arr, "character array")

    call fson_destroy_mpi(json)

  end subroutine test_fson_mpi_char

!------------------------------------------------------------------------

  subroutine test_fson_mpi_array_rank(test)

    ! Test fson_mpi_array_rank routine

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    PetscInt :: r

    json => fson_parse_mpi(filename)

    r = fson_mpi_array_rank(json, "missing")
    call test%assert(-1, r, "missing")

    r = fson_mpi_array_rank(json, "real")
    call test%assert(0, r, "scalar")

    r = fson_mpi_array_rank(json, "int_array")
    call test%assert(1, r, "rank 1")

    r = fson_mpi_array_rank(json, "double_array_2d")
    call test%assert(2, r, "rank 2")

    call fson_destroy_mpi(json)

  end subroutine test_fson_mpi_array_rank

!------------------------------------------------------------------------

  subroutine test_fson_get_name_mpi(test)

    ! fson_get_name_mpi()

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json, sub_json
    character(:), allocatable :: name
    PetscInt :: num_items

    json => fson_parse_mpi(str = '{"foo": 1, "bar": 2}')

    num_items = fson_value_count_mpi(json, ".")
    call test%assert(2, num_items, "num items")

    sub_json => fson_value_get_mpi(json, 1)
    name = fson_get_name_mpi(sub_json)
    call test%assert("foo", name, "foo name")

    sub_json => fson_value_get_mpi(json, 2)
    name = fson_get_name_mpi(sub_json)
    call test%assert("bar", name, "bar name")

    call fson_destroy_mpi(json)

  end subroutine test_fson_get_name_mpi

!------------------------------------------------------------------------

  subroutine test_fson_value_next_mpi(test)
    ! fson_value_next_mpi()

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json, arr_json, elt_json
    PetscInt, allocatable :: arr(:)
    PetscInt :: n, i
    PetscInt, parameter :: expected_arr(6) = [3, -1, 4, -1, 5, -9]

    json => fson_parse_mpi(filename)

    call fson_get_mpi(json, "int_array", arr_json)
    n = fson_value_count_mpi(arr_json, ".")
    allocate(arr(n))
    arr = 0
    elt_json => fson_value_children_mpi(arr_json)

    do i = 1, n
       call fson_get_mpi(elt_json, ".", val = arr(i))
       elt_json => fson_value_next_mpi(elt_json)
    end do

    call test%assert(expected_arr, arr, "array")

    deallocate(arr)
    call fson_destroy_mpi(json)
    
  end subroutine test_fson_value_next_mpi

!------------------------------------------------------------------------

end module fson_mpi_test
