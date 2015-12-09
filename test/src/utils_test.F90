module utils_test

  ! Tests for utils module

  use mpi_module
  use fruit
  use utils_module

  implicit none
  private 

#include <petsc/finclude/petscdef.h>

  public :: test_str_to_upper, test_str_to_lower, &
       test_split_filename, test_change_filename_extension

contains

!------------------------------------------------------------------------

  subroutine test_str_to_upper

    ! Test str_to_upper()
    
    PetscInt, parameter :: strlen = 8
    character(len = strlen) :: str, upper_str
    character(len = strlen), parameter :: expected = "ABCABC12"
    
    if (mpi%rank == mpi%output_rank) then
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
    
    if (mpi%rank == mpi%output_rank) then
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

    if (mpi%rank == mpi%output_rank) then
    
       filename = "model.h5"
       expected_base = "model"
       expected_ext = "h5"
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
       expected_ext = ""
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

    if (mpi%rank == mpi%output_rank) then
    
       filename = "model.h5"
       ext = "log"
       expected_filename = "model.log"
       new_filename = change_filename_extension(filename, ext)
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

end module utils_test
