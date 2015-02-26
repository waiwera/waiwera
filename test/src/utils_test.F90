module utils_test

  ! Tests for utils module

  use fruit
  use utils_module

  implicit none
  private 

  public :: test_str_to_upper, test_str_to_lower

contains

!------------------------------------------------------------------------

  subroutine test_str_to_upper

    ! Test str_to_upper() function
    
    integer, parameter :: strlen = 8
    character(len = strlen) :: str, upper_str
    character(len = strlen), parameter :: expected = "ABCABC12"
    
    str = "abcABC12"
    upper_str = str_to_upper(str)
    call assert_equals(expected, upper_str, 'str_to_upper')

  end subroutine test_str_to_upper

!------------------------------------------------------------------------

  subroutine test_str_to_lower

    ! Test str_to_lower() function
    
    integer, parameter :: strlen = 8
    character(len = strlen) :: str, lower_str
    character(len = strlen), parameter :: expected = "abcabc12"
    
    str = "abcABC12"
    lower_str = str_to_lower(str)
    call assert_equals(expected, lower_str, 'str_to_lower')

  end subroutine test_str_to_lower

!------------------------------------------------------------------------


end module utils_test
