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
  !! Utility functions for string handling, formatting, file names etc.

#include <petsc/finclude/petscsys.h>

  use petscsys

  implicit none
  private

  public :: str_to_upper, str_to_lower, &
       int_str_len, str_array_index, &
       split_filename, change_filename_extension, &
       date_time_str
  
contains

!------------------------------------------------------------------------

  function str_to_upper(strIn) result(strOut)
    !! Converts a string to all upper case.

    character(len=*), intent(in) :: strIn !! Input string
    character(len=len(strIn)) :: strOut !! Output uppercase string
    PetscInt :: i,j

    do i = 1, len(strIn)
       j = iachar(strIn(i:i))
       if (j>= iachar("a") .and. j<=iachar("z") ) then
          strOut(i:i) = achar(iachar(strIn(i:i))-32)
       else
          strOut(i:i) = strIn(i:i)
       end if
    end do

  end function str_to_upper

!------------------------------------------------------------------------

  function str_to_lower(strIn) result(strOut)
    !! Converts a string to all lower case.

    character(len=*), intent(in) :: strIn !! Input string
    character(len=len(strIn)) :: strOut !! Output lowercase string
    PetscInt :: i,j

    do i = 1, len(strIn)
       j = iachar(strIn(i:i))
       if (j>= iachar("A") .and. j<=iachar("Z") ) then
          strOut(i:i) = achar(iachar(strIn(i:i))+32)
       else
          strOut(i:i) = strIn(i:i)
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

end module utils_module
