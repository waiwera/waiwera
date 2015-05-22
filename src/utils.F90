module utils_module
  !! Various utilities

  implicit none
  private

#include <petsc-finclude/petscdef.h>

  public :: str_to_upper, str_to_lower
  
contains

!------------------------------------------------------------------------

  function str_to_upper(strIn) result(strOut)
    !! Converts a string to all upper case.

    character(len=*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut
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

    character(len=*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut
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

end module utils_module
