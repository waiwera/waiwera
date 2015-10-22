module kinds_test

  ! Test for kinds module

  use kinds_module
  use mpi_module
  use fruit

  implicit none
  private

#include <petsc/finclude/petsc.h90>

public :: test_nan

contains
  
!------------------------------------------------------------------------

  subroutine test_nan

    ! Test NaN constant

    PetscBool :: nan

    if (mpi%rank == mpi%output_rank) then

       nan = PetscIsInfOrNanReal(qnan_dp)
       call assert_equals(.true., nan, "NaN")

    end if

  end subroutine test_nan

!------------------------------------------------------------------------

end module kinds_test
