module fson_mpi_test

  ! Tests for fson_mpi module

  use kinds_module
  use fruit
  use mpi_module
  use fson_mpi_module

  implicit none
  private 

#include <petsc-finclude/petscsys.h>

  public :: test_fson_mpi_real

contains

!------------------------------------------------------------------------

  subroutine test_fson_mpi_real

    ! Test real number routines


  end subroutine test_fson_mpi_real

!------------------------------------------------------------------------

end module fson_mpi_test
