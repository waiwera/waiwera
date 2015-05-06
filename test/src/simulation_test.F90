module simulation_test

  ! Tests for simulation module

  use kinds_module
  use mpi_module
  use fruit
  use fson
  use simulation_module

  implicit none
  private

#include <petsc-finclude/petsc.h90>

public :: test_simulation_init

PetscReal, parameter :: tol = 1.e-6_dp

contains

!------------------------------------------------------------------------

  subroutine test_simulation_init

    ! test reading simulation title

    type(simulation_type) :: sim
    ! Locals:
    character(max_filename_length) :: filename = "data/simulation/init/test_init.json"
    character(20) :: title = "Test simulation init"

    call sim%init(filename)

    call assert_equals(title, sim%title, "Simulation title")

    call sim%destroy()

  end subroutine test_simulation_init

!------------------------------------------------------------------------

end module simulation_test
