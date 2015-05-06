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
    character(16) :: thermo = "IAPWS-97"
    character(1)  :: eos = "W"

    call sim%init(filename)

    if (mpi%rank == mpi%output_rank) then
       call assert_equals(title, sim%title, "Simulation title")
       call assert_equals(thermo, sim%thermo%name, "Simulation thermodynamics")
       call assert_equals(eos, sim%eos%name, "Simulation EOS")
    end if

    call sim%destroy()

  end subroutine test_simulation_init

!------------------------------------------------------------------------

end module simulation_test
