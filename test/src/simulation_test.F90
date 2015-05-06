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
    character(max_filename_length), parameter :: filename = "data/simulation/init/test_init.json"
    character(20), parameter :: expected_title = "Test simulation init"
    character(16), parameter :: expected_thermo = "IAPWS-97"
    character(1), parameter  :: expected_eos = "W"
    PetscInt, parameter :: expected_dim = 3, num_cells = 12, num_primary = 1
    PetscInt :: dim, global_solution_dof
    Vec :: x
    PetscErrorCode :: ierr

    call sim%init(filename)

    call DMGetDimension(sim%mesh%dm, dim, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(sim%mesh%dm, x, ierr); CHKERRQ(ierr)
    call VecGetSize(x, global_solution_dof, ierr); CHKERRQ(ierr)
    call VecDestroy(x, ierr); CHKERRQ(ierr)

    if (mpi%rank == mpi%output_rank) then
       call assert_equals(expected_title, sim%title, "Simulation title")
       call assert_equals(expected_thermo, sim%thermo%name, "Simulation thermodynamics")
       call assert_equals(expected_eos, sim%eos%name, "Simulation EOS")
       call assert_equals(expected_dim, dim, "Simulation mesh dimension")
       call assert_equals(num_cells * num_primary, global_solution_dof, "Simulation mesh dof")
    end if

    call sim%destroy()

  end subroutine test_simulation_init

!------------------------------------------------------------------------

end module simulation_test
