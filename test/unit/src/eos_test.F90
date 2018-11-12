module eos_test
  !! Tests for eos module.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use fson
  use fson_mpi_module
  use eos_module, only: max_component_name_length
  use eos_wge_module

  implicit none
  private

  public :: setup, teardown
  public :: test_eos_component_index

contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown()

    PetscErrorCode :: ierr

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

  subroutine test_eos_component_index(test)
    ! eos component_index() test

    use IAPWS_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IAPWS_type) :: thermo
    type(eos_wge_type) :: eos
    type(fson_value), pointer :: json
    character(max_component_name_length) :: name
    PetscMPIInt :: rank
    PetscInt :: ierr

    json => fson_parse_mpi(str = "{}")
    call thermo%init()
    call eos%init(json, thermo)

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       name = "water"
       call test%assert(1, eos%component_index(name), name)

       name = "Water"
       call test%assert(1, eos%component_index(name), name)

       name = "gas"
       call test%assert(2, eos%component_index(name), name)

       name = "energy"
       call test%assert(eos%num_primary_variables, &
            eos%component_index(name), name)

       name = "fred"
       call test%assert(-1, eos%component_index(name), name)

    end if

    call eos%destroy()
    call thermo%destroy()

  end subroutine test_eos_component_index

!------------------------------------------------------------------------

end module eos_test
