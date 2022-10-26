module eos_wse_test_module

  ! Tests for eos_wse module (non-isothermal salt water equation of state)

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use fluid_module
  use rock_module
  use relative_permeability_module
  use capillary_pressure_module
  use IAPWS_module
  use fson
  use fson_mpi_module
  use eos_wse_module
  use unit_test_utils_module, only: transition_compare

  implicit none
  private

  public :: setup, teardown, setup_test
  public :: test_eos_wse_transition

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

  subroutine setup_test(test)

    class(unit_test_type), intent(in out) :: test

    test%tolerance = 1.e-9

  end subroutine setup_test

!------------------------------------------------------------------------

  subroutine test_eos_wse_transition(test)

    ! eos_wse_transition() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: old_fluid, fluid
    PetscInt,  parameter :: offset = 1, num_primary_variables = 3
    PetscReal, pointer, contiguous :: old_fluid_data(:), fluid_data(:)
    PetscReal :: old_primary(num_primary_variables), primary(num_primary_variables)
    PetscReal :: expected_primary(num_primary_variables), temperature
    PetscInt :: expected_region
    PetscBool :: transition, expected_transition
    type(eos_wse_type) :: eos
    type(IAPWS_type) :: thermo
    type(fson_value), pointer :: json
    character(2) :: json_str = '{}'
    character(60) :: title
    PetscErrorCode :: err
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: small = 1.e-6_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = json_str)
    call thermo%init()
    call eos%init(json, thermo)
    call old_fluid%init(eos%num_components, eos%num_phases)
    call fluid%init(eos%num_components, eos%num_phases)
    allocate(old_fluid_data(old_fluid%dof), fluid_data(fluid%dof))
    old_fluid_data = 0._dp
    fluid_data = 0._dp
    call old_fluid%assign(old_fluid_data, offset)
    call fluid%assign(fluid_data, offset)

    if (rank == 0) then

       title = "Region 1 null transition Xs = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [1.e5_dp, 20._dp, 0._dp]
       expected_transition = PETSC_FALSE
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 1 to 4 Xs = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [16.647121334271149e5_dp, small, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [20.e5_dp, 210._dp, 0._dp]
       primary = [15.e5_dp, 200._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2 null transition Xs = 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [1.e5_dp, 120._dp, 0._dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2 to 4 Xs = 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [85.621455812056474e5_dp, 1._dp - small, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [84.0e5_dp, 302._dp, 0._dp]
       primary = [86.e5_dp, 299.27215502281706_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 null transition Xs = 0"
       old_fluid%region = dble(4)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.e5_dp, 0.5_dp, 0._dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 to 1 Xs = 0"
       temperature = 299.27215502281706_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [85.90917681818182e5_dp, 300.02645326107097_dp, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [85.e5_dp, 0.1_dp, 0._dp]
       primary = [86.e5_dp, -0.01_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 to 2 Xs = 0"
       temperature = 212.38453531849041_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [20.08331325e5_dp, 212.59487472987195_dp, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [20.e5_dp, 0.9_dp, 0._dp]
       primary = [20.1e5_dp, 1.02_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

    end if

    call old_fluid%destroy()
    call fluid%destroy()
    deallocate(old_fluid_data, fluid_data)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_eos_wse_transition

!------------------------------------------------------------------------

end module eos_wse_test_module
