module eos_wge_test_module

  ! Tests for eos_wge module (non-isothermal water-NCG equation of state)

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use fluid_module
  use rock_module
  use relative_permeability_module
  use IAPWS_module
  use fson
  use fson_mpi_module
  use eos_wge_module
  use unit_test_utils_module, only: transition_compare

  implicit none
  private

  public :: setup, teardown
  public :: test_eos_wge_transition

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

  subroutine test_eos_wge_transition(test)

    ! eos_wge_transition() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: old_fluid, fluid
    PetscInt,  parameter :: offset = 1, num_primary_variables = 3
    PetscReal, pointer, contiguous :: old_fluid_data(:), fluid_data(:)
    PetscReal :: old_primary(num_primary_variables), primary(num_primary_variables)
    PetscReal :: expected_primary(num_primary_variables), temperature
    PetscInt :: expected_region, expected_err
    PetscBool :: transition, expected_transition
    type(eos_wge_type) :: eos
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

       title = "Region 1 null transition, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [1.e5_dp, 20._dp, 0._dp]
       expected_transition = PETSC_FALSE
       expected_err = 0
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 null transition, Pg > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [1.e5_dp, 20._dp, 0.2e5_dp]
       expected_transition = PETSC_FALSE
       expected_err = 0
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to 4, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [16.647121334271149e5_dp, small, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [20.e5_dp, 210._dp, 0._dp]
       primary = [15.e5_dp, 200._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to 4, Pg > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [18.31769706741692e5_dp, small, 1.6705757331457702e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [21.e5_dp, 210._dp, 1.e5_dp]
       primary = [17.e5_dp, 200._dp, 2.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 null transition, Pg = 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [1.e5_dp, 120._dp, 0._dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       expected_err = 0
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 null transition, Pg > 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [1.e5_dp, 120._dp, 0.2e5_dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       expected_err = 0
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 4, Pg = 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [85.621455812056474e5_dp, 1._dp - small, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [84.0e5_dp, 302._dp, 0._dp]
       primary = [86.e5_dp, 299.27215502281706_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 4, Pg > 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [86.810727906028237e5_dp, 1._dp - small, &
            1.1892720939717567e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [86.0e5_dp, 302._dp, 2.e5_dp]
       primary = [87.e5_dp, 299.27215502281706_dp, 1.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 null transition, Pg = 0"
       old_fluid%region = dble(4)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.e5_dp, 0.5_dp, 0._dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       expected_err = 0
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 null transition, Pg > 0"
       old_fluid%region = dble(4)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.e5_dp, 0.5_dp, 0.2e5_dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       expected_err = 0
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to 1, Pg = 0"
       temperature = 299.27215502281706_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [85.909176818181816e5_dp, 300.02645326107097_dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [85.e5_dp, 0.1_dp, 0._dp]
       primary = [86.e5_dp, -0.01_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to 1, Pg > 0"
       temperature = 299.27215502281706_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [87.545540454545449e5_dp, 300.02645326107097_dp, &
            1.6363636363636365e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [88.e5_dp, 0.1_dp, 3.e5_dp]
       primary = [87.5e5_dp, -0.01_dp, 1.5e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to 2, Pg = 0"
       temperature = 212.38453531849041_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [20.08331325e5_dp, 212.59487472987195_dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [20.e5_dp, 0.9_dp, 0._dp]
       primary = [20.1e5_dp, 1.02_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to 2, Pg > 0"
       temperature = 212.38453531849041_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [23.749979916666667e5_dp, 212.59487472987195_dp, &
            3.6666666666666663e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [22.e5_dp, 0.9_dp, 2.e5_dp]
       primary = [24.1e5_dp, 1.02_dp, 4.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

    end if

    call old_fluid%destroy()
    call fluid%destroy()
    deallocate(old_fluid_data, fluid_data)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_eos_wge_transition

!------------------------------------------------------------------------
  
end module eos_wge_test_module
