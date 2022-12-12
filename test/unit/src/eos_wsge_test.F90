module eos_wsge_test_module

  ! Tests for eos_wsge module (non-isothermal salt water / NCG equation of state)

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use fluid_module
  use rock_module
  use relative_permeability_module
  use capillary_pressure_module
  use IFC67_module
  use IAPWS_module
  use fson
  use fson_mpi_module
  use eos_wsge_module
  use unit_test_utils_module, only: transition_compare

  implicit none
  private

  public :: setup, teardown, setup_test
  public :: test_eos_wsge_transition

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

  subroutine test_eos_wsge_transition(test)

    ! eos_wsge_transition() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: old_fluid, fluid
    PetscInt,  parameter :: offset = 1, num_primary_variables = 4
    PetscReal, pointer, contiguous :: old_fluid_data(:), fluid_data(:)
    PetscReal :: old_primary(num_primary_variables), primary(num_primary_variables)
    PetscReal :: expected_primary(num_primary_variables), temperature
    PetscInt :: expected_region
    PetscBool :: transition, expected_transition
    type(eos_wsge_type) :: eos
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

       ! no salt or NCG:

       title = "Region 1 null transition, Xs = 0, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [1.e5_dp, 20._dp, 0._dp, 0._dp]
       expected_transition = PETSC_FALSE
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 1 to 4, Xs = 0, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [16.647121334271149e5_dp, small, 0._dp, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [20.e5_dp, 210._dp, 0._dp, 0._dp]
       primary = [15.e5_dp, 200._dp, 0._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2 null transition, Xs = 0, Pg = 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [1.e5_dp, 120._dp, 0._dp, 0._dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2 to 4, Xs = 0, Pg = 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [85.621455812056474e5_dp, 1._dp - small, 0._dp, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [84.0e5_dp, 302._dp, 0._dp, 0._dp]
       primary = [86.e5_dp, 299.27215502281706_dp, 0._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 null transition, Xs = 0, Pg = 0"
       old_fluid%region = dble(4)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.e5_dp, 0.5_dp, 0._dp, 0._dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 to 1, Xs = 0, Pg = 0"
       temperature = 299.27215502281706_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [85.90917681818182e5_dp, 300.02645326107097_dp, 0._dp, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [85.e5_dp, 0.1_dp, 0._dp, 0._dp]
       primary = [86.e5_dp, -0.01_dp, 0._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 to 2, Xs = 0, Pg = 0"
       temperature = 212.38453531849041_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [20.08331325e5_dp, 212.59487472987195_dp, 0._dp, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [20.e5_dp, 0.9_dp, 0._dp, 0._dp]
       primary = [20.1e5_dp, 1.02_dp, 0._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       ! salt, no halite, no NCG:

       title = "Region 1 null transition, Xs > 0, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [1.e5_dp, 20._dp, 0.2_dp, 0._dp]
       expected_transition = PETSC_FALSE
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 1 to 4, Xs > 0, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.52428924e+06_dp, small, 4.80568610e-02_dp, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [20.e5_dp, 210._dp, 0.01_dp, 0._dp]
       primary = [15.e5_dp, 200._dp, 0.05_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2, Xs > 0, Pg = 0"
       ! Dry steam with salt will precipitate out
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 6
       expected_primary = [1.e5_dp, 120._dp, small, 0._dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_TRUE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2 to 4, Xs > 0, Pg = 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [8.501510451467e6_dp, 1._dp - small, 3.03013436e-02_dp, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [84.0e5_dp, 302._dp, 0.01_dp, 0._dp]
       primary = [86.e5_dp, 299.27215502281706_dp, 0.05_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 null transition, Xs > 0, Pg = 0"
       old_fluid%region = dble(4)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.e5_dp, 0.5_dp, 0.2_dp, 0._dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 to 1, Xs > 0, Pg = 0"
       temperature = 299.27215502281706_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [85.90917681818182e5_dp, 301.26248746444287_dp, &
            0.028181818181818_dp, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [85.e5_dp, 0.1_dp, 0.01_dp, 0._dp]
       primary = [86.e5_dp, -0.01_dp, 0.03_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       ! halite, no NCG:

       title = "Region 5 null transition, Pg = 0"
       old_fluid%region = dble(5)
       fluid%region = old_fluid%region
       expected_region = 5
       expected_primary = [1.e5_dp, 20._dp, 0.1_dp, 0._dp]
       expected_transition = PETSC_FALSE
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 6 null transition, Pg = 0"
       old_fluid%region = dble(6)
       fluid%region = old_fluid%region
       expected_region = 6
       expected_primary = [1.e5_dp, 120._dp, 0.2_dp, 0._dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 8 null transition, Pg = 0"
       old_fluid%region = dble(8)
       fluid%region = old_fluid%region
       expected_region = 8
       expected_primary = [1.e5_dp, 0.5_dp, 0.25_dp, 0._dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 1 to 5, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 5
       expected_primary = [20.e5_dp, 210._dp, small, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [20.e5_dp, 210._dp, 0.32_dp, 0._dp]
       primary = [20.e5_dp, 210._dp, 0.325_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 5 to 1, Pg = 0"
       old_fluid%region = dble(5)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [20.e5_dp, 210._dp, 0.3220677667197454_dp, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [20.e5_dp, 210._dp, 0.05_dp, 0._dp]
       primary = [20.e5_dp, 210._dp, -0.01_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2 to 6, Pg = 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 6
       expected_primary = [60.e5_dp, 302._dp, small, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [60.e5_dp, 302._dp, 0.1_dp, 0._dp]
       primary = [60.e5_dp, 302._dp, 0.1_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 6 to 2, Pg = 0"
       old_fluid%region = dble(6)
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [55.e5_dp, 302._dp, 0._dp, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [55.e5_dp, 302._dp, 0.1_dp, 0._dp]
       primary = [55.e5_dp, 302._dp, -0.05_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 to 8, Pg = 0"
       temperature = 200._dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 8
       expected_primary = [10.e5_dp, 0.1_dp, small, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [10.e5_dp, 0.1_dp, 0.25_dp, 0._dp]
       primary = [10.e5_dp, 0.1_dp, 0.33_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 8 to 4, Pg = 0"
       temperature = 200._dp
       old_fluid%region = dble(8)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.116895574534e6_dp, 0.1_dp, 0.3172414011477263_dp, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [1.116895574534e6_dp, 0.1_dp, 0.01_dp, 0._dp]
       primary = [1.116895574534e6_dp, 0.1_dp, -0.01_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       ! NCG with no salt:

       title = "Region 1 null transition, Xs = 0, Pg > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [1.e5_dp, 20._dp, 0._dp, 0.2e5_dp]
       expected_transition = PETSC_FALSE
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 1 to 4, Xs = 0, Pg > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [18.31769706741692e5_dp, small, 0._dp, 1.6705757331457702e5_dp]
       expected_transition = PETSC_TRUE
       old_primary = [21.e5_dp, 210._dp, 0._dp, 1.e5_dp]
       primary = [17.e5_dp, 200._dp, 0._dp, 2.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2 null transition, Xs = 0, Pg > 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [1.e5_dp, 120._dp, 0._dp, 0.2e5_dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2 to 4, Xs = 0, Pg > 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [86.810727906028237e5_dp, 1._dp - small, 0._dp, &
            1.1892720939717567e5_dp]
       expected_transition = PETSC_TRUE
       old_primary = [86.0e5_dp, 302._dp, 0._dp, 2.e5_dp]
       primary = [87.e5_dp, 299.27215502281706_dp, 0._dp, 1.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 null transition, Xs = 0, Pg > 0"
       old_fluid%region = dble(4)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.e5_dp, 0.5_dp, 0._dp, 0.2e5_dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 to 1, Xs = 0, Pg > 0"
       temperature = 299.27215502281706_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [87.545540454545449e5_dp, 300.02645326107097_dp, &
            0._dp, 1.6363636363636365e5_dp]
       expected_transition = PETSC_TRUE
       old_primary = [88.e5_dp, 0.1_dp, 0._dp, 3.e5_dp]
       primary = [87.5e5_dp, -0.01_dp, 0._dp, 1.5e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 to 2, Xs = 0, Pg > 0"
       temperature = 212.38453531849041_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [23.749979916666667e5_dp, 212.59487472987195_dp, &
            0._dp, 3.6666666666666663e5_dp]
       expected_transition = PETSC_TRUE
       old_primary = [22.e5_dp, 0.9_dp, 0._dp, 2.e5_dp]
       primary = [24.1e5_dp, 1.02_dp, 0._dp, 4.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       ! Salt + NCG:

       title = "Region 1 to 4, Xs > 0, Pg > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.62428924e+06_dp, small, 4.80568610e-02_dp, 1.e5_dp]
       expected_transition = PETSC_TRUE
       old_primary = [21.e5_dp, 210._dp, 0.01_dp, 1.e5_dp]
       primary = [16.e5_dp, 200._dp, 0.05_dp, 1.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2 to 4, Xs > 0, Pg > 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [8.601510451467e6_dp, 1._dp - small, 3.03013436e-02_dp, &
            1.e5_dp]
       expected_transition = PETSC_TRUE
       old_primary = [85.0e5_dp, 302._dp, 0.01_dp, 1.e5_dp]
       primary = [87.e5_dp, 299.27215502281706_dp, 0.05_dp, 1.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 to 1, Xs > 0, Pg > 0"
       temperature = 299.27215502281706_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [86.90917681818182e5_dp, 301.26248746444287_dp, &
            0.028181818181818_dp, 1.e5_dp]
       expected_transition = PETSC_TRUE
       old_primary = [86.e5_dp, 0.1_dp, 0.01_dp, 1.e5_dp]
       primary = [87.e5_dp, -0.01_dp, 0.03_dp, 1.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 1 to 5, Pg > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 5
       expected_primary = [21.e5_dp, 210._dp, small, 1.e5_dp]
       expected_transition = PETSC_TRUE
       old_primary = [21.e5_dp, 210._dp, 0.32_dp, 1.e5_dp]
       primary = [21.e5_dp, 210._dp, 0.325_dp, 1.e5_dp]
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

  end subroutine test_eos_wsge_transition

!------------------------------------------------------------------------

end module eos_wsge_test_module
