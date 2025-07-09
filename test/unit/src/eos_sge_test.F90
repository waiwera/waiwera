module eos_sge_test_module

  ! Tests for eos_sge module (supercritical water and NCG equation of
  ! state)

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
  use eos_sge_module
  use unit_test_utils_module, only: transition_compare

  implicit none
  private

  public :: setup, teardown, setup_test
  public :: test_eos_sge_transition, test_eos_sge_scale

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

  subroutine test_eos_sge_transition(test)

    ! eos_sge_transition() test

    use eos_wge_module, only: eos_wge_type

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: old_fluid, fluid
    PetscInt,  parameter :: offset = 1, num_primary_variables = 3
    PetscReal, pointer, contiguous :: old_fluid_data(:), fluid_data(:)
    PetscReal :: old_primary(num_primary_variables), primary(num_primary_variables)
    PetscReal :: expected_primary(num_primary_variables), temperature
    PetscReal :: d(num_primary_variables), Pg
    PetscInt :: expected_region, expected_err
    PetscBool :: transition, expected_transition
    type(eos_sge_type) :: eos
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
    allocate(eos_wge_type :: eos%eos_wge)
    call eos%eos_wge%init(json, thermo)
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

       title = "Region 2 to 3 subcritical, Pg = 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [123.304756880772_dp, 360._dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [12.e6_dp, 360._dp, 0._dp]
       primary = [18.e6_dp, 360._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 3 subcritical, Pg > 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [123.304756880772_dp, 360._dp, 4.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [12.2e6_dp, 360._dp, 2.e5_dp]
       primary = [18.4e6_dp, 360._dp, 4.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 4, T > 350, Pg = 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [18.666403421371095e6_dp, 1._dp - small, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [12.e6_dp, 360._dp, 0._dp]
       primary = [19.e6_dp, 360._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 4, T > 350, Pg > 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [18.866403421371095e6_dp, 1._dp - small, 2.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [12.2e6_dp, 360._dp, 2.e5_dp]
       primary = [19.2e6_dp, 360._dp, 2.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 3 supercritical, Pg = 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [215.1893440445802_dp, 390._dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [12.e6_dp, 360._dp, 0._dp]
       primary = [25.e6_dp, 390._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 3 supercritical, Pg > 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [215.1893440445802_dp, 390._dp, 3.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [12.2e6_dp, 360._dp, 2.e5_dp]
       primary = [25.3e6_dp, 390._dp, 3.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 3 Widom delta, Pg = 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [305.59127267242570_dp, 452.60284095132027_dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [30.e6_dp, 470._dp, 0._dp]
       primary = [45.e6_dp, 450._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 3 Widom delta, Pg > 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [305.59127267242570_dp, 452.60284095132027_dp, &
            3.739715904867971e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [30.2e6_dp, 470._dp, 2.e5_dp]
       primary = [45.4e6_dp, 450._dp, 4.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 3 near critical point, Pg = 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [278.40601191048711_dp, 373.99895032601586_dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [18.e6_dp, 380._dp, 0._dp]
       primary = [22.1e6_dp, thermo%critical%temperature + 1.e-2_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 3 near critical point, Pg > 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [278.40601191048711_dp, 373.99895032601586_dp, 4.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [18.4e6_dp, 380._dp, 4.e5_dp]
       primary = [22.5e6_dp, thermo%critical%temperature + 1.e-2_dp, 4.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to supercritical, Pg = 0"
       old_fluid%region = dble(4)
       old_fluid%temperature = 370._dp
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [263.32515822729528_dp, 374.59899361383282_dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [21.043367318975246e6_dp, 0.6_dp, 0._dp]
       primary = [22.2e6_dp, 0.8_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to supercritical, Pg > 0"
       old_fluid%region = dble(4)
       old_fluid%temperature = 370._dp
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [263.32515822729528_dp, 374.59899361383282_dp, 3.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [21.243367318975246e6_dp, 0.6_dp, 2.e5_dp]
       primary = [22.5e6_dp, 0.8_dp, 3.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to region 3 liquid, T > 350, Pg = 0"
       old_fluid%region = dble(4)
       old_fluid%temperature = 356.99181334434775_dp
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [510.49352557031114_dp, 362.91787808947413_dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [18.e6_dp, 0.2_dp, 0._dp]
       primary = [20.e6_dp, -0.1_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to region 3 liquid, T > 350, Pg > 0"
       old_fluid%region = dble(4)
       old_fluid%temperature = 356.99181334434775_dp
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [510.49352557031114_dp, 362.91787808947413_dp, 5.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [18.5e6_dp, 0.2_dp, 5.e5_dp]
       primary = [20.5e6_dp, -0.1_dp, 5.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to region 3 vapour, T > 350, Pg = 0"
       old_fluid%region = dble(4)
       old_fluid%temperature = 373.93854042827775_dp
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [314.99804745225049_dp, 373.94040546390397_dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [22.062e6_dp, 0.99_dp, 0._dp]
       primary = [22.063e6_dp, 1.01_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to region 3 vapour, T > 350, Pg > 0"
       old_fluid%region = dble(4)
       old_fluid%temperature = 373.93854042827775_dp
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [314.99804745225049_dp, 373.94040546390397_dp, 3.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [22.362e6_dp, 0.99_dp, 3.e5_dp]
       primary = [22.363e6_dp, 1.01_dp, 3.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to region 2, T > 350, Pg = 0"
       old_fluid%region = dble(4)
       old_fluid%temperature = 350.00000000000387_dp
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [16.529147723441228e6_dp, 350.00000000000387_dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [16.529164252605481e6_dp, 0.99_dp, 0._dp]
       primary = [16.529164252605481e6_dp, 1.01_dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to region 2, T > 350, Pg > 0"
       old_fluid%region = dble(4)
       old_fluid%temperature = 350.00000000000387_dp
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [17.129147723441228e6_dp, 350.00000000000387_dp, 6.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [17.129164252605481e6_dp, 0.99_dp, 6.e5_dp]
       primary = [17.129164252605481e6_dp, 1.01_dp, 6.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 null transition, Pg = 0"
       old_fluid%region = dble(3)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [450._dp, 390._dp, 0._dp]
       expected_transition = PETSC_FALSE
       expected_err = 0
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 null transition, Pg > 0"
       old_fluid%region = dble(3)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [450._dp, 390._dp, 0.2e5_dp]
       expected_transition = PETSC_FALSE
       expected_err = 0
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 to 1, Pg = 0"
       old_fluid%region = dble(3)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [23.43972168985525e6_dp, 340._dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [670._dp, 360._dp, 0._dp]
       primary = [650._dp, 340._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 to 1, Pg > 0"
       old_fluid%region = dble(3)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [23.93972168985525e6_dp, 340._dp, 5.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [670._dp, 360._dp, 4.e5_dp]
       primary = [650._dp, 340._dp, 5.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 to 2, T < 350, Pg = 0"
       old_fluid%region = dble(3)
       old_fluid%pressure = 17.856503054949517e6_dp
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [14.500492505102368e6_dp, 340._dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [120._dp, 360._dp, 0._dp]
       primary = [91._dp, 340._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 to 2, T < 350, Pg > 0"
       old_fluid%region = dble(3)
       old_fluid%pressure = 17.856503054949517e6_dp
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [14.700492505102368e6_dp, 340._dp, 2.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [120._dp, 360._dp, 3.e5_dp]
       primary = [91._dp, 340._dp, 2.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 liquid to 4, T < 350, Pg = 0"
       old_fluid%region = dble(3)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [15.540148054706076e6_dp, 0.29233698197382446_dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [600._dp, 360._dp, 0._dp]
       primary = [450._dp, 345._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 liquid to 4, T < 350, Pg > 0"
       old_fluid%region = dble(3)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [16.040148054706076e6_dp, 0.29233698197382446_dp, 5.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [600._dp, 360._dp, 4.e5_dp]
       primary = [450._dp, 345._dp, 5.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 vapour to 4, T < 350, Pg = 0"
       old_fluid%region = dble(3)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [14.600181056805235e6_dp, 0.9473529053320567_dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [150._dp, 380._dp, 0._dp]
       primary = [120._dp, 340._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 vapour to 4, T < 350, Pg > 0"
       old_fluid%region = dble(3)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [15.200181056805235e6_dp, 0.9473529053320567_dp, 6.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [150._dp, 380._dp, 7.e5_dp]
       primary = [120._dp, 340._dp, 6.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 to 2, T > 350, Pg = 0"
       old_fluid%region = dble(3)
       old_fluid%pressure = 20.18070653002556e6_dp
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [17.403804418775717e6_dp, 380._dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [150._dp, 370._dp, 0._dp]
       primary = [90._dp, 380._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 to 2, T > 350, Pg > 0"
       old_fluid%region = dble(3)
       old_fluid%pressure = 20.18070653002556e6_dp
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [17.603804418775717e6_dp, 380._dp, 2.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [150._dp, 370._dp, 4.e5_dp]
       primary = [90._dp, 380._dp, 2.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to 3, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [633.2633405333486_dp, 360._dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [40.e6_dp, 340._dp, 0._dp]
       primary = [35.e6_dp, 360._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to 3, Pg > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [633.2633405333486_dp, 360._dp, 2.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [40.5e6_dp, 340._dp, 5.e5_dp]
       primary = [35.2e6_dp, 360._dp, 2.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to subcritical region 3, T < Tc, P < Pc, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [563.75348104046191_dp, 360._dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [40.e6_dp, 340._dp, 0._dp]
       primary = [21.5e6_dp, 360._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to subcritical region 3, T < Tc, P < Pc, Pg > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [563.75348104046191_dp, 360._dp, 6.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [40.7e6_dp, 340._dp, 7.e5_dp]
       primary = [22.1e6_dp, 360._dp, 6.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to subcritical region 3, T > Tc, P < Pc, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [471.20844415025863_dp, 373.946_dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [30.e6_dp, 340._dp, 0._dp]
       primary = [21.8e6_dp, 380._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to subcritical region 3, T > Tc, P < Pc, Pg > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [471.20844415025863_dp, 373.946_dp, 4.84865e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [30.4e6_dp, 340._dp, 4.e5_dp]
       primary = [22.3e6_dp, 380._dp, 5.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to supercritical region 3, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [569.8416683921789_dp, 380._dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [40.e6_dp, 340._dp, 0._dp]
       primary = [35.e6_dp, 380._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to supercritical region 3, Pg > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [569.8416683921789_dp, 380._dp, 6.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [40.3e6_dp, 340._dp, 3.e5_dp]
       primary = [35.6e6_dp, 380._dp, 6.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to 3 Widom delta, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [475.23690666454377_dp, 453.36460370531074_dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [70.e6_dp, 340._dp, 0._dp]
       primary = [60.e6_dp, 465._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to 3 Widom delta, Pg > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [475.23690666454377_dp, 453.36460370531074_dp, &
            4.7207504889274575e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [70.2e6_dp, 340._dp, 2.e5_dp]
       primary = [60.5e6_dp, 465._dp, 5.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to 4 T > 350, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [19.888468913341358e6_dp, small, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [30.e6_dp, 340._dp, 0._dp]
       primary = [18.e6_dp, 370._dp, 0._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to 4 T > 350, Pg > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [20.388468913341358e6_dp, small, 5.e5_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [30.5e6_dp, 340._dp, 5.e5_dp]
       primary = [18.5e6_dp, 370._dp, 5.e5_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to region 3 through critical point, Pg = 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       old_primary = [20.e6_dp, 345._dp, 0._dp]
       d = [thermo%critical%pressure, thermo%critical%temperature, 0._dp] - old_primary
       primary = old_primary + 1.2 * d
       expected_region = 3
       expected_primary = [358.87867533602133_dp, 373.92641126518146_dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to region 3 through critical point, Pg > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       Pg = 3.e5_dp
       old_primary = [20.e6_dp, 345._dp, Pg]
       d(1:2) = [thermo%critical%pressure, thermo%critical%temperature] - old_primary(1:2)
       primary(1:2) = old_primary(1:2) + 1.2 * d(1:2)
       old_primary(1) = old_primary(1) + Pg
       primary(1) = primary(1) + Pg
       primary(3) = Pg
       expected_region = 3
       expected_primary = [358.87867533602133_dp, 373.92641126518146_dp, Pg]
       expected_transition = PETSC_TRUE
       expected_err = 0
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

  end subroutine test_eos_sge_transition

!------------------------------------------------------------------------

  subroutine test_eos_sge_scale(test)

    ! eos_sge_scale(), eos_sge_unscale() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt,  parameter :: num_primary_variables = 3
    type(eos_sge_type) :: eos
    type(IAPWS_type) :: thermo
    type(fson_value), pointer :: json
    character(2) :: json_str = '{}'
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = json_str)
    call thermo%init()
    call eos%init(json, thermo)

    if (rank == 0) then

       call scale_test('region 1', [1.e5_dp, 20._dp, 0.1e5_dp], 1, &
            [0.1_dp, 0.2_dp, 0.1_dp])
       call scale_test('region 2', [2.e5_dp, 300._dp, 0.4e5_dp], 2, &
            [0.2_dp, 3._dp, 0.2_dp])
       call scale_test('region 3', [500._dp, 480._dp, 2.e5_dp], 3, &
            [1.55279503106_dp, 4.8_dp, 0.2_dp])
       call scale_test('region 4', [100.e5_dp, 0.4_dp, 10.e5_dp], 4, &
            [10._dp, 0.4_dp, 0.1_dp])

    end if

    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  contains

    subroutine scale_test(title, primary, region, expected)

      character(*), intent(in) :: title
      PetscReal, intent(in) :: primary(num_primary_variables)
      PetscInt, intent(in) :: region
      PetscReal, intent(in) :: expected(num_primary_variables)
      ! Locals:
      PetscReal :: scaled(num_primary_variables), unscaled(num_primary_variables)

      scaled = eos%scale(primary, region)
      call test%assert(expected, scaled, trim(title) // ' scaled')
      unscaled = eos%unscale(scaled, region)
      call test%assert(primary, unscaled, trim(title) // ' unscaled')

    end subroutine scale_test

  end subroutine test_eos_sge_scale

end module eos_sge_test_module
