module eos_se_test_module

  ! Tests for eos_se module (non-isothermal pure water equation of state)

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
  use eos_se_module
  use unit_test_utils_module, only: transition_compare

  implicit none
  private

  public :: setup, teardown, setup_test
  public :: test_eos_se_fluid_properties, test_eos_se_transition, &
       test_eos_se_errors, test_eos_se_conductivity, &
       test_eos_se_phase_saturations, test_eos_se_check_primary_variables, &
       test_eos_se_convert_fluid

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

  subroutine test_eos_se_fluid_properties(test)

    ! eos_se fluid_properties() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: fluid
    type(rock_type) :: rock
    PetscInt,  parameter :: offset = 1, region = 3, phase_composition = int(b'100')
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal, allocatable:: primary(:), primary2(:)
    type(eos_se_type) :: eos
    type(IAPWS_type) :: thermo
    class(relative_permeability_type), allocatable :: rp
    class(capillary_pressure_type), allocatable :: cp
    type(fson_value), pointer :: json
    character(120) :: json_str = &
         '{"rock": {"relative_permeability": {"type": "linear", "liquid": [0.2, 0.8], "vapour": [0.2, 0.8]}}}'
    PetscErrorCode :: err
    PetscReal, parameter :: zero_phase(8) = 0._dp
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = json_str)
    call thermo%init()
    call eos%init(json, thermo)
    call setup_relative_permeabilities(json, rp)
    call setup_capillary_pressures(json, cp)

    call fluid%init(eos%num_components, eos%num_phases)
    call rock%init()
    allocate(primary(eos%num_primary_variables), primary2(eos%num_primary_variables))
    allocate(fluid_data(fluid%dof))
    fluid_data = 0._dp
    call fluid%assign(fluid_data, offset)

    call rock%assign_relative_permeability(rp)
    call rock%assign_capillary_pressure(cp)

    call properties_case([400._dp, 500._dp], 3, 0, &
         68.98524210481558e6_dp, 500._dp, int(b'100'), &
         zero_phase, zero_phase, &
         [400._dp, 5.293000292993047e-5_dp, 1._dp, 1._dp, 0._dp, &
         2474.9012953721726e3_dp, 2302.4381901101337e3_dp, 1._dp], 'case 1')

    call properties_case([22.0238e6_dp, 1.e-6_dp], 4, 0, &
         22.0238e6_dp, 373.79577673768074_dp, int(b'011'), &
         [355.80571268328060_dp, 4.247613738907617e-05_dp, 1._dp - 1.e-6_dp, 1._dp, 0._dp, &
         2033.7076693254362e3_dp, 1971.80928148363e3_dp, 1._dp], &
         [287.48968435146560_dp, 3.6291847165268826e-05_dp, 1.e-6_dp, 0._dp, 0._dp, &
         2149.0093965774586e3_dp, 2072.402126129886e3_dp, 1._dp], &
         zero_phase, 'case 2')

    call properties_case([470._dp, 370._dp], 3, 0, &
         21.348358089334372e6_dp, 370._dp, int(b'001'), &
         [470._dp, 5.409468698928427e-5_dp, 1._dp, 1._dp, 0._dp, &
         1869.1405399146092e3_dp, 1823.7185014266638e3_dp, 1._dp], &
         zero_phase, zero_phase, 'case 3')

    call properties_case([36.85e6_dp, 470._dp], 2, 0, &
         36.85e6_dp, 470._dp, int(b'100'), &
         zero_phase, zero_phase, &
         [187.41219462500467_dp, 3.448876799376717e-05_dp, 1._dp, 1._dp, 0._dp, &
         2779.57079988821e3_dp, 2582.945387791982e3_dp, 1._dp], &
         'case 4')

    call fluid%destroy()
    call rock%destroy()
    deallocate(primary, primary2, fluid_data)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(rp)
    deallocate(cp)

  contains

    subroutine properties_case(primary, region, expected_err, expected_pressure, &
         expected_temperature, expected_composition, expected_phase1, &
         expected_phase2, expected_phase3, name)

      PetscReal, intent(in) :: primary(eos%num_primary_variables)
      PetscInt, intent(in) :: region
      PetscErrorCode, intent(in) :: expected_err
      PetscReal, intent(in) :: expected_pressure, expected_temperature
      PetscInt, intent(in) :: expected_composition
      PetscReal, dimension(8), intent(in) :: expected_phase1, &
           expected_phase2, expected_phase3
      character(*), intent(in) :: name

      fluid%region = dble(region)
      fluid%permeability_factor = 1._dp
      call eos%fluid_properties(primary, rock, fluid, err)
      call eos%primary_variables(fluid, primary2)

      if (rank == 0) then
         call test%assert(expected_err, err, trim(name) // " error code")
         if (err == 0) then

            call test%assert(expected_pressure, fluid%pressure, &
                 trim(name) // " pressure")
            call test%assert(expected_temperature, fluid%temperature, &
                 trim(name) // " temperature")
            call test%assert(expected_composition, &
                 nint(fluid%phase_composition), trim(name) // " phase composition")

            call test%assert(expected_phase1, fluid%phase(1)%data, trim(name) // " phase 1")
            call test%assert(expected_phase2, fluid%phase(2)%data, trim(name) // " phase 2")
            call test%assert(expected_phase3, fluid%phase(3)%data, trim(name) // " phase 3")

            call test%assert(primary, primary2, trim(name) // " primary")

         end if
      end if

    end subroutine properties_case

  end subroutine test_eos_se_fluid_properties

!------------------------------------------------------------------------

  subroutine test_eos_se_transition(test)

    ! eos_se_transition() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: old_fluid, fluid
    PetscInt,  parameter :: offset = 1, num_primary_variables = 2
    PetscReal, pointer, contiguous :: old_fluid_data(:), fluid_data(:)
    PetscReal :: old_primary(num_primary_variables), primary(num_primary_variables)
    PetscReal :: expected_primary(num_primary_variables), temperature
    PetscInt :: expected_region, expected_err
    PetscBool :: transition, expected_transition
    type(eos_se_type) :: eos
    type(IAPWS_type) :: thermo
    type(fson_value), pointer :: json
    character(2) :: json_str = '{}'
    character(60) :: title
    PetscErrorCode :: err
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal :: d(2)
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

       title = "Region 1 null transition"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [1.e5_dp, 20._dp]
       expected_transition = PETSC_FALSE
       expected_err = 0
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to 4"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [16.647121334271149e5_dp, small]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [20.e5_dp, 210._dp]
       primary = [15.e5_dp, 200._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 null transition"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [1.e5_dp, 120._dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       expected_err = 0
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 4"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [85.621455812056474e5_dp, 1._dp - small]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [84.0e5_dp, 302._dp]
       primary = [86.e5_dp, 299.27215502281706_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 3 subcritical"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [123.304756880772_dp, 360._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [12.e6_dp, 360._dp]
       primary = [18.e6_dp, 360._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 4, T > 350"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [18.666403421371095e6_dp, 1._dp - small]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [12.e6_dp, 360._dp]
       primary = [19.e6_dp, 360._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 3 supercritical"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [215.1893440445802_dp, 390._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [12.e6_dp, 360._dp]
       primary = [25.e6_dp, 390._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 3 Widom delta"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [305.59127267242570_dp, 452.60284095132027_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [30.e6_dp, 470._dp]
       primary = [45.e6_dp, 450._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 2 to 3 near critical point"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [278.40601191048711_dp, 373.99895032601586_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [18.e6_dp, 380._dp]
       primary = [22.1e6_dp, thermo%critical%temperature + 1.e-2_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 null transition"
       old_fluid%region = dble(4)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.e5_dp, 0.5_dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       expected_err = 0
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to 1"
       temperature = 299.27215502281706_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [85.90917681818182e5_dp, 300.02645326107097_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [85.e5_dp, 0.1_dp]
       primary = [86.e5_dp, -0.01_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to 2"
       temperature = 212.38453531849041_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [20.08331325e5_dp, 212.59487472987195_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [20.e5_dp, 0.9_dp]
       primary = [20.1e5_dp, 1.02_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to supercritical"
       old_fluid%region = dble(4)
       old_fluid%temperature = 370._dp
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [thermo%critical%density, &
            (1._dp + 1.e-6_dp) * thermo%critical%temperature]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [21.043367318975246e6_dp, 0.6_dp]
       primary = [22.2e6_dp, 0.8_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to region 3 liquid, T > 350"
       old_fluid%region = dble(4)
       old_fluid%temperature = 356.99181334434775_dp
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [510.39402907666869_dp, 362.91787808947413_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [18.e6_dp, 0.2_dp]
       primary = [20.e6_dp, -0.1_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to region 3 vapour, T > 350"
       old_fluid%region = dble(4)
       old_fluid%temperature = 373.93854042827775_dp
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [315.09773231293144_dp, 373.94040546390397_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [22.062e6_dp, 0.99_dp]
       primary = [22.063e6_dp, 1.01_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 4 to region 2, T > 350"
       old_fluid%region = dble(4)
       old_fluid%temperature = 350.00000000000387_dp
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [16.529147723441228e6_dp, 350.00000000000387_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [16.529164252605481e6_dp, 0.99_dp]
       primary = [16.529164252605481e6_dp, 1.01_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 null transition"
       old_fluid%region = dble(3)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [450._dp, 390._dp]
       expected_transition = PETSC_FALSE
       expected_err = 0
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 to 1"
       old_fluid%region = dble(3)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [23.43972168985525e6_dp, 340._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [670._dp, 360._dp]
       primary = [650._dp, 340._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 to 2, T < 350"
       old_fluid%region = dble(3)
       old_fluid%pressure = 17.856503054949517e6_dp
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [14.500492505102368e6_dp, 340._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [120._dp, 360._dp]
       primary = [91._dp, 340._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 liquid to 4, T < 350"
       old_fluid%region = dble(3)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [15.540148054706076e6_dp, 0.29233698197382446_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [600._dp, 360._dp]
       primary = [450._dp, 345._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 vapour to 4, T < 350"
       old_fluid%region = dble(3)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [14.600181056805235e6_dp, 0.9473529053320567_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [150._dp, 380._dp]
       primary = [120._dp, 340._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 3 to 2, T > 350"
       old_fluid%region = dble(3)
       old_fluid%pressure = 20.18070653002556e6_dp
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [17.403804418775717e6_dp, 380._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [150._dp, 370._dp]
       primary = [90._dp, 380._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to 3"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [633.2633405333486_dp, 360._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [40.e6_dp, 340._dp]
       primary = [35.e6_dp, 360._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to subcritical region 3, T < Tc, P < Pc"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [563.75348104046191_dp, 360._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [40.e6_dp, 340._dp]
       primary = [21.5e6_dp, 360._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to subcritical region 3, T > Tc, P < Pc"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 0
       expected_primary = [0._dp, 0._dp]
       expected_transition = PETSC_TRUE
       expected_err = 1
       old_primary = [30.e6_dp, 340._dp]
       primary = [21.8e6_dp, 380._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to supercritical region 3"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [569.8416683921789_dp, 380._dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [40.e6_dp, 340._dp]
       primary = [35.e6_dp, 380._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to 3 Widom delta"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [475.23690666454377_dp, 453.36460370531074_dp]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [70.e6_dp, 340._dp]
       primary = [60.e6_dp, 465._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to 4 T > 350"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [19.888468913341358e6_dp, small]
       expected_transition = PETSC_TRUE
       expected_err = 0
       old_primary = [30.e6_dp, 340._dp]
       primary = [18.e6_dp, 370._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, expected_err, primary, fluid, transition, err, title)

       title = "Region 1 to region 3 through critical point"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       old_primary = [20.e6_dp, 345._dp]
       d = [thermo%critical%pressure, thermo%critical%temperature] - old_primary
       primary = old_primary + 1.2 * d
       expected_region = 3
       expected_primary = [358.87867533602133_dp, 373.92641126518146_dp]
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

  end subroutine test_eos_se_transition

!------------------------------------------------------------------------

  subroutine test_eos_se_errors(test)

    ! eos_se error handling

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: n = 2
    type(fluid_type) :: fluid
    type(rock_type) :: rock
    PetscInt, parameter :: num_components = 1
    PetscInt,  parameter :: offset = 1
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal, parameter :: data(n, num_components + 1) = reshape([ &
           20.e6_dp, 101.e6_dp, &
           360._dp, 20._dp], [n, num_components + 1])
    PetscInt, parameter :: region(n) = [1, 2]
    PetscReal :: primary(num_components + 1)
    type(eos_se_type) :: eos
    type(IAPWS_type) :: thermo
    class(relative_permeability_type), allocatable :: rp
    class(capillary_pressure_type), allocatable :: cp
    type(fson_value), pointer :: json
    character(120) :: json_str = &
         '{"rock": {"relative_permeability": {"type": "linear", "liquid": [0.2, 0.8], "vapour": [0.2, 0.8]}}}'
    PetscInt :: i
    PetscErrorCode :: err
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = json_str)
    call thermo%init()
    call eos%init(json, thermo)
    call setup_relative_permeabilities(json, rp)
    call setup_capillary_pressures(json, cp)

    call fluid%init(num_components, eos%num_phases)
    call rock%init()
    allocate(fluid_data(fluid%dof))
    fluid_data = 0._dp
    call fluid%assign(fluid_data, offset)

    call rock%assign_relative_permeability(rp)
    call rock%assign_capillary_pressure(cp)

    do i = 1, n
       primary = data(i, :)
       fluid%region = dble(region(i))
       call eos%fluid_properties(primary, rock, fluid, err)
       if (rank == 0) then
          call test%assert(1, err, "error code")
       end if
    end do

    call fluid%destroy()
    call rock%destroy()
    deallocate(fluid_data)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(rp)
    deallocate(cp)

  end subroutine test_eos_se_errors

!------------------------------------------------------------------------

  subroutine test_eos_se_conductivity(test)

    ! Test eos_se heat conductivity

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_se_type) :: eos
    PetscReal, pointer, contiguous :: fluid_data(:), rock_data(:)
    type(fluid_type) :: fluid
    type(rock_type) :: rock
    PetscInt :: offset = 1
    PetscReal :: cond, expected_cond
    PetscReal, parameter :: tol = 1.e-7_dp
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{}')
    call thermo%init()
    call eos%init(json, thermo)
    call fluid%init(eos%num_components, eos%num_phases)
    call rock%init()
    allocate(fluid_data(fluid%dof), rock_data(rock%dof))

    fluid_data = 0._dp
    rock_data = 0._dp
    rock_data(4:5) = [1.5_dp, 1.0_dp]

    call fluid%assign(fluid_data, offset)
    call rock%assign(rock_data, offset)

    if (rank == 0) then

       associate(sl => fluid_data(11))

         sl = 0.0_dp
         expected_cond = 1.0_dp
         cond = eos%conductivity(rock, fluid)
         call test%assert(expected_cond, cond, "sl = 0.0", tol)

         sl = 0.25_dp
         expected_cond = 1.25_dp
         cond = eos%conductivity(rock, fluid)
         call test%assert(expected_cond, cond, "sl = 0.25", tol)

         sl = 0.5_dp
         expected_cond = 1.3535534_dp
         cond = eos%conductivity(rock, fluid)
         call test%assert(expected_cond, cond, "sl = 0.5", tol)

         sl = 0.75_dp
         expected_cond = 1.4330127_dp
         cond = eos%conductivity(rock, fluid)
         call test%assert(expected_cond, cond, "sl = 0.75", tol)

         sl = 1.0_dp
         expected_cond = 1.5_dp
         cond = eos%conductivity(rock, fluid)
         call test%assert(expected_cond, cond, "sl = 1.0", tol)

       end associate

    end if

    call rock%destroy()
    call fluid%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(fluid_data, rock_data)

  end subroutine test_eos_se_conductivity

!------------------------------------------------------------------------

  subroutine test_eos_se_phase_saturations(test)
    ! Test eos_se phase saturations in region 3.

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_se_type) :: eos
    PetscReal, pointer, contiguous :: fluid_data(:)
    type(fluid_type) :: fluid
    PetscInt :: offset = 1
    PetscInt, parameter :: region = 3
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{}')
    call thermo%init()
    call eos%init(json, thermo)
    call fluid%init(eos%num_components, eos%num_phases)
    allocate(fluid_data(fluid%dof))
    fluid_data = 0._dp
    call fluid%assign(fluid_data, offset)
    fluid%region = dble(region)

    if (rank == 0) then

       call saturations_test([700._dp, 360._dp], [1._dp, 0._dp, 0._dp], "case 1")
       call saturations_test([550._dp, 370._dp], [1._dp, 0._dp, 0._dp], "case 2")
       call saturations_test([150._dp, 370._dp], [0._dp, 1._dp, 0._dp], "case 3")
       call saturations_test([150._dp, 380._dp], [0._dp, 1._dp, 0._dp], "case 4")
       call saturations_test([400._dp, 500._dp], [0._dp, 0._dp, 1._dp], "case 5")
       call saturations_test([600._dp, 450._dp], [0._dp, 0._dp, 1._dp], "case 6")

    end if

    call fluid%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(fluid_data)

  contains

    subroutine saturations_test(primary, expected_saturations, name)

      PetscReal, intent(in) :: primary(eos%num_primary_variables)
      PetscReal, intent(in) :: expected_saturations(eos%num_phases)
      character(*), intent(in) :: name
      ! Locals:
      PetscReal :: saturations(eos%num_phases), props(2)
      PetscInt :: i
      PetscErrorCode :: err

      call thermo%region(3)%ptr%properties(primary, props, err)
      fluid%pressure = props(1)
      fluid%temperature = primary(2)
      call fluid%update_phase_composition(thermo)
      call eos%phase_saturations(primary, fluid)

      do i = 1, eos%num_phases
         saturations(i) = fluid%phase(i)%saturation
      end do
      call test%assert(expected_saturations, saturations, name)

    end subroutine saturations_test

  end subroutine test_eos_se_phase_saturations

! ------------------------------------------------------------------------

  subroutine test_eos_se_check_primary_variables(test)
    ! Test eos_se check_primary_variables in region 3.

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_se_type) :: eos
    PetscReal, pointer, contiguous :: fluid_data(:)
    type(fluid_type) :: fluid
    PetscInt :: offset = 1
    PetscInt, parameter :: region = 3
    PetscMPIInt :: rank
    PetscReal, allocatable :: primary(:)
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{}')
    call thermo%init()
    call eos%init(json, thermo)
    call fluid%init(eos%num_components, eos%num_phases)
    allocate(fluid_data(fluid%dof))
    fluid_data = 0._dp
    call fluid%assign(fluid_data, offset)
    fluid%region = dble(region)
    allocate(primary(eos%num_primary_variables))

    if (rank == 0) then

       primary = [550._dp, 360._dp]
       call check_primary_test(primary, 0, "case 1")
       primary = [400._dp, 500._dp]
       call check_primary_test(primary, 0, "case 2")
       primary = [1000._dp, 360._dp]
       call check_primary_test(primary, 1, "case 3")
       primary = [200._dp, 810._dp]
       call check_primary_test(primary, 1, "case 4")

    end if

    call fluid%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(fluid_data, primary)

  contains

    subroutine check_primary_test(primary, expected_err, name)

      PetscReal, intent(in out) :: primary(eos%num_primary_variables)
      PetscErrorCode, intent(in) :: expected_err
      character(*), intent(in) :: name
      ! Locals:
      PetscErrorCode :: err
      PetscBool :: changed

      call eos%check_primary_variables(fluid, primary, changed, err)
      call test%assert(expected_err, err, name)

    end subroutine check_primary_test

  end subroutine test_eos_se_check_primary_variables

!------------------------------------------------------------------------

  subroutine test_eos_se_convert_fluid(test)
    ! Test eos_se convert_fluid ()

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_se_type) :: eos
    PetscReal, pointer, contiguous :: fluid_data(:), rock_data(:)
    type(fluid_type) :: fluid(2)
    type(rock_type) :: rock
    class(relative_permeability_type), allocatable :: rp
    class(capillary_pressure_type), allocatable :: cp
    PetscInt :: offset(2), rock_offset = 1
    PetscMPIInt :: rank
    PetscInt :: ierr, i
    PetscErrorCode :: err

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{}')
    call thermo%init()
    call eos%init(json, thermo)
    call setup_relative_permeabilities(json, rp)
    call setup_capillary_pressures(json, cp)

    do i = 1, 2
       call fluid(i)%init(eos%num_components, eos%num_phases)
    end do
    allocate(fluid_data(2 * fluid(1)%dof))
    fluid_data = 0._dp
    offset = [1, 1 + fluid(1)%dof]
    do i = 1, 2
       call fluid(i)%assign(fluid_data, offset(i))
    end do

    call rock%init()
    allocate(rock_data(rock%dof))
    call rock%assign(rock_data, rock_offset)
    call rock%assign_relative_permeability(rp)
    call rock%assign_capillary_pressure(cp)

    if (rank == 0) then

       call convert_fluid_test( &
            [1.e5_dp, 20._dp], 1, [1.e5_dp, 20._dp], 1, &
            [1, 1], [1._dp, 0._dp, 0._dp], "L-L")

       call convert_fluid_test( &
            [3.e5_dp, 120._dp], 1, [1.e5_dp, 120._dp], 2, &
            [1, 2], [0._dp, 1._dp, 0._dp], "L-V")

       call convert_fluid_test( &
            [3.e5_dp, 120._dp], 1, [1.e5_dp, 0.3_dp], 4, &
            [1, 3], [0.7_dp, 0.3_dp, 0._dp], "L-2P")

       call convert_fluid_test( &
            [600._dp, 400._dp], 3, [600._dp, 450._dp], 3, &
            [4, 4], [0._dp, 0._dp, 1._dp], "SL-SL")

       call convert_fluid_test( &
            [600._dp, 450._dp], 3, [200._dp, 400._dp], 3, &
            [4, 4], [0._dp, 0._dp, 1._dp], "SL-SV")

       call convert_fluid_test( &
            [600._dp, 450._dp], 3, [291._dp, 400._dp], 3, &
            [4, 4], [0._dp, 0._dp, 1._dp], "SL-S2P")

       call convert_fluid_test( &
            [40.e6_dp, 340._dp], 1, [650._dp, 360._dp], 3, &
            [1, 1], [1._dp, 0._dp, 0._dp], "L-3L")

       call convert_fluid_test( &
            [17.e6_dp, 360._dp], 2, [150._dp, 370._dp], 3, &
            [2, 2], [0._dp, 1._dp, 0._dp], "V-3V")

       call convert_fluid_test( &
            [50.e6_dp, 340._dp], 1, [600._dp, 400._dp], 3, &
            [1, 1], [1._dp, 0._dp, 0._dp], "L-SL")

       call convert_fluid_test( &
            [50.e6_dp, 340._dp], 1, [200._dp, 400._dp], 3, &
            [1, 2], [0._dp, 1._dp, 0._dp], "L-SV")

       call convert_fluid_test( &
            [50.e6_dp, 340._dp], 1, [280._dp, 400._dp], 3, &
            [1, 3], [0.34576855753421931_dp, 0.65423144246578069_dp, &
            0._dp], "L-S2P")

    end if

    do i = 1, 2
       call fluid(i)%destroy()
    end do
    call rock%destroy()
    deallocate(rp)
    deallocate(cp)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(fluid_data, rock_data)

  contains

    subroutine convert_fluid_test(primary1, region1, primary2, &
         region2, expected_phase_composition, expected_saturations2, name)

      PetscReal, intent(in) :: primary1(eos%num_primary_variables)
      PetscReal, intent(in) :: primary2(eos%num_primary_variables)
      PetscInt, intent(in) :: region1, region2
      PetscInt, intent(in) :: expected_phase_composition(2)
      PetscReal, intent(in) :: expected_saturations2(eos%num_mobile_phases)
      character(*), intent(in) :: name
      ! Locals:
      PetscReal :: sat2(eos%num_mobile_phases)

      fluid(1)%region = dble(region1)
      call eos%fluid_properties(primary1, rock, fluid(1), err)
      fluid(2)%region = dble(region2)
      call eos%fluid_properties(primary2, rock, fluid(2), err)

      call eos%convert_fluid(fluid(1), fluid(2))

      call test%assert(expected_phase_composition(1), &
           nint(fluid(1)%phase_composition), name // ' phases 1')
      call test%assert(expected_phase_composition(2), &
           nint(fluid(2)%phase_composition), name // ' phases 2')
      sat2 = [fluid(2)%phase(1)%saturation, fluid(2)%phase(2)%saturation, &
           fluid(2)%phase(3)%saturation]
      call test%assert(expected_saturations2, sat2, name // ' saturations 2')

    end subroutine convert_fluid_test

  end subroutine test_eos_se_convert_fluid

!------------------------------------------------------------------------

end module eos_se_test_module
