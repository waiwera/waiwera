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
  use eos_wse_module
  use eos_wsge_module
  use eos_wsce_module
  use unit_test_utils_module, only: transition_compare

  implicit none
  private

  public :: setup, teardown, setup_test
  public :: test_eos_wsge_fluid_properties, test_eos_wsge_transition

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

  subroutine test_eos_wsge_fluid_properties(test)

    ! eos_wsge fluid_properties() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: fluid_wse, fluid_wsce
    type(rock_type) :: rock
    PetscInt,  parameter :: offset = 1, region = 8, phase_composition = int(b'011')
    PetscReal, pointer, contiguous :: fluid_data_wse(:), fluid_data_wsce(:)
    PetscReal, allocatable:: primary_wse(:), primary2_wse(:)
    PetscReal, allocatable:: primary_wsce(:), primary2_wsce(:)
    type(eos_wse_type) :: eos_wse
    type(eos_wsce_type) :: eos_wsce
    type(IAPWS_type) :: thermo
    class(relative_permeability_type), allocatable :: rp
    class(capillary_pressure_type), allocatable :: cp
    type(fson_value), pointer :: json
    character(120) :: json_str = &
         '{"rock": {"relative_permeability": {"type": "linear", "liquid": [0.35, 1.0], "vapour": [0.0, 0.7]}}}'
    PetscErrorCode :: err
    PetscReal, parameter :: pressure = 33.7726e5_dp
    PetscReal, parameter :: temperature = 261.8535057357576_dp
    PetscReal, parameter :: vapour_saturation = 0.375914_dp
    PetscReal, parameter :: solid_saturation = 0.321895_dp
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = json_str)
    call thermo%init()
    call eos_wse%init(json, thermo)
    call eos_wsce%init(json, thermo)
    call setup_relative_permeabilities(json, rp)
    call setup_capillary_pressures(json, cp)

    call fluid_wse%init(eos_wse%num_components, eos_wse%num_phases)
    call fluid_wsce%init(eos_wsce%num_components, eos_wsce%num_phases)
    call rock%init()
    allocate(primary_wse(eos_wse%num_primary_variables), &
         primary2_wse(eos_wse%num_primary_variables))
    allocate(primary_wsce(eos_wsce%num_primary_variables), &
         primary2_wsce(eos_wsce%num_primary_variables))
    allocate(fluid_data_wse(fluid_wse%dof))
    allocate(fluid_data_wsce(fluid_wsce%dof))
    fluid_data_wse = 0._dp
    fluid_data_wsce = 0._dp
    call fluid_wse%assign(fluid_data_wse, offset)
    call fluid_wsce%assign(fluid_data_wsce, offset)

    call rock%assign_relative_permeability(rp)
    call rock%assign_capillary_pressure(cp)

    ! Check consistency between wse and wsce with zero CO2:
    primary_wse = [pressure, vapour_saturation, solid_saturation]
    primary_wsce = [pressure, vapour_saturation, solid_saturation, 0._dp]

    fluid_wse%region = dble(region)
    fluid_wsce%region = dble(region)
    call eos_wse%bulk_properties(primary_wse, fluid_wse, err)
    call eos_wse%phase_composition(fluid_wse, err)
    call eos_wse%phase_properties(primary_wse, rock, fluid_wse, err)
    call eos_wse%primary_variables(fluid_wse, primary2_wse)
    call eos_wsce%bulk_properties(primary_wsce, fluid_wsce, err)
    call eos_wsce%phase_composition(fluid_wsce, err)
    call eos_wsce%phase_properties(primary_wsce, rock, fluid_wsce, err)
    call eos_wsce%primary_variables(fluid_wsce, primary2_wsce)

    if (rank == 0) then

       call test%assert(fluid_wse%pressure, fluid_wsce%pressure, "Pressure")
       call test%assert(fluid_wse%temperature, fluid_wsce%temperature, "Temperature")
       call test%assert(nint(fluid_wse%phase_composition), &
            nint(fluid_wsce%phase_composition), "Phase composition")

       call test%assert(fluid_wse%phase(1)%density, fluid_wsce%phase(1)%density, &
            "Liquid density")
       call test%assert(fluid_wse%phase(1)%internal_energy, &
            fluid_wsce%phase(1)%internal_energy, "Liquid internal energy")
       call test%assert(fluid_wse%phase(1)%viscosity, fluid_wsce%phase(1)%viscosity, &
            "Liquid viscosity")
       call test%assert(fluid_wse%phase(1)%saturation, &
            fluid_wsce%phase(1)%saturation, "Liquid saturation")
       call test%assert(fluid_wse%phase(1)%relative_permeability, &
            fluid_wsce%phase(1)%relative_permeability, &
            "Liquid relative permeability")
       call test%assert(fluid_wse%phase(1)%capillary_pressure, &
            fluid_wsce%phase(1)%capillary_pressure, &
            "Liquid capillary pressure")
       call test%assert(fluid_wse%phase(1)%mass_fraction(1), &
            fluid_wsce%phase(1)%mass_fraction(1), "Liquid water mass fraction")
       call test%assert(fluid_wse%phase(1)%mass_fraction(2), &
            fluid_wsce%phase(1)%mass_fraction(2), "Liquid salt mass fraction")

       call test%assert(fluid_wse%phase(2)%density, fluid_wsce%phase(2)%density, &
            "Vapour density")
       call test%assert(fluid_wse%phase(2)%internal_energy, &
            fluid_wsce%phase(2)%internal_energy, "Vapour internal energy")
       ! (Vapour viscosity is calculated differently with NCG present)
       call test%assert(fluid_wse%phase(2)%saturation, fluid_wsce%phase(2)%saturation, &
            "Vapour saturation")
       call test%assert(fluid_wse%phase(2)%relative_permeability, &
            fluid_wsce%phase(2)%relative_permeability, &
            "Vapour relative permeability")
       call test%assert(fluid_wse%phase(2)%capillary_pressure, &
            fluid_wsce%phase(2)%capillary_pressure, &
            "Vapour capillary pressure")
       call test%assert(fluid_wse%phase(2)%mass_fraction(1), &
            fluid_wsce%phase(2)%mass_fraction(1), &
            "Vapour water mass fraction")
       call test%assert(fluid_wse%phase(2)%mass_fraction(2), &
            fluid_wsce%phase(2)%mass_fraction(2), &
            "Vapour salt mass fraction")

       call test%assert(fluid_wse%phase(3)%saturation, fluid_wsce%phase(3)%saturation, &
            "Solid saturation")
       call test%assert(fluid_wse%phase(3)%density, fluid_wsce%phase(3)%density, &
            "Solid density")
       call test%assert(fluid_wse%phase(3)%internal_energy, &
            fluid_wsce%phase(3)%internal_energy, "Solid internal energy")
       call test%assert(fluid_wse%phase(3)%mass_fraction(1), &
            fluid_wsce%phase(3)%mass_fraction(1), &
            "Solid water mass fraction")
       call test%assert(fluid_wse%phase(3)%mass_fraction(2), &
            fluid_wsce%phase(3)%mass_fraction(2), &
            "Solid salt mass fraction")

       call test%assert(primary2_wse, primary2_wsce(1:3), "Primary")

    end if

    call fluid_wse%destroy()
    call fluid_wsce%destroy()
    call rock%destroy()
    deallocate(primary_wse, primary2_wse, fluid_data_wse)
    deallocate(primary_wsce, primary2_wsce, fluid_data_wsce)
    call eos_wse%destroy()
    call eos_wsce%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(rp)
    deallocate(cp)

  end subroutine test_eos_wsge_fluid_properties

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
       expected_primary = [8.5621455812056474e6_dp, 1._dp - small, &
            4.2429116241129744e-2_dp, 0._dp]
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
       expected_primary = [8.6621455812056493e6_dp, 1._dp - small, 4.2429116241129744e-2_dp, &
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
