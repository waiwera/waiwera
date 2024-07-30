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
  use IFC67_module
  use IAPWS_module
  use fson
  use fson_mpi_module
  use eos_wse_module
  use unit_test_utils_module, only: transition_compare

  implicit none
  private

  public :: setup, teardown, setup_test
  public :: test_eos_wse_fluid_properties, test_eos_wse_transition

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

  subroutine test_eos_wse_fluid_properties(test)

    ! eos_wse fluid_properties() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: fluid
    type(rock_type) :: rock
    PetscInt,  parameter :: offset = 1, region = 8, phase_composition = int(b'011')
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal, allocatable:: primary(:), primary2(:)
    type(eos_wse_type) :: eos
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
    PetscReal, parameter :: expected_liquid_density = 1106.6518135561505_dp
    PetscReal, parameter :: expected_liquid_saturation = 1._dp - solid_saturation &
         - vapour_saturation
    PetscReal, parameter :: expected_liquid_internal_energy = 855553.6653034494_dp
    PetscReal, parameter :: expected_liquid_viscosity = 0.0003082594394253848_dp
    PetscReal, parameter :: expected_liquid_relative_permeability = 0.14713911448930356_dp
    PetscReal, parameter :: expected_liquid_capillary_pressure = 0._dp
    PetscReal, parameter :: expected_liquid_salt_mass_fraction = 0.3509746512916994_dp
    PetscReal, parameter :: expected_vapour_density = 15.649607177473_dp
    PetscReal, parameter :: expected_vapour_internal_energy = 2658784.4508994194_dp
    PetscReal, parameter :: expected_vapour_viscosity = 1.8139908004858663e-05_dp
    PetscReal, parameter :: expected_vapour_relative_permeability = 0.791942250831361_dp
    PetscReal, parameter :: expected_vapour_capillary_pressure = 0._dp
    PetscReal, parameter :: expected_solid_density = 2099.742571479692_dp
    PetscReal, parameter :: expected_solid_internal_energy = -326148.8593836402_dp
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

    primary = [pressure, vapour_saturation, solid_saturation]
    fluid%region = dble(region)
    call eos%fluid_properties(primary, rock, fluid, err)
    call eos%primary_variables(fluid, primary2)

    if (rank == 0) then

       call test%assert(0, err, "Error code")
       call test%assert(pressure, fluid%pressure, "Pressure")
       call test%assert(temperature, fluid%temperature, "Temperature")
       call test%assert(phase_composition, nint(fluid%phase_composition), "Phase composition")

       call test%assert(expected_liquid_density, fluid%phase(1)%density, &
            "Liquid density")
       call test%assert(expected_liquid_internal_energy, &
            fluid%phase(1)%internal_energy, "Liquid internal energy")
       call test%assert(expected_liquid_viscosity, fluid%phase(1)%viscosity, &
            "Liquid viscosity")
       call test%assert(expected_liquid_saturation, &
            fluid%phase(1)%saturation, "Liquid saturation")
       call test%assert(expected_liquid_relative_permeability, &
            fluid%phase(1)%relative_permeability, &
            "Liquid relative permeability")
       call test%assert(expected_liquid_capillary_pressure, &
            fluid%phase(1)%capillary_pressure, &
            "Liquid capillary pressure")
       call test%assert(1._dp - expected_liquid_salt_mass_fraction, &
            fluid%phase(1)%mass_fraction(1), "Liquid water mass fraction")
       call test%assert(expected_liquid_salt_mass_fraction, &
            fluid%phase(1)%mass_fraction(2), "Liquid salt mass fraction")

       call test%assert(expected_vapour_density, fluid%phase(2)%density, &
            "Vapour density")
       call test%assert(expected_vapour_internal_energy, &
            fluid%phase(2)%internal_energy, "Vapour internal energy")
       call test%assert(expected_vapour_viscosity, fluid%phase(2)%viscosity, &
            "Vapour viscosity")
       call test%assert(vapour_saturation, fluid%phase(2)%saturation, &
            "Vapour saturation")
       call test%assert(expected_vapour_relative_permeability, &
            fluid%phase(2)%relative_permeability, &
            "Vapour relative permeability")
       call test%assert(expected_vapour_capillary_pressure, &
            fluid%phase(2)%capillary_pressure, &
            "Vapour capillary pressure")
       call test%assert(1._dp, fluid%phase(2)%mass_fraction(1), &
            "Vapour water mass fraction")
       call test%assert(0._dp, fluid%phase(2)%mass_fraction(2), &
            "Vapour salt mass fraction")

       call test%assert(solid_saturation, fluid%phase(3)%saturation, &
            "Solid saturation")
       call test%assert(expected_solid_density, fluid%phase(3)%density, &
            "Solid density")
       call test%assert(expected_solid_internal_energy, &
            fluid%phase(3)%internal_energy, "Solid internal energy")
       call test%assert(0._dp, fluid%phase(3)%mass_fraction(1), &
            "Solid water mass fraction")
       call test%assert(1._dp, fluid%phase(3)%mass_fraction(2), &
            "Solid salt mass fraction")

       call test%assert(pressure, primary2(1), "Primary 1")
       call test%assert(vapour_saturation, primary2(2), "Primary 2")
       call test%assert(solid_saturation, primary2(3), "Primary 3")

    end if

    call fluid%destroy()
    call rock%destroy()
    deallocate(primary, primary2, fluid_data)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(rp)
    deallocate(cp)

  end subroutine test_eos_wse_fluid_properties

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

       ! no salt:

       title = "Region 1 null transition, Xs = 0"
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

       title = "Region 1 to 4, Xs = 0"
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

       title = "Region 2 null transition, Xs = 0"
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

       title = "Region 2 to 4, Xs = 0"
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

       title = "Region 4 null transition, Xs = 0"
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

       title = "Region 4 to 1, Xs = 0"
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

       title = "Region 4 to 2, Xs = 0"
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

       ! salt, no halite:

       title = "Region 1 null transition, Xs > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [1.e5_dp, 20._dp, 0.2_dp]
       expected_transition = PETSC_FALSE
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 1 to 4, Xs > 0"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.52428924e+06_dp, small, 4.80568610e-02_dp]
       expected_transition = PETSC_TRUE
       old_primary = [20.e5_dp, 210._dp, 0.01_dp]
       primary = [15.e5_dp, 200._dp, 0.05_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2, Xs > 0"
       ! Dry steam with salt will precipitate out
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 6
       expected_primary = [1.e5_dp, 120._dp, small]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_TRUE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2 to 4, Xs > 0"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [8.501510451467e6_dp, 1._dp - small, 3.03013436e-02_dp]
       expected_transition = PETSC_TRUE
       old_primary = [84.0e5_dp, 302._dp, 0.01_dp]
       primary = [86.e5_dp, 299.27215502281706_dp, 0.05_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 null transition, Xs > 0"
       old_fluid%region = dble(4)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.e5_dp, 0.5_dp, 0.2_dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 to 1, Xs > 0"
       temperature = 299.27215502281706_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [85.90917681818182e5_dp, 301.26248746444287_dp, &
            0.028181818181818_dp]
       expected_transition = PETSC_TRUE
       old_primary = [85.e5_dp, 0.1_dp, 0.01_dp]
       primary = [86.e5_dp, -0.01_dp, 0.03_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       ! halite:

       title = "Region 5 null transition"
       old_fluid%region = dble(5)
       fluid%region = old_fluid%region
       expected_region = 5
       expected_primary = [1.e5_dp, 20._dp, 0.1_dp]
       expected_transition = PETSC_FALSE
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 6 null transition"
       old_fluid%region = dble(6)
       fluid%region = old_fluid%region
       expected_region = 6
       expected_primary = [1.e5_dp, 120._dp, 0.2_dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 8 null transition"
       old_fluid%region = dble(8)
       fluid%region = old_fluid%region
       expected_region = 8
       expected_primary = [1.e5_dp, 0.5_dp, 0.25_dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 1 to 5"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 5
       expected_primary = [20.e5_dp, 210._dp, small]
       expected_transition = PETSC_TRUE
       old_primary = [20.e5_dp, 210._dp, 0.32_dp]
       primary = [20.e5_dp, 210._dp, 0.325_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 5 to 1"
       old_fluid%region = dble(5)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [20.e5_dp, 210._dp, 0.3220677667197454_dp]
       expected_transition = PETSC_TRUE
       old_primary = [20.e5_dp, 210._dp, 0.05_dp]
       primary = [20.e5_dp, 210._dp, -0.01_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2 to 6"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 6
       expected_primary = [60.e5_dp, 302._dp, small]
       expected_transition = PETSC_TRUE
       old_primary = [60.e5_dp, 302._dp, 0.1_dp]
       primary = [60.e5_dp, 302._dp, 0.1_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 6 to 2"
       old_fluid%region = dble(6)
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [55.e5_dp, 302._dp, 0._dp]
       expected_transition = PETSC_TRUE
       old_primary = [55.e5_dp, 302._dp, 0.1_dp]
       primary = [55.e5_dp, 302._dp, -0.05_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 to 8"
       temperature = 200._dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 8
       expected_primary = [10.e5_dp, 0.1_dp, small]
       expected_transition = PETSC_TRUE
       old_primary = [10.e5_dp, 0.1_dp, 0.25_dp]
       primary = [10.e5_dp, 0.1_dp, 0.33_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 8 to 4"
       temperature = 200._dp
       old_fluid%region = dble(8)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.116895574534e6_dp, 0.1_dp, 0.3172414011477263_dp]
       expected_transition = PETSC_TRUE
       old_primary = [1.116895574534e6_dp, 0.1_dp, 0.01_dp]
       primary = [1.116895574534e6_dp, 0.1_dp, -0.01_dp]
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
