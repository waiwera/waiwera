module eos_we_test_module

  ! Tests for eos_we module (non-isothermal pure water equation of state)

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
  use eos_we_module
  use unit_test_utils_module, only: transition_compare

  implicit none
  private

  public :: setup, teardown, setup_test
  public :: test_eos_we_fluid_properties, test_eos_we_transition, &
       test_eos_we_errors, test_eos_we_conductivity

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

  subroutine test_eos_we_fluid_properties(test)

    ! eos_we fluid_properties() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: fluid
    type(rock_type) :: rock
    PetscInt,  parameter :: offset = 1, region = 4, phase_composition = b'011'
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal, allocatable:: primary(:), primary2(:)
    type(eos_we_type) :: eos
    type(IAPWS_type) :: thermo
    class(relative_permeability_type), allocatable :: rp
    class(capillary_pressure_type), allocatable :: cp
    type(fson_value), pointer :: json
    character(120) :: json_str = &
         '{"rock": {"relative_permeability": {"type": "linear", "liquid": [0.2, 0.8], "vapour": [0.2, 0.8]}}}'
    PetscErrorCode :: err
    PetscReal, parameter :: temperature = 230._dp
    PetscReal, parameter :: pressure = 27.967924557686445e5_dp
    PetscReal, parameter :: vapour_saturation = 0.25
    PetscReal, parameter :: expected_liquid_density = 827.12247049977032_dp
    PetscReal, parameter :: expected_liquid_internal_energy = 986828.18916209263_dp
    PetscReal, parameter :: expected_liquid_specific_enthalpy = 990209.54144729744_dp
    PetscReal, parameter :: expected_liquid_viscosity = 1.1619412513757267e-4_dp
    PetscReal, parameter :: expected_liquid_relative_permeability = 11._dp / 12._dp
    PetscReal, parameter :: expected_liquid_capillary_pressure = 0._dp
    PetscReal, parameter :: expected_vapour_density = 13.984012253728331_dp
    PetscReal, parameter :: expected_vapour_internal_energy = 2603010.010356456_dp
    PetscReal, parameter :: expected_vapour_specific_enthalpy = 2803009.2956133024_dp
    PetscReal, parameter :: expected_vapour_viscosity = 1.6704837258831552e-5_dp
    PetscReal, parameter :: expected_vapour_relative_permeability = 1._dp / 12._dp
    PetscReal, parameter :: expected_vapour_capillary_pressure = 0._dp
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = json_str)
    call thermo%init()
    call eos%init(json, thermo)
    call setup_relative_permeabilities(json, rp, err = err)
    call setup_capillary_pressures(json, cp, err = err)

    call fluid%init(eos%num_components, eos%num_phases)
    call rock%init()
    allocate(primary(eos%num_primary_variables), primary2(eos%num_primary_variables))
    allocate(fluid_data(fluid%dof))
    fluid_data = 0._dp
    call fluid%assign(fluid_data, offset)

    call rock%assign_relative_permeability(rp)
    call rock%assign_capillary_pressure(cp)

    primary = [pressure, vapour_saturation]
    fluid%region = dble(region)
    call eos%bulk_properties(primary, fluid, err)
    call eos%phase_composition(fluid, err)
    call eos%phase_properties(primary, rock, fluid, err)
    call eos%primary_variables(fluid, primary2)

    if (rank == 0) then

       call test%assert(pressure, fluid%pressure, "Pressure")
       call test%assert(temperature, fluid%temperature, "Temperature")
       call test%assert(phase_composition, nint(fluid%phase_composition), "Phase composition")

       call test%assert(expected_liquid_density, fluid%phase(1)%density, &
            "Liquid density")
       call test%assert(expected_liquid_internal_energy, &
            fluid%phase(1)%internal_energy, "Liquid internal energy")
       call test%assert(expected_liquid_specific_enthalpy, &
            fluid%phase(1)%specific_enthalpy, "Liquid specific enthalpy")
       call test%assert(expected_liquid_viscosity, fluid%phase(1)%viscosity, &
            "Liquid viscosity")
       call test%assert(1._dp - vapour_saturation, &
            fluid%phase(1)%saturation, "Liquid saturation")
       call test%assert(expected_liquid_relative_permeability, &
            fluid%phase(1)%relative_permeability, &
            "Liquid relative permeability")
       call test%assert(expected_liquid_capillary_pressure, &
            fluid%phase(1)%capillary_pressure, &
            "Liquid capillary pressure")
       call test%assert(1._dp, fluid%phase(1)%mass_fraction(1), &
            "Liquid mass fraction")

       call test%assert(expected_vapour_density, fluid%phase(2)%density, &
            "Vapour density")
       call test%assert(expected_vapour_internal_energy, &
            fluid%phase(2)%internal_energy, "Vapour internal energy")
       call test%assert(expected_vapour_specific_enthalpy, &
            fluid%phase(2)%specific_enthalpy, "Vapour specific enthalpy")
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
            "Vapour mass fraction")

       call test%assert(pressure, primary2(1), "Primary 1")
       call test%assert(vapour_saturation, primary2(2), "Primary 2")

    end if

    call fluid%destroy()
    call rock%destroy()
    deallocate(primary, primary2, fluid_data)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(rp)
    deallocate(cp)

  end subroutine test_eos_we_fluid_properties

!------------------------------------------------------------------------

  subroutine test_eos_we_transition(test)

    ! eos_we_transition() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: old_fluid, fluid
    PetscInt,  parameter :: offset = 1, num_primary_variables = 2
    PetscReal, pointer, contiguous :: old_fluid_data(:), fluid_data(:)
    PetscReal :: old_primary(num_primary_variables), primary(num_primary_variables)
    PetscReal :: expected_primary(num_primary_variables), temperature
    PetscInt :: expected_region
    PetscBool :: transition, expected_transition
    type(eos_we_type) :: eos
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

       title = "Region 1 null transition"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [1.e5_dp, 20._dp]
       expected_transition = PETSC_FALSE
       old_primary = expected_primary
       primary = expected_primary
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 1 to 4"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [16.647121334271149e5_dp, small]
       expected_transition = PETSC_TRUE
       old_primary = [20.e5_dp, 210._dp]
       primary = [15.e5_dp, 200._dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2 null transition"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [1.e5_dp, 120._dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 2 to 4"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [85.621455812056474e5_dp, 1._dp - small]
       expected_transition = PETSC_TRUE
       old_primary = [84.0e5_dp, 302._dp]
       primary = [86.e5_dp, 299.27215502281706_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 null transition"
       old_fluid%region = dble(4)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.e5_dp, 0.5_dp]
       old_primary = expected_primary
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 to 1"
       temperature = 299.27215502281706_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [85.90917681818182e5_dp, 300.02645326107097_dp]
       expected_transition = PETSC_TRUE
       old_primary = [85.e5_dp, 0.1_dp]
       primary = [86.e5_dp, -0.01_dp]
       call eos%transition(old_primary, primary, old_fluid, fluid, transition, err)
       call transition_compare(test, expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, err, title)

       title = "Region 4 to 2"
       temperature = 212.38453531849041_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [20.08331325e5_dp, 212.59487472987195_dp]
       expected_transition = PETSC_TRUE
       old_primary = [20.e5_dp, 0.9_dp]
       primary = [20.1e5_dp, 1.02_dp]
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

  end subroutine test_eos_we_transition

!------------------------------------------------------------------------

  subroutine test_eos_we_errors(test)

    ! eos_we error handling

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
    type(eos_we_type) :: eos
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
    call setup_relative_permeabilities(json, rp, err = err)
    call setup_capillary_pressures(json, cp, err = err)

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
       call eos%bulk_properties(primary, fluid, err)
       call eos%phase_composition(fluid, err)
       call eos%phase_properties(primary, rock, fluid, err)
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

  end subroutine test_eos_we_errors

!------------------------------------------------------------------------

  subroutine test_eos_we_conductivity(test)

    ! Test eos_we heat conductivity

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
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

       associate(sl => fluid_data(8))

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

  end subroutine test_eos_we_conductivity

!------------------------------------------------------------------------
  
end module eos_we_test_module
