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
       test_eos_se_phase_saturations, test_eos_se_check_primary_variables

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
    PetscReal, parameter :: temperature = 500._dp
    PetscReal, parameter :: density = 400._dp
    PetscReal, parameter :: expected_pressure = 68.98524210481558e6_dp
    PetscReal, parameter :: expected_internal_energy = 2302.4381901101337e3_dp
    PetscReal, parameter :: expected_specific_enthalpy = 2474.9012953721726e3_dp
    PetscReal, parameter :: expected_viscosity = 5.293000292993047e-5_dp
    PetscReal, parameter :: expected_relative_permeability = 1._dp
    PetscReal, parameter :: expected_capillary_pressure = 0._dp
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

    primary = [density, temperature]
    fluid%region = dble(region)
    fluid%permeability_factor = 1._dp
    call eos%fluid_properties(primary, rock, fluid, err)
    call eos%primary_variables(fluid, primary2)

    if (rank == 0) then

       call test%assert(0, err, "Error code")
       call test%assert(expected_pressure, fluid%pressure, "Pressure")
       call test%assert(temperature, fluid%temperature, "Temperature")
       call test%assert(phase_composition, nint(fluid%phase_composition), "Phase composition")

       call test%assert(density, fluid%phase(3)%density, "density")
       call test%assert(expected_internal_energy, &
            fluid%phase(3)%internal_energy, "internal energy")
       call test%assert(expected_specific_enthalpy, &
            fluid%phase(3)%specific_enthalpy, "specific enthalpy")
       call test%assert(expected_viscosity, fluid%phase(3)%viscosity, "viscosity")
       call test%assert(1._dp, fluid%phase(3)%saturation, "saturation")
       call test%assert(expected_relative_permeability, &
            fluid%phase(3)%relative_permeability, "relative permeability")
       call test%assert(expected_capillary_pressure, &
            fluid%phase(3)%capillary_pressure, "capillary pressure")
       call test%assert(1._dp, fluid%phase(3)%mass_fraction(1), &
            "mass fraction")

       call test%assert(0._dp, fluid%phase(1)%density, "liquid density")
       call test%assert(0._dp, fluid%phase(2)%density, "vapour density")

       call test%assert(density, primary2(1), "Primary 1")
       call test%assert(temperature, primary2(2), "Primary 2")

    end if

    call fluid%destroy()
    call rock%destroy()
    deallocate(primary, primary2, fluid_data)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(rp)
    deallocate(cp)

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
       expected_primary = [305.66012571460158_dp, 452.59804760149916_dp]
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
       expected_region = 0
       expected_primary = [0._dp, 0._dp]
       expected_transition = PETSC_FALSE
       expected_err = 1
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

       title = "Region 1 to region 3 Widom delta"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 3
       expected_primary = [475.1891651872765_dp, 453.37857530542698_dp]
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
       expected_region = 0
       expected_primary = [0._dp, 0._dp]
       expected_transition = PETSC_FALSE
       expected_err = 1
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

       associate(sl => fluid_data(9))

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

end module eos_se_test_module
