module eos_test

  ! Tests for EOS module

  use kinds_module
  use mpi_module
  use fruit
  use eos_module
  use fluid_module
  use IAPWS_module
  use fson
  use fson_mpi_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>

public :: test_eos_w_fluid_properties, &
     test_eos_we_transition

contains

!------------------------------------------------------------------------

  subroutine test_eos_w_fluid_properties

    ! eos_w fluid_properties() test

    type(fluid_type) :: fluid
    PetscInt, parameter :: num_components = 1, num_phases = 1
    PetscInt,  parameter :: offset = 1, region = 1, phase_composition = b'01'
    PetscReal, allocatable :: fluid_data(:)
    PetscReal :: primary(1)
    type(eos_w_type) :: eos
    type(IAPWS_type) :: thermo
    type(fson_value), pointer :: json
    character(40) :: json_str = '{"eos":{"temperature": 20.0}}'
    PetscReal, parameter :: pressure = 1.e5_dp, tol = 1.e-8_dp
    PetscReal, parameter :: expected_density = 998.20548637769673_dp
    PetscReal, parameter :: expected_internal_energy = 83911.631393167205_dp
    PetscReal, parameter :: expected_specific_enthalpy = 84011.811167136271_dp
    PetscReal, parameter :: expected_viscosity = 1.0015972622270245e-3_dp

    json => fson_parse_mpi(str = json_str)
    call thermo%init()
    call eos%init(json, thermo)

    call fluid%init(num_components, num_phases)
    allocate(fluid_data(fluid%dof()))
    fluid_data = 0._dp
    call fluid%assign(fluid_data, offset)

    primary = pressure
    fluid%region = dble(region)
    call eos%bulk_properties(primary, fluid)
    call eos%phase_composition(fluid)
    call eos%phase_properties(primary, fluid)

    if (mpi%rank == mpi%output_rank) then

       call assert_equals(pressure, fluid%pressure, tol, "Pressure")
       call assert_equals(eos%temperature, fluid%temperature, tol, "Temperature")
       call assert_equals(phase_composition, nint(fluid%phase_composition), "Phase composition")

       call assert_equals(expected_density, fluid%phase(1)%density, &
            tol, "Density")
       call assert_equals(expected_internal_energy, fluid%phase(1)%internal_energy, &
            tol, "Internal energy")
       call assert_equals(expected_specific_enthalpy, fluid%phase(1)%specific_enthalpy, &
            tol, "Specific enthalpy")
       call assert_equals(expected_viscosity, fluid%phase(1)%viscosity, &
            tol, "Viscosity")

       call assert_equals(1._dp, fluid%phase(1)%saturation, tol, "Saturation")
       call assert_equals(1._dp, fluid%phase(1)%relative_permeability, tol, "Relative permeability")
       call assert_equals(1._dp, fluid%phase(1)%mass_fraction(1), tol, "Mass fraction")

    end if

    call fluid%destroy()
    deallocate(fluid_data)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_eos_w_fluid_properties

!------------------------------------------------------------------------

  subroutine transition_compare(expected_primary, expected_region, &
       primary, fluid, message)

    ! Runs asserts to test transition

    PetscReal, intent(in) :: expected_primary(:), primary(:)
    PetscInt, intent(in) :: expected_region
    type(fluid_type), intent(in) :: fluid
    character(60), intent(in) :: message
    ! Locals:
    PetscInt :: n, i
    character(4) :: istr
    PetscReal, parameter :: tol = 1.e-6_dp

    n = size(primary)
    do i = 1, n
       write(istr, '(1x, a1, i2)') '#', i
       call assert_equals(expected_primary(1), primary(1), tol, &
            trim(message) // istr)
    end do
    call assert_equals(expected_region, nint(fluid%region), &
         trim(message) // " region")

  end subroutine transition_compare

!------------------------------------------------------------------------

  subroutine test_eos_we_transition

    ! eos_we_transition() test

    type(fluid_type) :: fluid
    PetscInt, parameter :: num_components = 1, num_phases = 1
    PetscInt,  parameter :: offset = 1
    PetscReal, allocatable :: fluid_data(:)
    PetscReal :: primary(2), expected_primary(2), temperature
    PetscInt :: expected_region
    type(eos_we_type) :: eos
    type(IAPWS_type) :: thermo
    type(fson_value), pointer :: json
    character(2) :: json_str = '{}'
    character(60) :: title
    PetscReal, parameter :: small = 1.e-6_dp

    json => fson_parse_mpi(str = json_str)
    call thermo%init()
    call eos%init(json, thermo)
    call fluid%init(num_components, num_phases)
    allocate(fluid_data(fluid%dof()))
    fluid_data = 0._dp
    call fluid%assign(fluid_data, offset)

    if (mpi%rank == mpi%output_rank) then

       title = "Region 1 null transition"
       expected_region = 1
       expected_primary = [1.e5_dp, 20._dp]
       fluid%region = dble(expected_region)
       primary = expected_primary
       call eos%transition(primary, fluid)
       call transition_compare(expected_primary, expected_region, &
            primary, fluid, title)

       title = "Region 1 to 4"
       expected_region = 4
       expected_primary = [15.546718682698252e5_dp, small]
       primary = [15.e5_dp, 200._dp]
       fluid%region = dble(1)
       call eos%transition(primary, fluid)
       call transition_compare(expected_primary, expected_region, &
            primary, fluid, title)

       title = "Region 2 null transition"
       expected_region = 2
       expected_primary = [1.e5_dp, 120._dp]
       fluid%region = dble(expected_region)
       primary = expected_primary
       call eos%transition(primary, fluid)
       call transition_compare(expected_primary, expected_region, &
            primary, fluid, title)

       title = "Region 2 to 4"
       expected_region = 4
       expected_primary = [85.e5_dp, 1._dp - small]
       fluid%region = dble(2)
       primary = [85.01e5_dp, 299.27215502281706_dp]
       call eos%transition(primary, fluid)
       call transition_compare(expected_primary, expected_region, &
            primary, fluid, title)

       title = "Region 4 null transition"
       expected_region = 4
       expected_primary = [1.e5_dp, 0.5_dp]
       fluid%region = dble(expected_region)
       primary = expected_primary
       call eos%transition(primary, fluid)
       call transition_compare(expected_primary, expected_region, &
            primary, fluid, title)

       title = "Region 4 to 1"
       expected_region = 1
       temperature = 299.27215502281706_dp
       expected_primary = [85.000085e5_dp, temperature]
       fluid%region = dble(4)
       fluid%temperature = temperature
       primary = [85.e5_dp, -0.01_dp]
       call eos%transition(primary, fluid)
       call transition_compare(expected_primary, expected_region, &
            primary, fluid, title)

       title = "Region 4 to 2"
       expected_region = 2
       temperature = 212.38453531849041_dp
       expected_primary = [19.99998e5_dp, temperature]
       fluid%region = dble(4)
       fluid%temperature = temperature
       primary = [20.e5_dp, 1.02_dp]
       call eos%transition(primary, fluid)
       call transition_compare(expected_primary, expected_region, &
            primary, fluid, title)

    end if

    call fluid%destroy()
    deallocate(fluid_data)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_eos_we_transition

!------------------------------------------------------------------------

end module eos_test
