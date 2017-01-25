module eos_we_test_module

  ! Tests for eos_we module (non-isothermal pure water equation of state)

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fruit
  use fluid_module
  use rock_module
  use relative_permeability_module
  use IAPWS_module
  use fson
  use fson_mpi_module
  use eos_we_module

  implicit none
  private

public :: test_eos_we_fluid_properties, test_eos_we_transition, &
     test_eos_we_errors, test_eos_we_conductivity

contains

!------------------------------------------------------------------------

  subroutine transition_compare(expected_primary, expected_region, &
       expected_transition, primary, fluid, transition, message)

    ! Runs asserts to test EOS transition

    PetscReal, intent(in) :: expected_primary(:), primary(:)
    PetscInt, intent(in) :: expected_region
    type(fluid_type), intent(in) :: fluid
    PetscBool, intent(in) :: expected_transition, transition

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

    call assert_equals(expected_transition, transition, &
         trim(message) // " transition")

  end subroutine transition_compare

!------------------------------------------------------------------------

  subroutine test_eos_we_fluid_properties

    ! eos_we fluid_properties() test

    type(fluid_type) :: fluid
    type(rock_type) :: rock
    PetscInt,  parameter :: offset = 1, region = 4, phase_composition = b'011'
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal, allocatable:: primary(:), primary2(:)
    type(eos_we_type) :: eos
    type(IAPWS_type) :: thermo
    class(relative_permeability_type), allocatable, target :: rp
    type(fson_value), pointer :: json
    character(120) :: json_str = &
         '{"rock": {"relative permeability": {"type": "linear", "liquid": [0.2, 0.8], "vapour": [0.2, 0.8]}}}'
    PetscErrorCode :: err
    PetscReal, parameter :: temperature = 230._dp
    PetscReal, parameter :: pressure = 27.967924557686445e5_dp
    PetscReal, parameter :: vapour_saturation = 0.25
    PetscReal, parameter :: tol = 1.e-8_dp
    PetscReal, parameter :: expected_liquid_density = 827.12247049977032_dp
    PetscReal, parameter :: expected_liquid_internal_energy = 986828.18916209263_dp
    PetscReal, parameter :: expected_liquid_specific_enthalpy = 990209.54144729744_dp
    PetscReal, parameter :: expected_liquid_viscosity = 1.1619412513757267e-4_dp
    PetscReal, parameter :: expected_liquid_relative_permeability = 11._dp / 12._dp
    PetscReal, parameter :: expected_vapour_density = 13.984012253728331_dp
    PetscReal, parameter :: expected_vapour_internal_energy = 2603010.010356456_dp
    PetscReal, parameter :: expected_vapour_specific_enthalpy = 2803009.2956133024_dp
    PetscReal, parameter :: expected_vapour_viscosity = 1.6704837258831552e-5_dp
    PetscReal, parameter :: expected_vapour_relative_permeability = 1._dp / 12._dp

    json => fson_parse_mpi(str = json_str)
    call thermo%init()
    call eos%init(json, thermo)
    call setup_relative_permeabilities(json, rp)

    call fluid%init(eos%num_components, eos%num_phases)
    call rock%init()
    allocate(primary(eos%num_primary_variables), primary2(eos%num_primary_variables))
    allocate(fluid_data(fluid%dof))
    fluid_data = 0._dp
    call fluid%assign(fluid_data, offset)

    rock%relative_permeability => rp

    primary = [pressure, vapour_saturation]
    fluid%region = dble(region)
    call eos%bulk_properties(primary, fluid, err)
    call eos%phase_composition(fluid, err)
    call eos%phase_properties(primary, rock, fluid, err)
    call eos%primary_variables(fluid, primary2)

    if (mpi%rank == mpi%output_rank) then

       call assert_equals(pressure, fluid%pressure, tol, "Pressure")
       call assert_equals(temperature, fluid%temperature, tol, "Temperature")
       call assert_equals(phase_composition, nint(fluid%phase_composition), "Phase composition")

       call assert_equals(expected_liquid_density, fluid%phase(1)%density, &
            tol, "Liquid density")
       call assert_equals(expected_liquid_internal_energy, &
            fluid%phase(1)%internal_energy, tol, "Liquid internal energy")
       call assert_equals(expected_liquid_specific_enthalpy, &
            fluid%phase(1)%specific_enthalpy, tol, "Liquid specific enthalpy")
       call assert_equals(expected_liquid_viscosity, fluid%phase(1)%viscosity, &
            tol, "Liquid viscosity")
       call assert_equals(1._dp - vapour_saturation, &
            fluid%phase(1)%saturation, tol, "Liquid saturation")
       call assert_equals(expected_liquid_relative_permeability, &
            fluid%phase(1)%relative_permeability, &
            tol, "Liquid relative permeability")
       call assert_equals(1._dp, fluid%phase(1)%mass_fraction(1), &
            tol, "Liquid mass fraction")

       call assert_equals(expected_vapour_density, fluid%phase(2)%density, &
            tol, "Vapour density")
       call assert_equals(expected_vapour_internal_energy, &
            fluid%phase(2)%internal_energy, tol, "Vapour internal energy")
       call assert_equals(expected_vapour_specific_enthalpy, &
            fluid%phase(2)%specific_enthalpy, tol, "Vapour specific enthalpy")
       call assert_equals(expected_vapour_viscosity, fluid%phase(2)%viscosity, &
            tol, "Vapour viscosity")
       call assert_equals(vapour_saturation, fluid%phase(2)%saturation, &
            tol, "Vapour saturation")
       call assert_equals(expected_vapour_relative_permeability, &
            fluid%phase(2)%relative_permeability, tol, &
            "Vapour relative permeability")
       call assert_equals(1._dp, fluid%phase(2)%mass_fraction(1), &
            tol, "Vapour mass fraction")

       call assert_equals(pressure, primary2(1), tol, "Primary 1")
       call assert_equals(vapour_saturation, primary2(2), tol, "Primary 2")

    end if

    call fluid%destroy()
    call rock%destroy()
    deallocate(primary, primary2, fluid_data)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(rp)

  end subroutine test_eos_we_fluid_properties

!------------------------------------------------------------------------

  subroutine test_eos_we_transition

    ! eos_we_transition() test

    type(fluid_type) :: old_fluid, fluid
    PetscInt,  parameter :: offset = 1
    PetscReal, pointer, contiguous :: old_fluid_data(:), fluid_data(:)
    PetscReal :: primary(2), expected_primary(2), temperature
    PetscInt :: expected_region
    PetscBool :: transition, expected_transition
    type(eos_we_type) :: eos
    type(IAPWS_type) :: thermo
    type(fson_value), pointer :: json
    character(2) :: json_str = '{}'
    character(60) :: title
    PetscErrorCode :: err
    PetscReal, parameter :: small = 1.e-6_dp

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

    if (mpi%rank == mpi%output_rank) then

       title = "Region 1 null transition"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [1.e5_dp, 20._dp]
       expected_transition = PETSC_FALSE
       primary = expected_primary
       call eos%transition(primary, old_fluid, fluid, transition, err)
       call transition_compare(expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, title)

       title = "Region 1 to 4"
       old_fluid%region = dble(1)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [15.546718682698252e5_dp, small]
       expected_transition = PETSC_TRUE
       primary = [15.e5_dp, 200._dp]
       call eos%transition(primary, old_fluid, fluid, transition, err)
       call transition_compare(expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, title)

       title = "Region 2 null transition"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [1.e5_dp, 120._dp]
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(primary, old_fluid, fluid, transition, err)
       call transition_compare(expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, title)

       title = "Region 2 to 4"
       old_fluid%region = dble(2)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [85.e5_dp, 1._dp - small]
       expected_transition = PETSC_TRUE
       primary = [85.01e5_dp, 299.27215502281706_dp]
       call eos%transition(primary, old_fluid, fluid, transition, err)
       call transition_compare(expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, title)

       title = "Region 4 null transition"
       old_fluid%region = dble(4)
       fluid%region = old_fluid%region
       expected_region = 4
       expected_primary = [1.e5_dp, 0.5_dp]
       primary = expected_primary
       expected_transition = PETSC_FALSE
       call eos%transition(primary, old_fluid, fluid, transition, err)
       call transition_compare(expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, title)

       title = "Region 4 to 1"
       temperature = 299.27215502281706_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 1
       expected_primary = [85.000085e5_dp, temperature]
       expected_transition = PETSC_TRUE
       primary = [85.e5_dp, -0.01_dp]
       call eos%transition(primary, old_fluid, fluid, transition, err)
       call transition_compare(expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, title)

       title = "Region 4 to 2"
       temperature = 212.38453531849041_dp
       old_fluid%region = dble(4)
       old_fluid%temperature = temperature
       fluid%region = old_fluid%region
       expected_region = 2
       expected_primary = [19.99998e5_dp, temperature]
       expected_transition = PETSC_TRUE
       primary = [20.e5_dp, 1.02_dp]
       call eos%transition(primary, old_fluid, fluid, transition, err)
       call transition_compare(expected_primary, expected_region, &
            expected_transition, primary, fluid, transition, title)

    end if

    call old_fluid%destroy()
    call fluid%destroy()
    deallocate(old_fluid_data, fluid_data)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_eos_we_transition

!------------------------------------------------------------------------

  subroutine test_eos_we_errors

    ! eos_we error handling

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
    class(relative_permeability_type), allocatable, target :: rp
    type(fson_value), pointer :: json
    character(120) :: json_str = &
         '{"rock": {"relative permeability": {"type": "linear", "liquid": [0.2, 0.8], "vapour": [0.2, 0.8]}}}'
    PetscReal, parameter :: tol = 1.e-8_dp
    PetscInt :: i, p
    PetscErrorCode :: err

    json => fson_parse_mpi(str = json_str)
    call thermo%init()
    call eos%init(json, thermo)
    call setup_relative_permeabilities(json, rp)

    call fluid%init(num_components, eos%num_phases)
    call rock%init()
    allocate(fluid_data(fluid%dof))
    fluid_data = 0._dp
    call fluid%assign(fluid_data, offset)

    rock%relative_permeability => rp

    do i = 1, n
       p = i
       primary = data(i, :)
       fluid%region = dble(region(i))
       call eos%bulk_properties(primary, fluid, err)
       call eos%phase_composition(fluid, err)
       call eos%phase_properties(primary, rock, fluid, err)
       if (mpi%rank == mpi%output_rank) then
          call assert_equals(1, err, "error code")
       end if
    end do

    call fluid%destroy()
    call rock%destroy()
    deallocate(fluid_data)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(rp)

  end subroutine test_eos_we_errors

!------------------------------------------------------------------------

  subroutine test_eos_we_conductivity

    ! Test eos_we heat conductivity

    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    PetscReal, pointer, contiguous :: fluid_data(:), rock_data(:)
    type(fluid_type) :: fluid
    type(rock_type) :: rock
    PetscInt :: offset = 1
    PetscReal :: cond, expected_cond
    PetscReal, parameter :: tol = 1.e-6_dp

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

    if (mpi%rank == mpi%output_rank) then

       fluid_data(7) = 0.0_dp
       expected_cond = 1.0_dp
       cond = eos%conductivity(rock, fluid)
       call assert_equals(expected_cond, cond, tol, "sl = 0.0")

       fluid_data(7) = 0.25_dp
       expected_cond = 1.25_dp
       cond = eos%conductivity(rock, fluid)
       call assert_equals(expected_cond, cond, tol, "sl = 0.25")

       fluid_data(7) = 0.5_dp
       expected_cond = 1.3535534_dp
       cond = eos%conductivity(rock, fluid)
       call assert_equals(expected_cond, cond, tol, "sl = 0.5")

       fluid_data(7) = 0.75_dp
       expected_cond = 1.4330127_dp
       cond = eos%conductivity(rock, fluid)
       call assert_equals(expected_cond, cond, tol, "sl = 0.75")

       fluid_data(7) = 1.0_dp
       expected_cond = 1.5_dp
       cond = eos%conductivity(rock, fluid)
       call assert_equals(expected_cond, cond, tol, "sl = 1.0")

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
