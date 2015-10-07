module eos_test

  ! Tests for EOS module

  use kinds_module
  use mpi_module
  use fruit
  use eos_module
  use fluid_module
  use IAPWS_module
  use fson

  implicit none
  private

#include <petsc/finclude/petscdef.h>

public :: test_eos_w_fluid_properties

contains

!------------------------------------------------------------------------

  subroutine test_eos_w_fluid_properties

    ! eos_w fluid_properties() test

    type(fluid_type) :: fluid
    PetscInt, parameter :: num_components = 1, num_phases = 1
    PetscInt,  parameter :: offset = 1, region = 1
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

    if (mpi%rank == mpi%output_rank) then

       json => fson_parse(str = json_str)
       call thermo%init()
       call eos%init(json, thermo)

       call fluid%init(num_components, num_phases)
       allocate(fluid_data(fluid%dof()))
       fluid_data = 0._dp
       call fluid%assign(fluid_data, offset)

       primary = pressure
       fluid%region = dble(region)
       call eos%fluid_properties(primary, fluid)

       call assert_equals(pressure, fluid%pressure, tol, "Pressure")
       call assert_equals(eos%temperature, fluid%temperature, tol, "Temperature")
       call assert_equals(region, nint(fluid%phase_composition), "Phase composition")

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

       call fluid%destroy()
       deallocate(fluid_data)
       call eos%destroy()
       call thermo%destroy()
       call fson_destroy(json)

    end if

  end subroutine test_eos_w_fluid_properties

!------------------------------------------------------------------------

end module eos_test
