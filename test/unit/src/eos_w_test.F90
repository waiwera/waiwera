module eos_w_test

  ! Tests for eos_w module (isothermal pure water equation of state)

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use fluid_module
  use rock_module
  use relative_permeability_module
  use IAPWS_module
  use fson
  use fson_mpi_module
  use eos_w_module

  implicit none
  private

  public :: setup, teardown, setup_test
  public :: test_eos_w_fluid_properties

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

    test%tolerance = 1.e-8

  end subroutine setup_test

!------------------------------------------------------------------------

  subroutine test_eos_w_fluid_properties(test)

    ! eos_w fluid_properties() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: fluid
    type(rock_type) :: rock
    PetscInt,  parameter :: offset = 1, region = 1, phase_composition = int(b'01')
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal, allocatable :: primary(:), primary2(:)
    type(eos_w_type) :: eos
    type(IAPWS_type) :: thermo
    type(fson_value), pointer :: json
    character(40) :: json_str = '{"eos":{"temperature": 20.0}}'
    PetscErrorCode :: err
    PetscReal, parameter :: pressure = 1.e5_dp
    PetscReal, parameter :: expected_density = 998.20548637769673_dp
    PetscReal, parameter :: expected_internal_energy = 83911.631393167205_dp
    PetscReal, parameter :: expected_specific_enthalpy = 84011.811167136271_dp
    PetscReal, parameter :: expected_viscosity = 1.0015972622270245e-3_dp
    PetscMPIInt :: rank
    PetscInt :: ierr

    json => fson_parse_mpi(str = json_str)
    call thermo%init()
    call eos%init(json, thermo)

    call fluid%init(eos%num_components, eos%num_phases)
    call rock%init()
    allocate(primary(eos%num_components), primary2(eos%num_components), &
         fluid_data(fluid%dof))
    fluid_data = 0._dp
    call fluid%assign(fluid_data, offset)

    primary = pressure
    fluid%region = dble(region)
    call eos%bulk_properties(primary, fluid, err)
    call eos%phase_composition(fluid, err)
    call eos%phase_properties(primary, rock, fluid, err)
    call eos%primary_variables(fluid, primary2)

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call test%assert(pressure, fluid%pressure, "Pressure")
       call test%assert(eos%temperature, fluid%temperature, "Temperature")
       call test%assert(phase_composition, nint(fluid%phase_composition), "Phase composition")

       call test%assert(expected_density, fluid%phase(1)%density, &
            "Density")
       call test%assert(expected_internal_energy, fluid%phase(1)%internal_energy, &
            "Internal energy")
       call test%assert(expected_specific_enthalpy, fluid%phase(1)%specific_enthalpy, &
            "Specific enthalpy")
       call test%assert(expected_viscosity, fluid%phase(1)%viscosity, &
            "Viscosity")

       call test%assert(1._dp, fluid%phase(1)%saturation, "Saturation")
       call test%assert(1._dp, fluid%phase(1)%relative_permeability, "Relative permeability")
       call test%assert(1._dp, fluid%phase(1)%mass_fraction(1), "Mass fraction")

       call test%assert(pressure, primary2(1), "Primary")

    end if

    call fluid%destroy()
    call rock%destroy()
    deallocate(primary, primary2, fluid_data)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_eos_w_fluid_properties

!------------------------------------------------------------------------ 

end module eos_w_test
