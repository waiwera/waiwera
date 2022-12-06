module fluid_test

  ! Test for fluid module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use zofu
  use fluid_module

  implicit none
  private

  public :: setup, teardown
  public :: test_fluid_assign, test_fluid_component_density, &
     test_fluid_energy, test_fluid_enthalpy

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

  subroutine test_fluid_assign(test)

    ! Test fluid assign()

    class(unit_test_type), intent(in out) :: test
    type(fluid_type) :: fluid
    PetscInt, parameter :: num_components = 2, num_phases = 2
    PetscInt,  parameter :: offset = 7
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscInt :: i, ip, nc, phase_dof
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscInt, parameter :: expected_dof = 7 + 2 * 9

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call fluid%init(num_components, num_phases)
       call test%assert(expected_dof, fluid%dof, "dof")

       allocate(fluid_data(offset-1 + fluid%dof))
       do i = 1, size(fluid_data)
          fluid_data(i) = dble(i) ! arbitrary unique values
       end do
       call fluid%assign(fluid_data, offset)

       phase_dof = num_phase_variables + num_components - 1

       call test%assert(fluid_data(offset), fluid%pressure, "pressure")
       call test%assert(fluid_data(offset+1), fluid%temperature, "temperature")
       call test%assert(fluid_data(offset+2), fluid%region, "region")
       call test%assert(fluid_data(offset+3), fluid%phase_composition, "phase composition")
       call test%assert(fluid_data(offset+4), fluid%permeability_factor, "permeability factor")
       call test%assert(fluid_data(offset+5: offset+5 + num_components - 1), &
            fluid%partial_pressure, "partial pressure")

       i = offset + num_fluid_variables + num_components - 1
       do ip = 1, num_phases
          call test%assert(fluid_data(i), &
               fluid%phase(ip)%density, "density")
          call test%assert(fluid_data(i+1), &
               fluid%phase(ip)%viscosity, "viscosity")
          call test%assert(fluid_data(i+2), &
               fluid%phase(ip)%saturation, "saturation")
          call test%assert(fluid_data(i+3), &
               fluid%phase(ip)%relative_permeability, "relative permeability")
          call test%assert(fluid_data(i+4), &
               fluid%phase(ip)%capillary_pressure, "capillary pressure")
          call test%assert(fluid_data(i+5), &
               fluid%phase(ip)%specific_enthalpy, "specific enthalpy")
          call test%assert(fluid_data(i+6), &
               fluid%phase(ip)%internal_energy, "internal energy")
          nc = size(fluid%phase(ip)%mass_fraction)
          call test%assert(0._dp, norm2(fluid_data(i+7: i + 7 + nc-1) - &
               fluid%phase(ip)%mass_fraction), "mass fraction")
          i = i + phase_dof
       end do

       call fluid%destroy()
       deallocate(fluid_data)

    end if

  end subroutine test_fluid_assign

!------------------------------------------------------------------------

  subroutine test_fluid_component_density(test)
    ! Test fluid component_density()

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: fluid
    PetscInt, parameter :: num_components = 2, num_phases = 2
    PetscInt,  parameter :: offset = 1
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal :: cd(num_components)
    PetscReal, parameter :: expected_cd(num_components) = [523.72_dp, 224.58_dp]
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       call fluid%init(num_components, num_phases)
       allocate(fluid_data(offset - 1 + fluid%dof))

       fluid_data = [2.7e5_dp, 130._dp, 4._dp, 3._dp, 1._dp, 0._dp, 0._dp, &
            935._dp, 0._dp, 0.8_dp, 0._dp, 0._dp, 0._dp, 5.461e5_dp, 0.7_dp, 0.3_dp, &
            1.5_dp,  0._dp, 0.2_dp, 0._dp, 0._dp, 0._dp, 2.540e6_dp, 0.4_dp, 0.6_dp]

       call fluid%assign(fluid_data, offset)

       cd = fluid%component_density()

       call test%assert(expected_cd, cd, "Fluid component density")

       call fluid%destroy()
       deallocate(fluid_data)

    end if

  end subroutine test_fluid_component_density

!------------------------------------------------------------------------

  subroutine test_fluid_energy(test)
    ! Test fluid energy()

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: fluid
    PetscInt, parameter :: num_components = 2, num_phases = 2
    PetscInt,  parameter :: offset = 1
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal :: ef
    PetscReal, parameter :: expected_ef = 4.092448e8_dp
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call fluid%init(num_components, num_phases)

       allocate(fluid_data(offset - 1 + fluid%dof))
       fluid_data = [2.7e5_dp, 130._dp, 4._dp, 3._dp, 1._dp, 0._dp, 0._dp, &
            935._dp, 0._dp, 0.8_dp, 0._dp, 0._dp, 0._dp, 5.461e5_dp, 0.7_dp, 0.3_dp, &
            1.5_dp,  0._dp, 0.2_dp, 0._dp, 0._dp, 0._dp, 2.540e6_dp, 0.4_dp, 0.6_dp]

       call fluid%assign(fluid_data, offset)

       ef = fluid%energy()

       call test%assert(expected_ef, ef, "Fluid energy")

       call fluid%destroy()
       deallocate(fluid_data)

    end if

  end subroutine test_fluid_energy

!------------------------------------------------------------------------

  subroutine test_fluid_enthalpy(test)
    ! Test fluid phase_mobilities(), phase_flow_fractions(),
    ! component_flow_fractions() and specific_enthalpy()

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: fluid
    PetscInt, parameter :: num_components = 2, num_phases = 2
    PetscInt,  parameter :: offset = 1
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal :: mob(num_phases), ff(num_phases), cff(num_components)
    PetscReal :: h
    PetscReal, parameter :: expected_mob(num_phases) = [654500000._dp, 2250000._dp]
    PetscReal, parameter :: expected_ff(num_phases) = [0.9965740388_dp, 0.0034259612_dp]
    PetscReal, parameter :: expected_cff(num_phases) = [0.6989722116_dp, 0.3010277884_dp]
    PetscReal, parameter :: expected_h = 86353.3307955843_dp
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call fluid%init(num_components, num_phases)

       allocate(fluid_data(offset - 1 + fluid%dof))
       fluid_data = [2.7e5_dp, 130._dp, 4._dp, 3._dp, 1._dp, 0._dp, 0._dp, &
            935._dp, 1.e-6_dp, 0.8_dp, 0.7_dp, 0._dp, 83.9e3_dp, 5.461e5_dp, 0.7_dp, 0.3_dp, &
            1.5_dp,  2.e-7_dp, 0.2_dp, 0.3_dp, 0._dp, 800.e3_dp, 2.540e6_dp, 0.4_dp, 0.6_dp]

       call fluid%assign(fluid_data, offset)

       mob = fluid%phase_mobilities()
       call test%assert(expected_mob, mob, "Fluid phase mobilities")

       ff = fluid%phase_flow_fractions()
       call test%assert(expected_ff, ff, "Fluid phase flow fractions")

       cff = fluid%component_flow_fractions(ff)

       call test%assert(expected_cff, cff, "Fluid component flow fractions")

       h = fluid%specific_enthalpy(ff)
       call test%assert(expected_h, h, "Fluid enthalpy")

       call fluid%destroy()
       deallocate(fluid_data)

    end if

  end subroutine test_fluid_enthalpy

!------------------------------------------------------------------------

end module fluid_test
