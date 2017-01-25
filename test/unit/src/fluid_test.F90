module fluid_test

  ! Test for fluid module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use fruit
  use fluid_module

  implicit none
  private

public :: test_fluid_assign, test_fluid_component_density, &
     test_fluid_energy, test_fluid_enthalpy

contains
  
!------------------------------------------------------------------------

  subroutine test_fluid_assign

    ! Test fluid assign()

    type(fluid_type) :: fluid
    PetscInt, parameter :: num_components = 2, num_phases = 2
    PetscInt,  parameter :: offset = 7
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscInt :: i, ip, nc, phase_dof
    PetscReal, parameter :: tol = 1.e-6_dp

    if (mpi%rank == mpi%output_rank) then

       call fluid%init(num_components, num_phases)

       allocate(fluid_data(offset-1 + fluid%dof))
       do i = 1, size(fluid_data)
          fluid_data(i) = dble(i)
       end do
       call fluid%assign(fluid_data, offset)

       phase_dof = num_phase_variables + num_components

       call assert_equals(fluid_data(offset), fluid%pressure, tol, "pressure")
       call assert_equals(fluid_data(offset+1), fluid%temperature, tol, "temperature")
       call assert_equals(fluid_data(offset+2), fluid%region, tol, "region")
       call assert_equals(fluid_data(offset+3), fluid%phase_composition, tol, "phase composition")

       i = offset + num_fluid_variables
       do ip = 1, num_phases
          call assert_equals(fluid_data(i), &
               fluid%phase(ip)%density, tol, "density")
          call assert_equals(fluid_data(i+1), &
               fluid%phase(ip)%viscosity, tol, "viscosity")
          call assert_equals(fluid_data(i+2), &
               fluid%phase(ip)%saturation, tol, "saturation")
          call assert_equals(fluid_data(i+3), &
               fluid%phase(ip)%relative_permeability, tol, "relative permeability")
          call assert_equals(fluid_data(i+4), &
               fluid%phase(ip)%specific_enthalpy, tol, "specific enthalpy")
          call assert_equals(fluid_data(i+5), &
               fluid%phase(ip)%internal_energy, tol, "internal energy")
          nc = size(fluid%phase(ip)%mass_fraction)
          call assert_equals(0._dp, norm2(fluid_data(i+6: i + 6 + nc-1) - &
               fluid%phase(ip)%mass_fraction), tol, "mass fraction")
          i = i + phase_dof
       end do

       call fluid%destroy()
       deallocate(fluid_data)

    end if

  end subroutine test_fluid_assign

!------------------------------------------------------------------------

  subroutine test_fluid_component_density
    ! Test fluid component_density()

    type(fluid_type) :: fluid
    PetscInt, parameter :: num_components = 2, num_phases = 2
    PetscInt,  parameter :: offset = 1
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal :: cd(num_components)
    PetscReal, parameter :: expected_cd(num_components) = [523.72_dp, 224.58_dp]
    PetscReal, parameter :: tol = 1.e-3_dp

    if (mpi%rank == mpi%output_rank) then

       call fluid%init(num_components, num_phases)
       allocate(fluid_data(offset - 1 + fluid%dof))

       fluid_data = [2.7e5_dp, 130._dp, 4._dp, 3._dp, &
            935._dp, 0.0_dp, 0.8_dp, 0.0_dp, 0._dp, 5.461e5_dp, 0.7_dp, 0.3_dp, &
            1.5_dp,  0.0_dp, 0.2_dp, 0.0_dp, 0._dp, 2.540e6_dp, 0.4_dp, 0.6_dp]

       call fluid%assign(fluid_data, offset)

       cd = fluid%component_density()

       call assert_equals(expected_cd, cd, num_components, tol, "Fluid component density")

       call fluid%destroy()
       deallocate(fluid_data)

    end if

  end subroutine test_fluid_component_density

!------------------------------------------------------------------------

  subroutine test_fluid_energy
    ! Test fluid energy()

    type(fluid_type) :: fluid
    PetscInt, parameter :: num_components = 2, num_phases = 2
    PetscInt,  parameter :: offset = 1
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal :: ef
    PetscReal, parameter :: expected_ef = 4.092448e8_dp
    PetscReal, parameter :: tol = 1.e-3_dp

    if (mpi%rank == mpi%output_rank) then

       call fluid%init(num_components, num_phases)

       allocate(fluid_data(offset - 1 + fluid%dof))
       fluid_data = [2.7e5_dp, 130._dp, 4._dp, 3._dp, &
            935._dp, 0.0_dp, 0.8_dp, 0.0_dp, 0._dp, 5.461e5_dp, 0.7_dp, 0.3_dp, &
            1.5_dp,  0.0_dp, 0.2_dp, 0.0_dp, 0._dp, 2.540e6_dp, 0.4_dp, 0.6_dp]

       call fluid%assign(fluid_data, offset)

       ef = fluid%energy()

       call assert_equals(expected_ef, ef, tol, "Fluid energy")

       call fluid%destroy()
       deallocate(fluid_data)

    end if

  end subroutine test_fluid_energy

!------------------------------------------------------------------------

  subroutine test_fluid_enthalpy
    ! Test fluid phase_mobilities(), phase_flow_fractions(),
    ! component_flow_fractions() and specific_enthalpy()

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
    PetscReal, parameter :: tol = 1.e-6_dp

    if (mpi%rank == mpi%output_rank) then

       call fluid%init(num_components, num_phases)

       allocate(fluid_data(offset - 1 + fluid%dof))
       fluid_data = [2.7e5_dp, 130._dp, 4._dp, 3._dp, &
            935._dp, 1.e-6_dp, 0.8_dp, 0.7_dp, 83.9e3_dp, 5.461e5_dp, 0.7_dp, 0.3_dp, &
            1.5_dp,  2.e-7_dp, 0.2_dp, 0.3_dp, 800.e3_dp, 2.540e6_dp, 0.4_dp, 0.6_dp]

       call fluid%assign(fluid_data, offset)

       mob = fluid%phase_mobilities()
       call assert_equals(expected_mob, mob, num_phases, tol, "Fluid phase mobilities")

       ff = fluid%phase_flow_fractions()
       call assert_equals(expected_ff, ff, num_phases, tol, "Fluid phase flow fractions")

       cff = fluid%component_flow_fractions(ff)

       call assert_equals(expected_cff, cff, num_components, tol, &
            "Fluid component flow fractions")

       h = fluid%specific_enthalpy(ff)
       call assert_equals(expected_h, h, tol, "Fluid enthalpy")

       call fluid%destroy()
       deallocate(fluid_data)

    end if

  end subroutine test_fluid_enthalpy

!------------------------------------------------------------------------

end module fluid_test
