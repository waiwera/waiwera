module fluid_test

  ! Test for fluid module

  use kinds_module
  use mpi_module
  use fruit
  use fluid_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>

  public :: test_fluid_assign, test_fluid_component_density, &
       test_fluid_component_mass_fraction, &
       test_fluid_energy, test_fluid_energy_production

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

  subroutine test_fluid_component_mass_fraction
    ! Test fluid component_mass_fraction()

    type(fluid_type) :: fluid
    PetscInt, parameter :: num_components = 2, num_phases = 2
    PetscInt,  parameter :: offset = 1
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal :: x
    PetscInt :: c
    character(60) :: cstr
    PetscReal :: expected_x(num_components)
    PetscReal, parameter :: tol = 1.e-6_dp

    if (mpi%rank == mpi%output_rank) then

       call fluid%init(num_components, num_phases)
       allocate(fluid_data(offset - 1 + fluid%dof))

       fluid_data = [2.7e5_dp, 130._dp, 4._dp, 3._dp, &
            935._dp, 0.0_dp, 0.8_dp, 0.0_dp, 0._dp, 5.461e5_dp, 0.7_dp, 0.3_dp, &
            1.5_dp,  0.0_dp, 0.2_dp, 0.0_dp, 0._dp, 2.540e6_dp, 0.4_dp, 0.6_dp]

       call fluid%assign(fluid_data, offset)
       expected_x = [0.699879727382_dp, 0.300120272618_dp]

       do c = 1, num_components
          write(cstr, '(a, i1, a)') "Fluid component mass fraction ", c, " case 1"
          x = fluid%component_mass_fraction(c)
          call assert_equals(expected_x(c), x, tol, trim(cstr))
       end do

       fluid_data = [2.7e5_dp, 130._dp, 4._dp, 3._dp, &
            935._dp, 0.0_dp, 0.8_dp, 0.0_dp, 0._dp, 5.461e5_dp, 1.0_dp, 0.0_dp, &
            1.5_dp,  0.0_dp, 0.2_dp, 0.0_dp, 0._dp, 2.540e6_dp, 1.0_dp, 0.0_dp]

       call fluid%assign(fluid_data, offset)
       expected_x = [1._dp, 0._dp]

       do c = 1, num_components
          write(cstr, '(a, i1, a)') "Fluid component mass fraction ", c, " case 2"
          x = fluid%component_mass_fraction(c)
          call assert_equals(expected_x(c), x, tol, trim(cstr))
       end do

       call fluid%destroy()
       deallocate(fluid_data)

    end if

  end subroutine test_fluid_component_mass_fraction

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

  subroutine test_fluid_energy_production
    ! Test fluid energy_production()

    type(fluid_type) :: fluid
    PetscInt, parameter :: num_components = 1, num_phases = 2
    PetscBool, parameter :: isothermal = .false.
    PetscInt, parameter :: offset = 1
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal :: source(num_components + 1), ff(num_phases)
    PetscReal, parameter :: tol = 1.e-6_dp
    PetscReal, parameter :: expected_production = -3190284.011514185_dp
    PetscReal, parameter :: expected_flow_fractions(num_phases) = &
         [0.98532727106605555_dp, 0.014672728933944498_dp]

    if (mpi%rank == mpi%output_rank) then

       call fluid%init(num_components, num_phases)

       allocate(fluid_data(offset - 1 + fluid%dof))
       fluid_data = [3346651.871510162_dp, 240._dp, 4._dp, 3._dp, &
            813.36485916981576_dp, 0.00011105570007981882_dp, 0.8_dp, &
            0.9_dp, 1037522.7548256445_dp, 1033408.1784042757_dp, 1._dp, &
            16.747578872158215_dp, 1.7062182385337129e-05_dp, 0.2_dp, &
            0.1_dp, 2803059.9721381501_dp, 2603230.9761348078_dp, 1._dp]

       call fluid%assign(fluid_data, offset)

       ff = fluid%flow_fractions()
       call assert_equals(expected_flow_fractions(1), ff(1), tol, "Flow fraction 1")
       call assert_equals(expected_flow_fractions(2), ff(2), tol, "Flow fraction 2")

       source = [-3._dp, 0._dp]
       call fluid%energy_production(source, isothermal)

       call assert_equals(expected_production, source(num_components + 1), &
            tol, "Energy production")

       call fluid%destroy()
       deallocate(fluid_data)

    end if

  end subroutine test_fluid_energy_production

!------------------------------------------------------------------------

end module fluid_test
