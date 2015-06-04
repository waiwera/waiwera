module fluid_test

  ! Test for fluid module

  use kinds_module
  use mpi_module
  use fruit
  use fluid_module

  implicit none
  private

#include <petsc-finclude/petscdef.h>

public :: test_fluid_assign, test_fluid_component_density

PetscReal, parameter :: tol = 1.e-6_dp

contains
  
!------------------------------------------------------------------------

  subroutine test_fluid_assign

    ! Test fluid assign()

    type(fluid_type) :: fluid
    PetscInt, parameter :: num_components = 2, num_phases = 2
    PetscInt,  parameter :: offset = 7
    PetscReal, allocatable :: fluid_data(:)
    PetscInt :: i, p, nc

    if (mpi%rank == mpi%output_rank) then

       call fluid%init(num_components, num_phases)

       allocate(fluid_data(offset-1 + fluid%dof()))
       do i = 1, size(fluid_data)
          fluid_data(i) = dble(i)
       end do
       call fluid%assign(fluid_data, offset)

       call assert_equals(fluid_data(offset), fluid%pressure, tol, "pressure")
       call assert_equals(fluid_data(offset+1), fluid%temperature, tol, "temperature")
       call assert_equals(fluid_data(offset+2), fluid%region, tol, "region")

       i = offset + num_fluid_variables
       do p = 1, num_phases
          call assert_equals(fluid_data(i), &
               fluid%phase(p)%density, tol, "density")
          call assert_equals(fluid_data(i+1), &
               fluid%phase(p)%viscosity, tol, "viscosity")
          call assert_equals(fluid_data(i+2), &
               fluid%phase(p)%saturation, tol, "saturation")
          call assert_equals(fluid_data(i+3), &
               fluid%phase(p)%relative_permeability, tol, "relative permeability")
          call assert_equals(fluid_data(i+4), &
               fluid%phase(p)%specific_enthalpy, tol, "specific enthalpy")
          call assert_equals(fluid_data(i+5), &
               fluid%phase(p)%internal_energy, tol, "internal energy")
          nc = size(fluid%phase(p)%mass_fraction)
          call assert_equals(0._dp, norm2(fluid_data(i+6: i + 6 + nc-1) - &
               fluid%phase(p)%mass_fraction), tol, "mass fraction")
          i = i + fluid%phase(p)%dof()
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
    PetscReal, allocatable :: fluid_data(:)
    PetscReal :: cd(num_components)
    PetscReal, parameter :: expected_cd(num_components) = [523.72_dp, 224.58_dp]
    PetscReal, parameter :: tol = 1.e-3_dp

    if (mpi%rank == mpi%output_rank) then

       call fluid%init(num_components, num_phases)

       fluid_data = [2.7e5_dp, 130._dp, 4._dp, &
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

end module fluid_test
