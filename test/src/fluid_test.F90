module fluid_test

  ! Test for fluid module

  use kinds_module
  use mpi_module
  use fruit
  use fluid_module

  implicit none
  private

#include <petsc-finclude/petscdef.h>

public :: test_fluid

PetscReal, parameter :: tol = 1.e-6_dp

contains
  
!------------------------------------------------------------------------

  subroutine test_fluid

    ! Fluid test

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
          nc = size(fluid%phase(p)%mass_fraction)
          call assert_equals(0._dp, norm2(fluid_data(i+5: i + 5 + nc-1) - &
               fluid%phase(p)%mass_fraction), tol, "mass fraction")
          i = i + fluid%phase(p)%dof()
       end do

       deallocate(fluid_data)

    end if

  end subroutine test_fluid

!------------------------------------------------------------------------

end module fluid_test
