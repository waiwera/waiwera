module rock_test

  ! Test for rock module

  use kinds_module
  use mpi_module
  use fruit
  use rock_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>

public :: test_rock_assign, test_rock_energy

contains

!------------------------------------------------------------------------

  subroutine test_rock_assign

    ! Rock assign() test

    type(rock_type) :: rock
    PetscReal, parameter :: permeability(3) = [1.e-12_dp, 1.e-13_dp, 1.e-14_dp]
    PetscReal, parameter :: porosity = 0.1_dp
    PetscReal, parameter :: wet_conductivity = 2.5_dp
    PetscReal, parameter :: dry_conductivity = 1.5_dp
    PetscReal, parameter :: density = 2200._dp
    PetscReal, parameter :: specific_heat = 1000._dp
    PetscInt,  parameter :: offset = 5
    PetscReal, parameter :: tol = 1.e-6_dp
    PetscReal :: offset_padding(offset-1) = 0._dp
    PetscReal, pointer, contiguous :: rock_data(:)

    if (mpi%rank == mpi%output_rank) then

       call rock%init()
       allocate(rock_data(offset - 1 + rock%dof))
       rock_data = [offset_padding, permeability, wet_conductivity, &
            dry_conductivity, porosity, density, specific_heat]

       call assert_equals(rock%dof, size(rock_data) - (offset-1), "rock dof")

       call rock%assign(rock_data, offset)

       call assert_equals(0._dp, norm2(permeability - rock%permeability), tol, "rock permeability")
       call assert_equals(porosity, rock%porosity, tol, "rock porosity")
       call assert_equals(wet_conductivity, rock%wet_conductivity, tol, "rock wet heat conductivity")
       call assert_equals(dry_conductivity, rock%dry_conductivity, tol, "rock dry heat conductivity")
       call assert_equals(density, rock%density, tol, "rock density")
       call assert_equals(specific_heat, rock%specific_heat, tol, "rock specific heat")
       
       call rock%destroy()
       deallocate(rock_data)

    end if

  end subroutine test_rock_assign

!------------------------------------------------------------------------

  subroutine test_rock_energy

    ! Rock energy() test

    type(rock_type) :: rock
    PetscReal, pointer, contiguous :: rock_data(:)
    PetscReal :: er
    PetscReal, parameter :: temp = 130._dp
    PetscReal, parameter :: expected_er = 2.717e8
    PetscReal, parameter :: tol = 1.e-3_dp
    PetscInt,  parameter :: offset = 1

    if (mpi%rank == mpi%output_rank) then

       call rock%init()
       allocate(rock_data(offset - 1 + rock%dof))
       rock_data = [0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
            0.1_dp, 2200._dp, 950._dp]
       call rock%assign(rock_data, offset)

       er = rock%energy(temp)
       call assert_equals(expected_er, er, tol, "rock energy")

       call rock%destroy()
       deallocate(rock_data)

    end if

  end subroutine test_rock_energy

!------------------------------------------------------------------------

end module rock_test
