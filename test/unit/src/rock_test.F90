module rock_test

  ! Test for rock module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use zofu
  use rock_module

  implicit none
  private

  public :: setup, teardown, setup_test
  public :: test_rock_assign, test_rock_energy

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

    test%tolerance = 1.e-9

  end subroutine setup_test

!------------------------------------------------------------------------

  subroutine test_rock_assign(test)

    ! Rock assign() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(rock_type) :: rock
    PetscReal, parameter :: permeability(3) = [1.e-12_dp, 1.e-13_dp, 1.e-14_dp]
    PetscReal, parameter :: porosity = 0.1_dp
    PetscReal, parameter :: wet_conductivity = 2.5_dp
    PetscReal, parameter :: dry_conductivity = 1.5_dp
    PetscReal, parameter :: density = 2200._dp
    PetscReal, parameter :: specific_heat = 1000._dp
    PetscInt,  parameter :: offset = 5
    PetscReal :: offset_padding(offset-1) = 0._dp
    PetscReal, pointer, contiguous :: rock_data(:)
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call rock%init()
       allocate(rock_data(offset - 1 + rock%dof))
       rock_data = [offset_padding, permeability, wet_conductivity, &
            dry_conductivity, porosity, density, specific_heat]

       call test%assert(rock%dof, size(rock_data) - (offset-1), "rock dof")

       call rock%assign(rock_data, offset)

       call test%assert(permeability, rock%permeability, "rock permeability")
       call test%assert(porosity, rock%porosity, "rock porosity")
       call test%assert(wet_conductivity, rock%wet_conductivity, "rock wet heat conductivity")
       call test%assert(dry_conductivity, rock%dry_conductivity, "rock dry heat conductivity")
       call test%assert(density, rock%density, "rock density")
       call test%assert(specific_heat, rock%specific_heat, "rock specific heat")
       
       call rock%destroy()
       deallocate(rock_data)

    end if

  end subroutine test_rock_assign

!------------------------------------------------------------------------

  subroutine test_rock_energy(test)

    ! Rock energy() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(rock_type) :: rock
    PetscReal, pointer, contiguous :: rock_data(:)
    PetscReal :: er
    PetscReal, parameter :: temp = 130._dp
    PetscReal, parameter :: expected_er = 2.717e8
    PetscInt,  parameter :: offset = 1
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call rock%init()
       allocate(rock_data(offset - 1 + rock%dof))
       rock_data = [0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
            0.1_dp, 2200._dp, 950._dp]
       call rock%assign(rock_data, offset)

       er = rock%energy(temp)
       call test%assert(expected_er, er, "rock energy")

       call rock%destroy()
       deallocate(rock_data)

    end if

  end subroutine test_rock_energy

!------------------------------------------------------------------------

end module rock_test
