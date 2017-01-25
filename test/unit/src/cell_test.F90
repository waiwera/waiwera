module cell_test

  ! Test for cell module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use fruit
  use cell_module

  implicit none
  private

public :: test_cell_assign_geometry, test_cell_balance

contains
  
!------------------------------------------------------------------------

  subroutine test_cell_assign_geometry

    ! Cell assign_geometry() test

    type(cell_type) :: cell
    PetscInt, parameter :: num_components = 1, num_phases = 2
    PetscReal, parameter :: volume = 1.e3_dp
    PetscReal, parameter :: centroid(3) = [20._dp, 30._dp, 75._dp]
    PetscInt, parameter :: offset = 6
    PetscReal :: offset_padding(offset-1) = 0._dp
    PetscReal, pointer, contiguous :: cell_data(:)
    PetscReal, parameter :: tol = 1.e-6_dp

    if (mpi%rank == mpi%output_rank) then

       allocate(cell_data(offset - 1 + sum(cell_variable_num_components)))
       cell_data = [offset_padding, centroid, volume]

       call cell%init(num_components, num_phases)

       call assert_equals(cell%dof, size(cell_data) - (offset-1), "cell dof")

       call cell%assign_geometry(cell_data, offset)

       call assert_equals(volume, cell%volume, tol, "volume")
       call assert_equals(0._dp, norm2(cell%centroid - centroid), tol, "centroid")

       call cell%destroy()
       deallocate(cell_data)

    end if

  end subroutine test_cell_assign_geometry

!------------------------------------------------------------------------

  subroutine test_cell_balance
    !! Test cell mass and energy balance routine

    use rock_module, only: rock_variable_num_components
    use fluid_module, only: num_fluid_variables, num_phase_variables
    type(cell_type) :: cell
    PetscReal, pointer, contiguous :: rock_data(:), fluid_data(:)
    PetscInt, parameter :: rock_offset = 1, fluid_offset = 1
    PetscInt, parameter :: num_components = 2, num_phases = 2
    PetscInt, parameter :: num_primary = 3
    PetscReal :: bal(num_primary)
    PetscReal, parameter :: expected_bal(num_components + 1) = &
         [52.372_dp, 22.458_dp, 2.8545448e8_dp]
    PetscReal, parameter :: tol = 1.e-6_dp

    if (mpi%rank == mpi%output_rank) then

       allocate(rock_data(sum(rock_variable_num_components)), &
            fluid_data(num_fluid_variables + num_phases * &
            (num_phase_variables + num_components)))
       rock_data = [0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0.1_dp, 2200._dp, 950._dp]
       fluid_data = [2.7e5_dp, 130._dp, 4._dp, 3._dp, &
            935._dp, 0.0_dp, 0.8_dp, 0.0_dp, 0._dp, 5.461e5_dp, 0.7_dp, 0.3_dp, &
            1.5_dp,  0.0_dp, 0.2_dp, 0.0_dp, 0._dp, 2.540e6_dp, 0.4_dp, 0.6_dp]

       call cell%init(num_components, num_phases)

       call cell%rock%assign(rock_data, rock_offset)
       call cell%fluid%assign(fluid_data, fluid_offset)

       bal = cell%balance(num_primary)

       call assert_equals(expected_bal, bal, num_primary, tol, "cell balance")

       call cell%destroy()
       deallocate(rock_data, fluid_data)

    end if

  end subroutine test_cell_balance

!------------------------------------------------------------------------

end module cell_test
