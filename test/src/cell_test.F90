module cell_test

  ! Test for cell module

  use kinds_module
  use mpi_module
  use fruit
  use cell_module

  implicit none
  private

#include <petsc-finclude/petscdef.h>

public :: test_cell

PetscReal, parameter :: tol = 1.e-6_dp

contains
  
!------------------------------------------------------------------------

  subroutine test_cell

    ! Cell test

    type(cell_type) :: cell
    PetscInt, parameter :: num_components = 1, num_phases = 2
    PetscReal, parameter :: volume = 1.e3_dp
    PetscReal, parameter :: centroid(3) = [20._dp, 30._dp, 75._dp]
    PetscInt, parameter :: offset = 6
    PetscReal :: offset_padding(offset-1) = 0._dp
    PetscReal, allocatable :: cell_data(:)

    if (mpi%rank == mpi%output_rank) then

       cell_data = [offset_padding, centroid, volume]

       call cell%init(num_components, num_phases)

       call assert_equals(cell%dof(), size(cell_data) - (offset-1), "cell dof")

       call cell%assign(cell_data, offset)

       call assert_equals(volume, cell%volume, tol, "volume")
       call assert_equals(0._dp, norm2(cell%centroid - centroid), tol, "centroid")

       call cell%destroy()
       deallocate(cell_data)

    end if

  end subroutine test_cell

!------------------------------------------------------------------------

end module cell_test
