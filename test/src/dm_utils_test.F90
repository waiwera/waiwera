module dm_utils_test

  ! Tests for dm_utils module

  use kinds_module
  use mpi_module
  use fruit
  use dm_utils_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

public :: test_dm_cell_normal_face

contains

!------------------------------------------------------------------------

  subroutine test_dm_cell_normal_face

    ! dm_cell_normal_face() test

    DM :: dm
    character(len = 80) :: filename
    PetscErrorCode :: ierr
    PetscInt :: f, dim

    filename = "data/mesh/block3.exo"

    call DMPlexCreateFromFile(mpi%comm, filename, PETSC_TRUE, dm, ierr)
    CHKERRQ(ierr)
    call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)

    if (mpi%rank == mpi%input_rank) then

       call dm_cell_normal_face(dm, 0, [0._dp, 0._dp, 1._dp], f)
       call assert_equals(20, f, "top")

       call dm_cell_normal_face(dm, 2, [0._dp, 0._dp, -1._dp], f)
       call assert_equals(30, f, "bottom")

       call dm_cell_normal_face(dm, 1, [1._dp, 0._dp, 0._dp], f)
       call assert_equals(26, f, "middle +x side")

       call dm_cell_normal_face(dm, 1, [-2._dp, 0._dp, 0.1_dp], f)
       call assert_equals(27, f, "middle -x side")

       call dm_cell_normal_face(dm, 0, [0.1_dp, 0.9_dp, 0._dp], f)
       call assert_equals(23, f, "top +y side")

       call dm_cell_normal_face(dm, 2, [2._dp, -8._dp, 0.1_dp], f)
       call assert_equals(34, f, "bottom -y side")

    end if

    call DMDestroy(dm, ierr); CHKERRQ(ierr)

  end subroutine test_dm_cell_normal_face

!------------------------------------------------------------------------

end module dm_utils_test
