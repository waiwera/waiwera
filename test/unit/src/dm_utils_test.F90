module dm_utils_test

  ! Tests for dm_utils module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fruit
  use dm_utils_module

  implicit none
  private

public :: test_vec_reorder, test_dm_cell_normal_face

contains

!------------------------------------------------------------------------

  subroutine test_vec_reorder

    ! vec_reorder() test

    Vec :: v
    IS :: is1, is2
    PetscInt :: istart, iend, num_local, i, j, ij
    PetscInt, allocatable :: global(:), ind1_array(:), ind2_array(:)
    PetscReal, pointer :: v_array(:)
    PetscReal, allocatable :: expected_local(:)
    PetscErrorCode :: ierr
    PetscInt, parameter :: n = 8, bs = 3
    PetscInt, parameter :: ind1(n) = [3, 6, 1, 0, 2, 5, 4, 7]
    PetscInt, parameter :: ind2(n) = [4, 0, 7, 3, 6, 1, 2, 5]
    PetscInt, parameter :: expected(n) = [6, 5, 4, 0, 3, 7, 2, 1]
    PetscReal, parameter :: tol = 1.e-6_dp

    call VecCreate(PETSC_COMM_WORLD, v, ierr); CHKERRQ(ierr)
    call VecSetBlockSize(v, bs, ierr); CHKERRQ(ierr)
    call VecSetType(v, VECMPI, ierr); CHKERRQ(ierr)
    call VecSetSizes(v, PETSC_DECIDE, n * bs, ierr); CHKERRQ(ierr)
    call VecGetOwnershipRange(v, istart, iend, ierr); CHKERRQ(ierr)

    num_local = (iend - istart) / bs
    allocate(global(num_local))
    global(1) = istart / bs + 1
    do i = 1, num_local - 1
       global(i + 1) = global(i) + 1
    end do

    allocate(expected_local(num_local * bs))
    allocate(ind1_array(num_local), ind2_array(num_local))
    call VecGetArrayF90(v, v_array, ierr); CHKERRQ(ierr)
    do i = 1, num_local
       do j = 1, bs
          ij = bs * (i-1) + j
          v_array(ij) = dble(global(i) - 1)
          expected_local(ij) = dble(expected(global(i)))
       end do
       ind1_array(i) = ind1(global(i))
       ind2_array(i) = ind2(global(i))
    end do

    call ISCreateGeneral(PETSC_COMM_WORLD, num_local, ind1_array, &
         PETSC_COPY_VALUES, is1, ierr); CHKERRQ(ierr)
    call ISCreateGeneral(PETSC_COMM_WORLD, num_local, ind2_array, &
         PETSC_COPY_VALUES, is2, ierr); CHKERRQ(ierr)
    deallocate(ind1_array, ind2_array)

    call vec_reorder(v, is1, is2)
    call assert_equals(expected_local, v_array, num_local, tol)

    call VecRestoreArrayF90(v, v_array, ierr); CHKERRQ(ierr)
    call VecDestroy(v, ierr); CHKERRQ(ierr)
    call VecDestroy(v, ierr); CHKERRQ(ierr)
    deallocate(global, expected_local)

  end subroutine test_vec_reorder

!------------------------------------------------------------------------

  subroutine test_dm_cell_normal_face

    ! dm_cell_normal_face() test

    DM :: dm
    character(len = 80) :: filename
    PetscErrorCode :: ierr
    PetscInt :: f
    PetscMPIInt :: rank

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    ! 2D tests:
    filename = "data/mesh/2D.msh"

    call DMPlexCreateFromFile(PETSC_COMM_WORLD, filename, &
         PETSC_TRUE, dm, ierr); CHKERRQ(ierr)

    if (rank == 0) then

       call dm_cell_normal_face(dm, 0, [0._dp, -1._dp, 0._dp], f)
       call assert_equals(216, f, "2D cell 0 down")

       call dm_cell_normal_face(dm, 0, [-1._dp, 0._dp, 0._dp], f)
       call assert_equals(215, f, "2D cell 0 left")

       call dm_cell_normal_face(dm, 9, [0._dp, -1._dp, 0._dp], f)
       call assert_equals(243, f, "2D cell 9 down")

       call dm_cell_normal_face(dm, 59, [1._dp, 0._dp, 0._dp], f)
       call assert_equals(348, f, "2D cell 59 right")

       call dm_cell_normal_face(dm, 95, [0._dp, 1._dp, 0._dp], f)
       call assert_equals(424, f, "2D cell 95 up")

    end if

    call DMDestroy(dm, ierr); CHKERRQ(ierr)

    ! 3D tests:
    filename = "data/mesh/block3.exo"

    call DMPlexCreateFromFile(PETSC_COMM_WORLD, filename, &
         PETSC_TRUE, dm, ierr); CHKERRQ(ierr)

    if (rank == 0) then

       call dm_cell_normal_face(dm, 0, [0._dp, 0._dp, 1._dp], f)
       call assert_equals(20, f, "3D top")

       call dm_cell_normal_face(dm, 2, [0._dp, 0._dp, -1._dp], f)
       call assert_equals(30, f, "3D bottom")

       call dm_cell_normal_face(dm, 1, [1._dp, 0._dp, 0._dp], f)
       call assert_equals(26, f, "3D middle +x side")

       call dm_cell_normal_face(dm, 1, [-2._dp, 0._dp, 0.1_dp], f)
       call assert_equals(27, f, "3D middle -x side")

       call dm_cell_normal_face(dm, 0, [0.1_dp, 0.9_dp, 0._dp], f)
       call assert_equals(23, f, "3D top +y side")

       call dm_cell_normal_face(dm, 2, [2._dp, -8._dp, 0.1_dp], f)
       call assert_equals(34, f, "3S bottom -y side")

    end if

    call DMDestroy(dm, ierr); CHKERRQ(ierr)

  end subroutine test_dm_cell_normal_face

!------------------------------------------------------------------------

end module dm_utils_test
