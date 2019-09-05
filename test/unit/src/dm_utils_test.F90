module dm_utils_test

  ! Tests for dm_utils module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use dm_utils_module

  implicit none
  private

  character(len = 512) :: data_path

  public :: setup, teardown
  public :: test_vec_reorder, test_dm_cell_normal_face, test_field_subvector

contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: ios

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()

    call get_environment_variable('WAIWERA_TEST_DATA_PATH', &
         data_path, status = ios)
    if (ios /= 0) data_path = ''

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown()

    PetscErrorCode :: ierr

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

  subroutine test_vec_reorder(test)

    ! vec_reorder() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
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
    call test%assert(expected_local, v_array)

    call VecRestoreArrayF90(v, v_array, ierr); CHKERRQ(ierr)
    call VecDestroy(v, ierr); CHKERRQ(ierr)
    call VecDestroy(v, ierr); CHKERRQ(ierr)
    deallocate(global, expected_local)

  end subroutine test_vec_reorder

!------------------------------------------------------------------------

  subroutine test_dm_cell_normal_face(test)

    ! dm_cell_normal_face() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    DM :: dm
    character(:), allocatable :: filename
    PetscErrorCode :: ierr
    PetscInt :: f
    PetscMPIInt :: rank

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    ! 2D tests:
    filename = trim(adjustl(data_path)) // "mesh/2D.msh"

    call DMPlexCreateFromFile(PETSC_COMM_WORLD, filename, &
         PETSC_TRUE, dm, ierr); CHKERRQ(ierr)

    if (rank == 0) then

       call dm_cell_normal_face(dm, 0, [0._dp, -1._dp, 0._dp], f)
       call test%assert(216, f, "2D cell 0 down")

       call dm_cell_normal_face(dm, 0, [-1._dp, 0._dp, 0._dp], f)
       call test%assert(215, f, "2D cell 0 left")

       call dm_cell_normal_face(dm, 9, [0._dp, -1._dp, 0._dp], f)
       call test%assert(243, f, "2D cell 9 down")

       call dm_cell_normal_face(dm, 59, [1._dp, 0._dp, 0._dp], f)
       call test%assert(348, f, "2D cell 59 right")

       call dm_cell_normal_face(dm, 95, [0._dp, 1._dp, 0._dp], f)
       call test%assert(424, f, "2D cell 95 up")

    end if

    call DMDestroy(dm, ierr); CHKERRQ(ierr)

    ! 3D tests:
    filename = trim(adjustl(data_path)) // "mesh/block3.exo"

    call DMPlexCreateFromFile(PETSC_COMM_WORLD, filename, &
         PETSC_TRUE, dm, ierr); CHKERRQ(ierr)

    if (rank == 0) then

       call dm_cell_normal_face(dm, 0, [0._dp, 0._dp, 1._dp], f)
       call test%assert(20, f, "3D top")

       call dm_cell_normal_face(dm, 2, [0._dp, 0._dp, -1._dp], f)
       call test%assert(30, f, "3D bottom")

       call dm_cell_normal_face(dm, 1, [1._dp, 0._dp, 0._dp], f)
       call test%assert(26, f, "3D middle +x side")

       call dm_cell_normal_face(dm, 1, [-2._dp, 0._dp, 0.1_dp], f)
       call test%assert(27, f, "3D middle -x side")

       call dm_cell_normal_face(dm, 0, [0.1_dp, 0.9_dp, 0._dp], f)
       call test%assert(23, f, "3D top +y side")

       call dm_cell_normal_face(dm, 2, [2._dp, -8._dp, 0.1_dp], f)
       call test%assert(34, f, "3D bottom -y side")

    end if

    call DMDestroy(dm, ierr); CHKERRQ(ierr)

  end subroutine test_dm_cell_normal_face

!------------------------------------------------------------------------

  subroutine test_field_subvector(test)
    ! get_field_subvector test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscErrorCode :: ierr
    PetscInt :: num_vertices, f, i, val
    PetscInt :: local_size, sub_local_size
    DM :: dm
    Vec :: v, sub_v
    PetscReal, pointer :: v_array(:)
    IS :: index_set
    PetscInt, parameter :: num_fields = 3
    PetscInt, parameter :: num_field_components(num_fields) = 1
    PetscInt, parameter :: field_dim(num_fields) = 0
    character(3), parameter :: field_names(num_fields) = ["foo", "bar", "baz"]

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    num_vertices = min(rank, 4)
    call create_path_dm(num_vertices, dm)
    call dm_set_data_layout(dm, num_field_components, field_dim, field_names)

    call DMCreateGlobalVector(dm, v, ierr); CHKERRQ(ierr)
    call VecGetLocalSize(v, local_size, ierr); CHKERRQ(ierr)
    call test%assert(num_vertices * num_fields, local_size, "vec size")

    ! Fill each field with its field index:
    call VecGetArrayF90(v, v_array, ierr); CHKERRQ(ierr)
    val = 0
    do i = 1, size(v_array)
       v_array(i) = dble(val)
       val = val + 1
       if (val >= num_fields) val = 0
    end do
    call VecRestoreArrayF90(v, v_array, ierr); CHKERRQ(ierr)

    do f = 0, num_fields - 1
       call get_field_subvector(v, f, index_set, sub_v)
       call VecGetLocalSize(sub_v, sub_local_size, ierr); CHKERRQ(ierr)
       call test%assert(num_vertices, sub_local_size, "subvec size")
       call VecGetArrayReadF90(sub_v, v_array, ierr); CHKERRQ(ierr)
       call test%assert(all(nint(v_array) == f), "subvec value")
       call VecRestoreArrayReadF90(sub_v, v_array, ierr); CHKERRQ(ierr)
       call VecRestoreSubVector(v, index_set, sub_v, ierr); CHKERRQ(ierr)
       call ISDestroy(index_set, ierr); CHKERRQ(ierr)
    end do

    call VecDestroy(v, ierr); CHKERRQ(ierr)
    call DMDestroy(dm, ierr); CHKERRQ(ierr)

  end subroutine test_field_subvector

!------------------------------------------------------------------------

end module dm_utils_test
