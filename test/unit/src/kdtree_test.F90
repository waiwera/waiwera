module kdtree_test

  ! Tests for k-d tree module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use kdtree_module

  implicit none
  private

  character(len = 512) :: data_path

  public :: setup, teardown
  public :: test_kdtree_linear, test_kdtree_mesh

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

  subroutine test_kdtree_linear(test)

    !! Test k-d tree on linear mesh

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call linear_test_case(2)
       call linear_test_case(3)

    end if

  contains

    subroutine linear_test_case(dim)

      !! test case for dimension dim

      PetscInt, intent(in) :: dim
      ! Locals:
      type(kdtree_type) :: kdt
      PetscReal, allocatable :: data(:,:), x(:)
      PetscInt :: n, i
      PetscInt, allocatable :: index(:), expected_index(:)
      PetscErrorCode, allocatable :: err(:)
      character(12) :: msg

      write(msg, '(a, i0)') 'linear dim ', dim
      n = 1000
      allocate(data(dim, n), index(n), expected_index(n), err(n))
      do i = 1, n
         data(1, i) = i - 0.5_dp
         data(2:, i) = 0.5_dp
      end do

      call kdt%init(data)

      do i = 1, n
         x = data(:, i)
         call kdt%search(x, index(i), err = err(i))
         expected_index(i) = i
      end do

      call test%assert(index, expected_index, msg // ' indices')
      call test%assert(all(err == 0), msg // ' error')

      call kdt%destroy()
      deallocate(data, x, index, expected_index, err)

    end subroutine linear_test_case

  end subroutine test_kdtree_linear

!------------------------------------------------------------------------

  subroutine test_kdtree_mesh(test)

    !! Test k-d tree on mesh centroids

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    character(24), parameter :: filenames(5) = &
         ['   col100.exo', '  7x7grid.exo', '       2D.msh', &
         '   triopt.msh', '3Drefined.msh']
    PetscInt :: i

    do i = 1, size(filenames)
       call mesh_test_case(trim(adjustl(data_path)) // 'mesh/' // &
            trim(adjustl(filenames(i))))
    end do

  contains

    subroutine mesh_test_case(filename)

      !! Test case for mesh centroids

      use dm_utils_module, only: local_vec_section, section_offset, &
           dm_set_fv_adjacency, dm_set_default_data_layout
      use cell_module, only: cell_type

      character(*), intent(in) :: filename
      ! Locals:
      DM :: dm, dist_dm
      PetscSF :: dist_sf
      PetscErrorCode :: ierr
      Vec :: cell_geom, face_geom
      PetscReal, allocatable :: centroids(:, :), x(:)
      PetscInt :: c, start_cell, end_cell, offset, i, n, dim, np
      PetscSection :: cell_section
      type(cell_type) :: cell
      PetscReal, contiguous, pointer :: cell_geom_array(:)
      type(kdtree_type) :: kdt
      PetscInt, allocatable :: index(:), expected_index(:)
      PetscErrorCode, allocatable :: err(:)
      character(:), allocatable :: msg

      call MPI_comm_size(PETSC_COMM_WORLD, np, ierr)

      call DMPlexCreateFromFile(PETSC_COMM_WORLD, filename, PETSC_TRUE, &
           dm, ierr); CHKERRQ(ierr)
      call DMSetBasicAdjacency(dm, PETSC_TRUE, PETSC_FALSE, ierr)
      if (np > 1) then
         call DMPlexDistribute(dm, 1, dist_sf, dist_dm, ierr)
         if (dist_dm .ne. PETSC_NULL_DM) then
            call DMDestroy(dm, ierr)
            dm = dist_dm
            call PetscSFDestroy(dist_sf, ierr)
         end if
      end if

      call DMPlexComputeGeometryFVM(dm, cell_geom, face_geom, ierr); CHKERRQ(ierr)

      call local_vec_section(cell_geom, cell_section)
      call VecGetArrayReadF90(cell_geom, cell_geom_array, ierr); CHKERRQ(ierr)

      call cell%init(1, 1)
      call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
      CHKERRQ(ierr)
      n = end_cell - start_cell
      call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)
      allocate(centroids(dim, n))
      do c = start_cell, end_cell - 1
         offset = section_offset(cell_section, c)
         call cell%assign_geometry(cell_geom_array, offset)
         i = c - start_cell + 1
         centroids(:, i) = cell%centroid(1: dim)
      end do

      call cell%destroy()
      call VecRestoreArrayReadF90(cell_geom, cell_geom_array, ierr); CHKERRQ(ierr)

      call VecDestroy(cell_geom, ierr); CHKERRQ(ierr)
      call VecDestroy(face_geom, ierr); CHKERRQ(ierr)
      call DMDestroy(dm, ierr); CHKERRQ(ierr)

      call kdt%init(centroids)
      allocate(x(dim), index(n), expected_index(n), err(n))
      msg = 'mesh ' // filename

      do i = 1, n
         x = centroids(:, i)
         call kdt%search(x, index(i), err = err(i))
         expected_index(i) = i
      end do

      call test%assert(index, expected_index, msg // ' indices')
      call test%assert(all(err == 0), msg // ' error')

      call kdt%destroy()
      deallocate(centroids, x, index, expected_index, err)

    end subroutine mesh_test_case

  end subroutine test_kdtree_mesh

!------------------------------------------------------------------------

end module kdtree_test
