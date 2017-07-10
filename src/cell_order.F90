!   Copyright 2017 University of Auckland.

!   This file is part of Waiwera.

!   Waiwera is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.

!   Waiwera is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.

!   You should have received a copy of the GNU Lesser General Public License
!   along with Waiwera.  If not, see <http://www.gnu.org/licenses/>.

module cell_order_module
  !! Module for setting up cell index on a DM.

#include <petsc/finclude/petsc.h>

  use petsc

  implicit none

  private

  character(len = 13), parameter, public :: &
       cell_order_label_name = "cell_order" !! Name of DMLabel for generating cell index set

  public :: dm_setup_cell_order_label, dm_setup_cell_index
  
contains

!------------------------------------------------------------------------

  subroutine dm_setup_cell_order_label(dm)
    !! Sets up cell ordering label on mesh DM cells. This assumes the DM
    !! has not yet been distributed.

    DM, intent(in out) :: dm
    ! Locals:
    PetscBool :: has_label
    PetscInt :: start_cell, end_cell, c
    PetscErrorCode :: ierr

    call DMHasLabel(dm, cell_order_label_name, has_label, &
         ierr); CHKERRQ(ierr)

    if (.not. has_label) then
       call DMCreateLabel(dm, cell_order_label_name, ierr)
       CHKERRQ(ierr)
       call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, &
            ierr); CHKERRQ(ierr)
       do c = start_cell, end_cell - 1
          call DMSetLabelValue(dm, cell_order_label_name, c, &
               c, ierr); CHKERRQ(ierr)
       end do
    end if

  end subroutine dm_setup_cell_order_label

!------------------------------------------------------------------------

  subroutine dm_setup_cell_index(dm, cell_index, viewer)
    !! Sets up cell index set from cell order label on DM.  This index
    !! set corresponds to a block size of 1.
    !! Also writes the cell interior index set to HDF5 output. This is
    !! used for post-processing and is similar to the cell index set
    !! except that the global indices apply to vectors containing only
    !! interior cells (not boundary ghost cells).

    DM, intent(in out) :: dm !! Mesh DM
    IS, intent(out) :: cell_index !! Output cell index set
    PetscViewer, intent(in out), optional :: viewer !! Viewer for optional output
    ! Locals:
    PetscInt :: total_count, local_count
    PetscInt :: total_allocate_count, allocate_size
    PetscInt :: c, i, ghost, order
    DMLabel :: ghost_label, order_label
    PetscInt, allocatable :: global_index(:), natural_index(:)
    PetscInt, allocatable :: global_index_all(:), natural_index_all(:)
    PetscInt, allocatable :: index_array_all(:), index_array(:)
    PetscInt, allocatable :: local_counts(:), displacements(:)
    Vec :: v
    IS :: cell_interior_index
    PetscInt :: blocksize, start_global_index, i_global
    PetscInt :: cmax, fmax, emax, vmax
    PetscInt :: start_cell, end_cell, end_interior_cell
    PetscMPIInt :: rank, num_procs
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)
    call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call DMPlexGetHybridBounds(dm, cmax, fmax, emax, vmax, ierr)
    CHKERRQ(ierr)
    end_interior_cell = cmax
    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)

    ! Count interior cells:
    local_count = 0
    do c = start_cell, end_interior_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) local_count = local_count + 1
    end do
    allocate(global_index(local_count), natural_index(local_count))

    ! Get starting global index for each process:
    call DMGetGlobalVector(dm, v, ierr); CHKERRQ(ierr)
    call VecGetOwnershipRange(v, start_global_index, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    call VecGetBlockSize(v, blocksize, ierr); CHKERRQ(ierr)
    call DMRestoreGlobalVector(dm, v, ierr); CHKERRQ(ierr)
    start_global_index = start_global_index / blocksize

    ! Set up global and natural index arrays on each process:
    call DMGetLabel(dm, cell_order_label_name, order_label, ierr)
    CHKERRQ(ierr)
    i = 1
    do c = start_cell, end_interior_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          call DMLabelGetValue(order_label, c, order, ierr); CHKERRQ(ierr)
          i_global = start_global_index + c - start_cell
          global_index(i) = i_global
          natural_index(i) = order
          i = i + 1
       end if
    end do

    ! Gather arrays to root process:
    if (rank == 0) then
       allocate_size = num_procs
    else ! have to allocate non-zero size, even if not actually used:
       allocate_size = 1
    end if
    allocate(local_counts(allocate_size), displacements(allocate_size))
    call MPI_gather(local_count, 1, MPI_INTEGER, local_counts, 1, &
         MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       total_count = sum(local_counts)
       displacements(1) = 0
       do i = 2, num_procs
          displacements(i) = displacements(i-1) + local_counts(i-1)
       end do
       total_allocate_count = total_count
    else
       total_allocate_count = 1
    end if
    allocate(global_index_all(total_allocate_count), &
         natural_index_all(total_allocate_count))
    call MPI_gatherv(global_index, local_count, MPI_INTEGER, &
         global_index_all, local_counts, displacements, MPI_INTEGER, &
         0, PETSC_COMM_WORLD, ierr)
    call MPI_gatherv(natural_index, local_count, MPI_INTEGER, &
         natural_index_all, local_counts, displacements, MPI_INTEGER, &
         0, PETSC_COMM_WORLD, ierr)
    deallocate(global_index, natural_index)

    ! Set up index array on root process, and scatter:
    allocate(index_array_all(total_allocate_count))
    if (rank == 0) then
       do i = 1, total_count
          index_array_all(natural_index_all(i) + 1) = global_index_all(i)
       end do
    end if
    deallocate(global_index_all)
    allocate(index_array(local_count))
    call MPI_scatterv(index_array_all, local_counts, displacements, &
         MPI_INTEGER, index_array, local_count, MPI_INTEGER, &
         0, PETSC_COMM_WORLD, ierr)
    call ISCreateGeneral(PETSC_COMM_WORLD, local_count, index_array, &
         PETSC_COPY_VALUES, cell_index, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(cell_index, "cell_index", ierr)

    ! Set up cell interior index set:
    if (rank == 0) then
       do i = 1, total_count
          index_array_all(natural_index_all(i) + 1) = i - 1
       end do
    end if
    deallocate(natural_index_all)
    call MPI_scatterv(index_array_all, local_counts, displacements, &
         MPI_INTEGER, index_array, local_count, MPI_INTEGER, &
         0, PETSC_COMM_WORLD, ierr)
    deallocate(index_array_all, local_counts, displacements)

    call ISCreateGeneral(PETSC_COMM_WORLD, local_count, index_array, &
         PETSC_COPY_VALUES, cell_interior_index, ierr); CHKERRQ(ierr)
    deallocate(index_array)
    call PetscObjectSetName(cell_interior_index, &
         "cell_interior_index", ierr)

    if (present(viewer)) then
       call ISView(cell_index, viewer, ierr); CHKERRQ(ierr)
       call ISView(cell_interior_index, viewer, ierr); CHKERRQ(ierr)
    end if

    call ISDestroy(cell_interior_index, ierr); CHKERRQ(ierr)

  end subroutine dm_setup_cell_index

!------------------------------------------------------------------------  

end module cell_order_module
