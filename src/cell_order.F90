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

    use minc_module, only: minc_level_label_name

    DM, intent(in out) :: dm !! Mesh DM
    IS, intent(out) :: cell_index !! Output cell index set
    PetscViewer, intent(in out), optional :: viewer !! Viewer for optional output
    ! Locals:
    PetscInt :: local_count, total_count, local_level_count
    PetscInt :: i, m, inatural
    DMLabel :: ghost_label, order_label, minc_level_label
    PetscInt, allocatable :: global_index(:), natural_index(:)
    PetscInt, allocatable :: global_index_all(:), natural_index_all(:)
    PetscInt, allocatable :: index_array_all(:), index_array(:)
    PetscInt, allocatable :: local_counts(:), displacements(:)
    PetscInt, allocatable :: local_level_counts(:)
    IS :: cell_interior_index
    PetscInt :: start_global_index, max_level
    PetscInt :: cmax, fmax, emax, vmax
    PetscInt :: start_cell, end_cell, end_interior_cell
    PetscMPIInt :: rank, num_procs
    PetscErrorCode :: ierr
    PetscBool :: minc

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)
    call DMPlexGetHybridBounds(dm, cmax, fmax, emax, vmax, ierr)
    end_interior_cell = cmax
    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    CHKERRQ(ierr)
    call DMGetLabel(dm, cell_order_label_name, order_label, ierr)
    CHKERRQ(ierr)
    call DMHasLabel(dm, minc_level_label_name, minc, ierr); CHKERRQ(ierr)

    start_global_index = get_local_start_global_index()
    local_count = get_local_count()
    call get_displacements(local_count, local_counts, displacements, &
         total_count)
    allocate(global_index_all(0: total_count - 1), &
         natural_index_all(0: total_count - 1))

    if (minc) then

       call DMGetLabel(dm, minc_level_label_name, minc_level_label, ierr)
       CHKERRQ(ierr)
       max_level = get_max_minc_level()
       call allocate_processor_array(local_level_counts)
       do m = 0, max_level
          call get_minc_global_and_natural_indices(m, global_index, &
               natural_index, local_level_count)
          call MPI_gather(local_level_count, 1, MPI_INTEGER, &
               local_level_counts, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
          call MPI_gatherv(global_index, local_level_count, MPI_INTEGER, &
               global_index_all, local_level_counts, displacements, MPI_INTEGER, &
               0, PETSC_COMM_WORLD, ierr)
          call MPI_gatherv(natural_index, local_level_count, MPI_INTEGER, &
               natural_index_all, local_level_counts, displacements, MPI_INTEGER, &
               0, PETSC_COMM_WORLD, ierr)
          if (m == 0) then
             inatural = get_new_natural_index()
          else
             call assign_minc_natural_indices(inatural)
          end if
          deallocate(global_index, natural_index)
          displacements = displacements + local_level_counts
       end do

    else

       call get_global_and_natural_indices(global_index, natural_index)
       call MPI_gatherv(global_index, local_count, MPI_INTEGER, &
            global_index_all, local_counts, displacements, MPI_INTEGER, &
            0, PETSC_COMM_WORLD, ierr)
       call MPI_gatherv(natural_index, local_count, MPI_INTEGER, &
            natural_index_all, local_counts, displacements, MPI_INTEGER, &
            0, PETSC_COMM_WORLD, ierr)
       deallocate(global_index, natural_index)

    end if

    ! Set up index array on root process, and scatter:
    allocate(index_array_all(0: total_count - 1))
    if (rank == 0) then
       do i = 0, total_count - 1
          index_array_all(natural_index_all(i)) = global_index_all(i)
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
       do i = 0, total_count - 1
          index_array_all(natural_index_all(i)) = i
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

  contains

!........................................................................

    PetscInt function get_local_count() result(count)
      !! Count interior cells on current processor

      ! Locals:
      PetscInt :: c, ghost

      count = 0
      do c = start_cell, end_interior_cell - 1
         call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
         if (ghost < 0) count = count + 1
      end do

    end function get_local_count

!........................................................................

    PetscInt function get_local_start_global_index() result(start_index)

      !! Returns starting global index on current processor.

      ! Locals:
      Vec :: v
      PetscInt :: blocksize

      call DMGetGlobalVector(dm, v, ierr); CHKERRQ(ierr)
      call VecGetOwnershipRange(v, start_index, &
           PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
      call VecGetBlockSize(v, blocksize, ierr); CHKERRQ(ierr)
      call DMRestoreGlobalVector(dm, v, ierr); CHKERRQ(ierr)
      start_index = start_index / blocksize

    end function get_local_start_global_index

!........................................................................

    subroutine allocate_processor_array(array)

      !! Allocates array with size num_procs on root process, and 1 on
      !! all others.

      PetscInt, allocatable, intent(out) :: array(:)
      ! Locals:
      PetscInt :: size

      if (rank == 0) then
         size = num_procs
      else ! have to allocate non-zero size, even if not actually used:
         size = 1
      end if
      allocate(array(size))

    end subroutine allocate_processor_array

!........................................................................

    subroutine displacements_from_counts(counts, displacements)

      !! Returns displacements array from counts array.

      PetscInt, intent(in) :: counts(:)
      PetscInt, intent(out) :: displacements(:)
      ! Locals:
      PetscInt :: i

      associate(n => size(counts))
        displacements(1) = 0
        do i = 2, n
           displacements(i) = displacements(i - 1) + counts(i - 1)
        end do
      end associate

    end subroutine displacements_from_counts

!........................................................................

    subroutine get_displacements(local_count, local_counts, &
         displacements, total_count)

      !! Gathers local counts from all processors, calculates
      !! displacements in gathered array and total count for gathered
      !! array.

      PetscInt, intent(in) :: local_count
      PetscInt, allocatable, intent(out) :: local_counts(:), displacements(:)
      PetscInt, intent(out) :: total_count

      call allocate_processor_array(local_counts)
      call allocate_processor_array(displacements)

      call MPI_gather(local_count, 1, MPI_INTEGER, local_counts, 1, &
         MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      if (rank == 0) then
         total_count = sum(local_counts)
         call displacements_from_counts(local_counts, displacements)
      else
         total_count = 1
      end if
      
    end subroutine get_displacements

!........................................................................

    subroutine get_global_and_natural_indices(global_index, natural_index)

      !! Sets up global and natural index arrays on each process.

      PetscInt, allocatable, intent(out) :: global_index(:), natural_index(:)
      ! Locals:
      PetscInt :: i, c, ghost, i_global, order

      allocate(global_index(local_count), natural_index(local_count))
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

    end subroutine get_global_and_natural_indices

!........................................................................

    subroutine get_minc_global_and_natural_indices(level, global_index, &
         natural_index, level_count)

      !! Sets up global and natural index arrays on each process for
      !! the given MINC level, and returns count of local MINC cells
      !! for that level.

      PetscInt, intent(in) :: level
      PetscInt, allocatable, intent(out) :: global_index(:), natural_index(:)
      PetscInt, intent(out) :: level_count
      ! Locals:
      PetscInt :: i, c, ghost, i_global, order
      IS :: level_IS
      PetscInt, pointer :: level_cells(:)

      call DMGetStratumSize(dm, minc_level_label_name, level, &
           level_count, ierr); CHKERRQ(ierr)
      if (level_count > 0) then
         allocate(global_index(level_count), natural_index(level_count))
         call DMGetStratumIS(dm, minc_level_label_name, &
              level, level_IS, ierr); CHKERRQ(ierr)
         call ISGetIndicesF90(level_IS, level_cells, ierr); CHKERRQ(ierr)
         do i = 1, level_count
            c = level_cells(i)
            call DMLabelGetValue(ghost_label, c, ghost, ierr)
            if (ghost < 0) then
               call DMLabelGetValue(order_label, c, order, ierr); CHKERRQ(ierr)
               i_global = start_global_index + c - start_cell
               global_index(i) = i_global
               natural_index(i) = order
            end if
         end do
      else
         allocate(global_index(0), natural_index(0))
      end if

    end subroutine get_minc_global_and_natural_indices

!........................................................................

    PetscInt function get_max_minc_level() result(max_level)

      !! Returns maximum MINC level on all processors.

      ! Locals:
      PetscInt :: local_max_level, c, level

      local_max_level = 0
      do c = start_cell, end_interior_cell - 1
         call DMLabelGetValue(minc_level_label, c, level, ierr); CHKERRQ(ierr)
         local_max_level = max(level, local_max_level)
      end do
      call MPI_reduce(local_max_level, max_level, 1, MPI_INTEGER, MPI_MAX, &
           0, PETSC_COMM_WORLD, ierr)
      call MPI_bcast(max_level, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)

    end function get_max_minc_level

!........................................................................

    PetscInt function get_new_natural_index() result(i)

      !! Initialises new natural index for MINC cells.

      ! Locals:
      PetscInt :: local_max_natural, max_natural

      local_max_natural = maxval(natural_index)
      call MPI_reduce(local_max_natural, max_natural, 1, MPI_INTEGER, &
           MPI_MAX, 0, PETSC_COMM_WORLD, ierr)
      if (rank == 0) then
         i = max_natural + 1
      else
         i = 0
      end if

    end function get_new_natural_index

!........................................................................

    subroutine assign_minc_natural_indices(inatural)

      !! Assigns natural indices to MINC cells for a particular level.

      PetscInt, intent(in out) :: inatural
      ! Locals:
      PetscInt :: total_level_count, p, dl, i
      PetscInt, allocatable :: level_natural_index(:), isort(:)
      PetscInt :: level_displacements(num_procs)

      if (rank == 0) then

         call displacements_from_counts(local_level_counts, &
              level_displacements)
         total_level_count = sum(local_level_counts)
         allocate(level_natural_index(total_level_count))

         ! Assemble array with natural indices (currently from
         ! associated fracture cells) for this MINC level from all
         ! processors:
         dl = 0
         do p = 1, num_procs
            associate(d => level_displacements(p), &
                 n => local_level_counts(p))
              level_natural_index(dl: dl + n) = natural_index_all(d: d + n)
              dl = dl + n
            end associate
         end do

         allocate(isort(total_level_count))
         isort = [(i - 1, i = 1, total_level_count)]
         call PetscSortIntWithPermutation(total_level_count, &
              level_natural_index, isort, ierr); CHKERRQ(ierr)
         isort = isort + 1

         ! Assign new natural indices to MINC cells:
         do i = 1, total_level_count
            level_natural_index(isort(i)) = inatural
            inatural = inatural + 1
         end do

         ! Copy new natural indices back:
         dl = 0
         do p = 1, num_procs
            associate(d => level_displacements(p), &
                 n => local_level_counts(p))
              natural_index_all(d: d + n) = level_natural_index(dl: dl + n)
              dl = dl + n
            end associate
         end do

         deallocate(level_natural_index, isort)

      end if

    end subroutine assign_minc_natural_indices

  end subroutine dm_setup_cell_index

!------------------------------------------------------------------------

end module cell_order_module
