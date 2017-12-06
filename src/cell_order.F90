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
  
contains

!------------------------------------------------------------------------

  subroutine dm_setup_cell_index(dm, cell_index, viewer)
    !! Sets up cell index set from cell order label on DM.  This index
    !! set corresponds to a block size of 1, and applies to vectors
    !! which contain data for boundary ghost cells.
    !! Also writes the cell interior index set to HDF5 output. This is
    !! used for post-processing and is similar to the cell index set
    !! except that the global indices apply to vectors containing only
    !! interior cells (not boundary ghost cells).
    !! For MINC meshes, the cell order label is also updated in the
    !! MINC cells, giving each one a unique cell order value.

    use minc_module, only: minc_level_label_name
    use utils_module, only: array_cumulative_sum

    DM, intent(in out) :: dm !! Mesh DM
    IS, intent(out) :: cell_index !! Output cell index set
    PetscViewer, intent(in out), optional :: viewer !! Viewer for optional output
    ! Locals:
    PetscInt :: local_count, total_count, local_level_count, local_minc_count
    PetscInt :: i, m, ip, inatural, iglobal
    DMLabel :: ghost_label, order_label, minc_level_label
    PetscInt, allocatable :: natural_index(:), natural_index_all(:)
    PetscInt, allocatable :: index_array_all(:), index_array(:), local_counts(:)
    PetscInt, allocatable :: bdy_counts(:), displacements(:)
    PetscInt, allocatable :: cumulative_local_counts(:)
    PetscInt, allocatable :: level_displacements(:)
    PetscInt, allocatable :: local_level_counts(:)
    IS :: cell_interior_index
    PetscInt :: max_level
    PetscInt :: cmax, fmax, emax, vmax
    PetscInt :: start_cell, end_cell, end_interior_cell, bdy_count
    PetscMPIInt :: rank, num_procs
    PetscErrorCode :: ierr
    PetscBool :: minc

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)
    call DMPlexGetHybridBounds(dm, cmax, fmax, emax, vmax, ierr)
    end_interior_cell = cmax
    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    bdy_count = max(end_cell - end_interior_cell, 0)
    call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    CHKERRQ(ierr)
    call DMGetLabel(dm, cell_order_label_name, order_label, ierr)
    CHKERRQ(ierr)
    call DMHasLabel(dm, minc_level_label_name, minc, ierr); CHKERRQ(ierr)

    local_count = get_local_count()
    call get_displacements(local_count, local_counts, displacements, &
         bdy_counts, total_count)
    allocate(natural_index_all(0: total_count - 1))

    if (minc) then

       call DMGetLabel(dm, minc_level_label_name, minc_level_label, ierr)
       CHKERRQ(ierr)
       max_level = get_max_minc_level()
       call allocate_process_array(local_level_counts)
       call allocate_process_array(level_displacements)
       level_displacements = displacements
       local_minc_count = 0
       do m = 0, max_level
          call get_minc_natural_indices(m, natural_index, &
               local_level_count)
          call MPI_gather(local_level_count, 1, MPI_INTEGER, &
               local_level_counts, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
          call MPI_gatherv(natural_index, local_level_count, MPI_INTEGER, &
               natural_index_all, local_level_counts, level_displacements, &
               MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
          if (m == 0) then
             inatural = get_new_natural_index()
          else
             call assign_minc_natural_indices(inatural)
          end if
          deallocate(natural_index)
          level_displacements = level_displacements + local_level_counts
          if (m > 0) then
             local_minc_count = local_minc_count + local_level_count
          end if
       end do
       call update_minc_cell_order_label()

    else

       call get_natural_indices(natural_index)
       call MPI_gatherv(natural_index, local_count, MPI_INTEGER, &
            natural_index_all, local_counts, displacements, MPI_INTEGER, &
            0, PETSC_COMM_WORLD, ierr)
       deallocate(natural_index)

    end if

    ! Set up index array on root process, and scatter:
    allocate(index_array_all(0: total_count - 1))
    if (rank == 0) then
       iglobal = 0
       ip = 1
       cumulative_local_counts = array_cumulative_sum(local_counts)
       do i = 0, total_count - 1
          index_array_all(natural_index_all(i)) = iglobal
          iglobal = iglobal + 1
          if (i >= cumulative_local_counts(ip) - 1) then
             ! Add space for boundary ghost cells, and increment
             ! process index ip:
             iglobal = iglobal + bdy_counts(ip)
             ip = ip + 1
          end if
       end do
    end if
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
      !! Count interior cells on current process.

      ! Locals:
      PetscInt :: c, ghost

      count = 0
      do c = start_cell, end_interior_cell - 1
         call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
         if (ghost < 0) count = count + 1
      end do

    end function get_local_count

!........................................................................

    subroutine allocate_process_array(array)

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

    end subroutine allocate_process_array

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
         displacements, bdy_counts, total_count)

      !! Gathers local counts and boundary cell counts from all
      !! processes, calculates displacements in gathered array and
      !! total count for gathered array.

      PetscInt, intent(in) :: local_count
      PetscInt, allocatable, intent(out) :: local_counts(:), &
           bdy_counts(:), displacements(:)
      PetscInt, intent(out) :: total_count

      call allocate_process_array(local_counts)
      call allocate_process_array(bdy_counts)
      call allocate_process_array(displacements)

      call MPI_gather(local_count, 1, MPI_INTEGER, local_counts, 1, &
         MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      call MPI_gather(bdy_count, 1, MPI_INTEGER, bdy_counts, 1, &
         MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      if (rank == 0) then
         total_count = sum(local_counts)
         call displacements_from_counts(local_counts, displacements)
      else
         total_count = 1
      end if
      
    end subroutine get_displacements

!........................................................................

    subroutine get_natural_indices(natural_index)

      !! Sets up natural index array on each process.

      PetscInt, allocatable, intent(out) :: natural_index(:)
      ! Locals:
      PetscInt :: i, c, ghost, order

      allocate(natural_index(0: local_count - 1))
      i = 0
      do c = start_cell, end_interior_cell - 1
         call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
         if (ghost < 0) then
            call DMLabelGetValue(order_label, c, order, ierr); CHKERRQ(ierr)
            natural_index(i) = order
            i = i + 1
         end if
      end do

    end subroutine get_natural_indices

!........................................................................

    subroutine get_minc_natural_indices(level, natural_index, &
         local_level_count)

      !! Sets up natural index array on each process for the given
      !! MINC level, and returns count of local MINC cells for that
      !! level.

      PetscInt, intent(in) :: level
      PetscInt, allocatable, intent(out) :: natural_index(:)
      PetscInt, intent(out) :: local_level_count
      ! Locals:
      PetscInt :: i, c, ghost, order, stratum_size, count
      IS :: level_IS
      PetscInt, pointer :: level_cells(:)

      call DMGetStratumSize(dm, minc_level_label_name, level, &
           stratum_size, ierr); CHKERRQ(ierr)
      if (stratum_size > 0) then
         allocate(natural_index(stratum_size))
         call DMGetStratumIS(dm, minc_level_label_name, &
              level, level_IS, ierr); CHKERRQ(ierr)
         call ISGetIndicesF90(level_IS, level_cells, ierr); CHKERRQ(ierr)
         count = 0
         do i = 1, stratum_size
            c = level_cells(i)
            call DMLabelGetValue(ghost_label, c, ghost, ierr)
            if ((ghost < 0) .and. (c < end_interior_cell)) then
               count = count + 1
               call DMLabelGetValue(order_label, c, order, ierr); CHKERRQ(ierr)
               natural_index(count) = order
            end if
         end do
      else
         count = 0
         allocate(natural_index(0))
      end if

      local_level_count = count
      natural_index = natural_index(1: count)

    end subroutine get_minc_natural_indices

!........................................................................

    PetscInt function get_max_minc_level() result(max_level)

      !! Returns maximum MINC level on all processes.

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
      PetscInt :: total_level_count, p, i
      PetscInt, allocatable :: level_natural_index(:), isort(:)
      PetscInt :: local_level_displacements(num_procs)

      if (rank == 0) then

         call displacements_from_counts(local_level_counts, &
              local_level_displacements)
         total_level_count = sum(local_level_counts)
         allocate(level_natural_index(0: total_level_count - 1))

         ! Assemble array with natural indices (currently from
         ! associated fracture cells) for this MINC level from all
         ! processes:
         do p = 1, num_procs
            associate(dl => local_level_displacements(p), &
                 d => level_displacements(p), &
                 n => local_level_counts(p))
              level_natural_index(dl: dl + n - 1) = &
                   natural_index_all(d: d + n - 1)
            end associate
         end do

         allocate(isort(0: total_level_count - 1))
         isort = [(i, i = 0, total_level_count - 1)]
         call PetscSortIntWithPermutation(total_level_count, &
              level_natural_index, isort, ierr); CHKERRQ(ierr)

         ! Assign new natural indices to MINC cells:
         do i = 0, total_level_count - 1
            level_natural_index(isort(i)) = inatural
            inatural = inatural + 1
         end do

         ! Copy new natural indices back:
         do p = 1, num_procs
            associate(dl => local_level_displacements(p), &
                 d => level_displacements(p), &
                 n => local_level_counts(p))
              natural_index_all(d: d + n - 1) = &
                   level_natural_index(dl: dl + n - 1)
            end associate
         end do

         deallocate(level_natural_index, isort)

      end if

    end subroutine assign_minc_natural_indices

!........................................................................

    subroutine update_minc_cell_order_label()
      !! Updates cell order label in MINC cells on current process.

      ! Locals:
      PetscInt, allocatable :: minc_counts(:), minc_displacements(:)
      PetscInt :: c, i, m, stratum_size, ghost, count, old_order
      IS :: level_IS
      PetscInt, pointer :: level_cells(:)
      DMLabel :: ghost_label, cell_order_label
      PetscErrorCode :: ierr

      call allocate_process_array(minc_counts)
      call MPI_gather(local_minc_count, 1, MPI_INTEGER, minc_counts, &
           1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      call allocate_process_array(minc_displacements)
      minc_displacements = displacements + local_counts - minc_counts

      allocate(natural_index(local_minc_count))

      call MPI_scatterv(natural_index_all, minc_counts, &
           minc_displacements, MPI_INTEGER, natural_index, &
           local_minc_count, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)

      call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
      call DMGetLabel(dm, cell_order_label_name, cell_order_label, &
           ierr); CHKERRQ(ierr)

      count = 0
      do m = 1, max_level
         call DMGetStratumSize(dm, minc_level_label_name, m, &
              stratum_size, ierr); CHKERRQ(ierr)
         if (stratum_size > 0) then
            call DMGetStratumIS(dm, minc_level_label_name, &
                 m, level_IS, ierr); CHKERRQ(ierr)
            call ISGetIndicesF90(level_IS, level_cells, ierr); CHKERRQ(ierr)
            do i = 1, stratum_size
               c = level_cells(i)
               call DMLabelGetValue(ghost_label, c, ghost, ierr)
               if ((ghost < 0) .and. (c < end_interior_cell)) then
                  count = count + 1
                  call DMLabelGetValue(cell_order_label, c, &
                       old_order, ierr); CHKERRQ(ierr)
                  call DMLabelClearValue(cell_order_label, c, &
                       old_order, ierr); CHKERRQ(ierr)
                  call DMLabelSetValue(cell_order_label, c, &
                       natural_index(count), ierr); CHKERRQ(ierr)
               end if
            end do
         end if
      end do

      deallocate(natural_index, minc_counts, minc_displacements)

    end subroutine update_minc_cell_order_label

  end subroutine dm_setup_cell_index

!------------------------------------------------------------------------

end module cell_order_module
