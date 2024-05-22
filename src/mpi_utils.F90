!   Copyright 2016 University of Auckland.

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

module mpi_utils_module
  !! MPI utilities.

#include <petsc/finclude/petsc.h>

  use petsc

  implicit none
  private

  public :: mpi_broadcast_error_flag, mpi_broadcast_logical
  public :: get_mpi_int_gather_array, mpi_comm_root_world_rank
  public :: mpi_comm_send, mpi_int_gatherv, invert_indices

contains

!------------------------------------------------------------------------

  subroutine mpi_broadcast_error_flag(err)
    !! If integer error flag err is greater than zero on any rank,
    !! broadcast value 1 to all ranks.

    PetscErrorCode, intent(in out) :: err
    ! Locals:
    PetscBool :: any_err
    PetscErrorCode :: ierr

    call MPI_allreduce(err > 0, any_err, 1, MPI_LOGICAL, MPI_LOR, &
         PETSC_COMM_WORLD, ierr)
    if (any_err) then
       err = 1
    else
       err = 0
    end if

  end subroutine mpi_broadcast_error_flag

!------------------------------------------------------------------------

  subroutine mpi_broadcast_logical(val)
    !! If logical val is true on any rank, broadcast to all ranks.

    PetscBool, intent(in out) :: val
    ! Locals:
    PetscBool :: any_val
    PetscErrorCode :: ierr

    call MPI_allreduce(val, any_val, 1, MPI_LOGICAL, MPI_LOR, &
         PETSC_COMM_WORLD, ierr)
    val = any_val

  end subroutine mpi_broadcast_logical

!------------------------------------------------------------------------

  function get_mpi_int_gather_array(comm) result(array)
    !! Returns integer array for use in MPI gather call. This is of
    !! size equal to the number of processes on the root rank, and
    !! size 1 on other ranks (needs to be allocated even though it is
    !! not actually used.) If comm is not specified, the world
    !! communicator is assumed.

    MPI_Comm, intent(in), optional :: comm !! Communicator used for gathering
    PetscInt, allocatable :: array(:)
    ! Locals:
    MPI_Comm :: gather_comm
    PetscMPIInt :: rank, num_procs
    PetscInt :: size
    PetscErrorCode :: ierr

    if (present(comm)) then
       gather_comm = comm
    else
       gather_comm = PETSC_COMM_WORLD
    end if
    call MPI_comm_rank(gather_comm, rank, ierr)
    call MPI_comm_size(gather_comm, num_procs, ierr)

    if (rank == 0) then
       size = num_procs
    else ! have to allocate non-zero size, even if not actually used:
       size = 1
    end if
    allocate(array(size))

  end function get_mpi_int_gather_array

!------------------------------------------------------------------------

  PetscMPIInt function mpi_comm_root_world_rank(comm)
    !! Returns rank in world communicator of the root of the specified
    !! communicator comm, and broadcasts it to all world ranks.

    MPI_Comm, intent(in) :: comm
    ! Locals:
    PetscMPIInt :: rank, comm_rank, root_world_rank, max_root_world_rank
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    root_world_rank = -1
    if (comm /= MPI_COMM_NULL) then
       call MPI_comm_rank(comm, comm_rank, ierr)
       if (comm_rank == 0) then
          root_world_rank = rank
       end if
    end if

    call MPI_allreduce(root_world_rank, max_root_world_rank, 1, &
         MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)
    mpi_comm_root_world_rank = max_root_world_rank

  end function mpi_comm_root_world_rank

!------------------------------------------------------------------------

  subroutine mpi_comm_send(comm, data, from_rank, to_rank)
    !! Sends PetscReal data from from_rank to to_rank, using the
    !! specified communicator.  (Note that from_rank and to_rank must
    !! be the same on all processes.)

    MPI_Comm, intent(in) :: comm
    PetscReal, intent(in out) :: data !! PetscReal data to send
    PetscMPIInt, intent(in) :: from_rank !! Rank to send from
    PetscMPIInt, intent(in) :: to_rank !! Rank to send to
    ! Locals:
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_comm_rank(comm, rank, ierr)

    if (from_rank /= to_rank) then
       if (rank == from_rank) then
          call MPI_send(data, 1, MPI_DOUBLE_PRECISION, to_rank, 0, comm, ierr)
       else if (rank == to_rank) then
          call MPI_recv(data, 1, MPI_DOUBLE_PRECISION, from_rank, 0, comm, &
               MPI_STATUS_IGNORE, ierr)
       end if
    end if

  end subroutine mpi_comm_send

!------------------------------------------------------------------------

  function mpi_int_gatherv(comm, rank, v) result(vall)
    !! Gathers integer array v from all processes in the communicator
    !! onto the root rank. A convenience wrapper for MPI_gatherv().

    use utils_module, only: array_cumulative_sum

    MPI_Comm, intent(in) :: comm
    PetscMPIInt, intent(in) :: rank
    PetscInt, intent(in) :: v(:)
    PetscInt, allocatable :: vall(:)
    ! Locals:
    PetscMPIInt :: comm_size
    PetscInt :: local_count, count
    PetscInt, allocatable :: counts(:), displacements(:)
    PetscErrorCode :: ierr

    call mpi_comm_size(comm, comm_size, ierr)
    counts = get_mpi_int_gather_array(comm)
    displacements = get_mpi_int_gather_array(comm)
    local_count = size(v)

    call MPI_gather(local_count, 1, MPI_INTEGER, &
         counts, 1, MPI_INTEGER, 0, comm, ierr)
    if (rank == 0) then
       displacements = [[0], &
            array_cumulative_sum(counts(1: comm_size - 1))]
       count = sum(counts)
    else
       count = 1
    end if

    allocate(vall(count))
    call MPI_gatherv(v, local_count, MPI_INTEGER, &
         vall, counts, displacements, &
         MPI_INTEGER, 0, comm, ierr)

  end function mpi_int_gatherv

!------------------------------------------------------------------------

  IS function invert_indices(indices, name)
    !! Returns index set containing the inverse of the mapping
    !! represented by the values contained in the indices array on all
    !! processes. The resulting index set contains values only on the
    !! root process and is given the specified object name.

    use utils_module, only: array_cumulative_sum

    PetscInt, intent(in) :: indices(:)
    character(*), intent(in) :: name
    ! Locals:
    PetscInt :: i, local_size
    PetscMPIInt :: rank, num_procs, num_all, is_count
    PetscInt, allocatable :: counts(:), displacements(:)
    PetscInt, allocatable :: indices_all(:), global_indices(:)
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)
    local_size = size(indices)

    counts = get_mpi_int_gather_array()
    displacements = get_mpi_int_gather_array()
    call MPI_gather(local_size, 1, MPI_INTEGER, counts, 1, &
         MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       displacements = [[0], &
            array_cumulative_sum(counts(1: num_procs - 1))]
       num_all = sum(counts)
       is_count = num_all
    else
       num_all = 1
       is_count = 0
    end if

    allocate(indices_all(0: num_all - 1), global_indices(0: num_all - 1))
    global_indices = -1
    call MPI_gatherv(indices, local_size, MPI_INTEGER, &
         indices_all, counts, displacements, &
         MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       do i = 0, num_all - 1
          global_indices(indices_all(i)) = i
       end do
    end if
    deallocate(indices_all, counts, displacements)

    call ISCreateGeneral(PETSC_COMM_WORLD, is_count, &
         global_indices, PETSC_COPY_VALUES, invert_indices, ierr)
    CHKERRQ(ierr)
    call PetscObjectSetName(invert_indices, name, ierr)
    deallocate(global_indices)

  end function invert_indices

!------------------------------------------------------------------------

end module mpi_utils_module
