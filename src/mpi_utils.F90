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

#include <petsc/finclude/petscsys.h>

  use petscsys

  implicit none
  private

  public :: mpi_broadcast_error_flag, mpi_broadcast_logical
  public :: get_mpi_int_gather_array, mpi_comm_root_world_rank
  public :: mpi_comm_send

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

end module mpi_utils_module
