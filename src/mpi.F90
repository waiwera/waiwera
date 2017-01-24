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

module mpi_module
  !! MPI variables.

#include <petsc/finclude/petscsys.h>

  use petscsys

  implicit none
  private

  type, public :: mpi_type
     MPI_Comm, public :: comm !! MPI communicator
     PetscMPIInt, public :: size !! number of processors
     PetscMPIInt, public :: rank !! processor rank
     PetscMPIInt, public :: input_rank = 0 !! Processor rank for input handling
     PetscMPIInt, public :: output_rank = 0 !! Processor rank for output handling
   contains
     private
     procedure, public :: init => mpi_initialize
     procedure, public :: broadcast_error_flag => mpi_broadcast_error_flag
     procedure, public :: broadcast_logical => mpi_broadcast_logical
  end type mpi_type

  type(mpi_type), public :: mpi

contains

!------------------------------------------------------------------------

  subroutine mpi_initialize(self, comm)
    !! Initializes MPI object. (Note this does not call MPI_INIT- that is
    !! already done by PETSC_INITIALIZE.)

    class(mpi_type), intent(in out) :: self
    MPI_Comm, intent(in) :: comm !! MPI communicator
    ! Locals:
    PetscInt :: ierr

    self%comm = comm
    call MPI_COMM_RANK(self%comm, self%rank, ierr)
    call MPI_COMM_SIZE(self%comm, self%size, ierr)
    
  end subroutine mpi_initialize

!------------------------------------------------------------------------

  subroutine mpi_broadcast_error_flag(self, err)
    !! If integer error flag err is greater than zero on any rank,
    !! broadcast value 1 to all ranks.

    class(mpi_type), intent(in out) :: self
    PetscErrorCode, intent(in out) :: err
    ! Locals:
    PetscBool :: any_err
    PetscErrorCode :: ierr

    call MPI_reduce(err > 0, any_err, 1, MPI_LOGICAL, MPI_LOR, &
         self%input_rank, self%comm, ierr)
    if (self%rank == self%input_rank) then
       if (any_err) then
          err = 1
       else
          err = 0
       end if
    end if
    call MPI_bcast(err, 1, MPI_INTEGER, self%input_rank, self%comm, ierr)

  end subroutine mpi_broadcast_error_flag

!------------------------------------------------------------------------

  subroutine mpi_broadcast_logical(self, val)
    !! If logical val is true on any rank, broadcast to all ranks.

    class(mpi_type), intent(in out) :: self
    PetscBool, intent(in out) :: val
    ! Locals:
    PetscBool :: any_val
    PetscErrorCode :: ierr

    call MPI_reduce(val, any_val, 1, MPI_LOGICAL, MPI_LOR, &
         self%input_rank, self%comm, ierr)
    if (self%rank == self%input_rank) then
       val = any_val
    end if
    call MPI_bcast(val, 1, MPI_LOGICAL, self%input_rank, self%comm, ierr)

  end subroutine mpi_broadcast_logical

!------------------------------------------------------------------------

end module mpi_module
