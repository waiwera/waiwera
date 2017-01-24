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
  !! MPI variables.

#include <petsc/finclude/petscsys.h>

  use petscsys

  implicit none
  private

contains

!------------------------------------------------------------------------

  subroutine mpi_broadcast_error_flag(err)
    !! If integer error flag err is greater than zero on any rank,
    !! broadcast value 1 to all ranks.

    PetscErrorCode, intent(in out) :: err
    ! Locals:
    PetscBool :: any_err
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_reduce(err > 0, any_err, 1, MPI_LOGICAL, MPI_LOR, &
         0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       if (any_err) then
          err = 1
       else
          err = 0
       end if
    end if
    call MPI_bcast(err, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)

  end subroutine mpi_broadcast_error_flag

!------------------------------------------------------------------------

  subroutine mpi_broadcast_logical(val)
    !! If logical val is true on any rank, broadcast to all ranks.

    PetscBool, intent(in out) :: val
    ! Locals:
    PetscBool :: any_val
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    call MPI_reduce(val, any_val, 1, MPI_LOGICAL, MPI_LOR, &
         0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       val = any_val
    end if
    call MPI_bcast(val, 1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr)

  end subroutine mpi_broadcast_logical

!------------------------------------------------------------------------

end module mpi_utils_module
