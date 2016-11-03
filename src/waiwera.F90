!   Copyright 2016 University of Auckland.

!   This file is part of Waiwera.

!   Waiwera is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.

!   Waiwera is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.

!   You should have received a copy of the GNU General Public License
!   along with Waiwera.  If not, see <http://www.gnu.org/licenses/>.

program waiwera
  !! Main Waiwera driver program.

  use mpi_module
  use fson
  use fson_mpi_module
  use flow_simulation_module
  use timestepper_module
  use profiling_module

  implicit none

#include <petsc/finclude/petscsys.h>

  type(fson_value), pointer :: json !! JSON object for simulation input
  type(flow_simulation_type) :: simulation !! Flow simulation object
  type(timestepper_type) :: timestepper !! Timestepper for time-stepping the simulation
  character(max_flow_simulation_filename_length) :: filename !! JSON input filename
  PetscErrorCode :: ierr !! Error code for PETSc

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
  call mpi%init(PETSC_COMM_WORLD)

  call init_profiling()

  call get_filename(filename)
  json => fson_parse_mpi(filename)

  call simulation%init(json, filename)
  call timestepper%init(json, simulation)
  call fson_destroy_mpi(json)
  call simulation%input_summary()
  call timestepper%input_summary()

  call timestepper%run()

  call timestepper%destroy()
  call simulation%destroy()

  call PetscFinalize(ierr); CHKERRQ(ierr)

contains

!------------------------------------------------------------------------

  subroutine get_filename(filename)
    !! Gets filename from program argument or input.

    character(*), intent(out) :: filename !! JSON input filename
    ! Locals:
    PetscInt :: num_args, filename_length, ierr

    if (mpi%rank == mpi%input_rank) then

       num_args = command_argument_count()
       if (num_args == 0) then
          write (*,*) 'Input file:'
          read (*,'(a120)') filename
       else 
          call get_command_argument(1, filename)
          filename = trim(filename)
       end if

    end if
    filename_length = len_trim(filename)
    call MPI_bcast(filename, filename_length, MPI_CHARACTER, &
         mpi%input_rank, mpi%comm, ierr)

  end subroutine get_filename

!------------------------------------------------------------------------

end program waiwera
