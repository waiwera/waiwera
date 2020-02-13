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

program waiwera
  !! Main Waiwera driver program.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use fson
  use fson_mpi_module
  use flow_simulation_module
  use timestepper_module
  use profiling_module

  implicit none

  type(fson_value), pointer :: json !! JSON object for simulation input
  type(flow_simulation_type) :: simulation !! Flow simulation object
  type(timestepper_type) :: timestepper !! Timestepper for time-stepping the simulation
  character(:), allocatable :: filename !! JSON input filename
  PetscErrorCode :: ierr !! Error code

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRA(ierr)
  call init_profiling()

  call process_arguments(filename)

  if (allocated(filename)) then

     json => fson_parse_mpi(filename)

     call simulation%init(json, filename, ierr)

     if (ierr == 0) then
        call timestepper%init(json, simulation)
        call fson_destroy_mpi(json)
        call simulation%input_summary()
        call timestepper%input_summary()
        call simulation%log_statistics()
        call timestepper%run()
        call timestepper%destroy()
     else
        call fson_destroy_mpi(json)
     end if

     call simulation%destroy()

  end if

  call PetscFinalize(ierr); CHKERRA(ierr)

contains

!------------------------------------------------------------------------

  subroutine waiwera_help()
    !! Waiwera command line help.

    write(*,'(a)') 'Usage:'
    write(*,'(a)') '  waiwera --help | -h | --version | -v'
    write(*,'(a)') '  waiwera <filename> [PETSc-option...]'

    write(*,'(a)')
    write(*,'(a)') 'Options:'
    write(*,'(a)') '  --help          Show this help.'
    write(*,'(a)') '  -h              Show this help, together with PETSc help.'
    write(*,'(a)') '  --version       Show Waiwera version.'
    write(*,'(a)') '  -v              Show Waiwera version, togther with PETSc version.'
    write(*,'(a)') '  filename        Name of JSON input file.'
    write(*,'(a)') '  PETSc-option    PETSc command line options.'

  end subroutine waiwera_help

!------------------------------------------------------------------------

  subroutine process_arguments(filename)
    !! Process command line arguments, returning usage message if no
    !! arguments are specified, otherwise returning filename.

    use version_module, only: waiwera_version

    character(:), allocatable, intent(out) :: filename !! JSON input filename
    ! Locals:
    PetscInt :: num_args, len_arg, len_filename, ierr
    character(:), allocatable :: arg
    PetscMPIInt :: rank

    len_filename = 0
    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       num_args = command_argument_count()
       if (num_args == 0) then
          call waiwera_help()
       else
          call get_command_argument(1, length = len_arg)
          allocate(character(len_arg) :: arg)
          call get_command_argument(1, value = arg)
          select case (arg)
          case ('--help', '-h')
             call waiwera_help()
          case ('--version', '-v')
             write(*, '(a)') waiwera_version
          case default
             filename = arg
             len_filename = len_arg
          end select
          deallocate(arg)
       end if

    end if

    ! Broadcast filename:
    call MPI_bcast(len_filename, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    if (len_filename > 0) then
       if (rank > 0) then
          allocate(character(len_filename) :: filename)
       end if
       call MPI_bcast(filename, len_filename, MPI_CHARACTER, &
            0, PETSC_COMM_WORLD, ierr)
    end if

  end subroutine process_arguments

!------------------------------------------------------------------------

end program waiwera
