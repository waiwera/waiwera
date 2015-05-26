program supermodel
  !! Main supermodel driver program.

  use mpi_module
  use fson
  use flow_simulation_module
  use timestepper_module

  implicit none

#include <petsc-finclude/petscsys.h>

  type(fson_value), pointer :: json
  type(flow_simulation_type) :: sim
  type(timestepper_type) :: ts
  PetscInt, parameter, public :: max_filename_length = 200
  character(max_filename_length) :: filename
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
  call mpi%init(PETSC_COMM_WORLD)

  call output_program_info()

  call get_filename(filename)
  if (mpi%rank == mpi%input_rank) json => fson_parse(filename)
  call sim%init(json)
  call ts%init(json, sim)
  if (mpi%rank == mpi%input_rank) call fson_destroy(json)

  call ts%run()

  call ts%destroy()
  call sim%destroy()

  call PetscFinalize(ierr); CHKERRQ(ierr)

contains

!------------------------------------------------------------------------

  subroutine output_program_info()
    !! Outputs program information.

    if (mpi%rank == mpi%output_rank) then
       write (*,*) 'Supermodel version 0.001'
    end if

  end subroutine output_program_info

!------------------------------------------------------------------------

  subroutine get_filename(filename)
    !! Gets filename from program argument or input.

    character(*), intent(out) :: filename
    ! Locals:
    PetscInt :: num_args

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

  end subroutine get_filename

!------------------------------------------------------------------------

end program supermodel
