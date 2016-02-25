program supermodel
  !! Main supermodel driver program.

  use mpi_module
  use fson
  use fson_mpi_module
  use flow_simulation_module
  use timestepper_module
  use logfile_module

  implicit none

#include <petsc/finclude/petscsys.h>

  type(fson_value), pointer :: json
  type(flow_simulation_type) :: simulation
  type(timestepper_type) :: timestepper
  character(max_flow_simulation_filename_length) :: filename
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
  call mpi%init(PETSC_COMM_WORLD)

  call get_filename(filename)
  json => fson_parse_mpi(filename)

  call simulation%init(json, filename)
  call timestepper%init(json, simulation)
  call fson_destroy_mpi(json)
  call simulation%input_summary()

  call timestepper%run()

  call timestepper%destroy()
  call simulation%destroy()

  call PetscFinalize(ierr); CHKERRQ(ierr)

contains

!------------------------------------------------------------------------

  subroutine get_filename(filename)
    !! Gets filename from program argument or input.

    character(*), intent(out) :: filename
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

end program supermodel
