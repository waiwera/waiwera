program supermodel
  !! Main supermodel driver program.

  use mpi_module
  use simulation_module

  implicit none

#include <petsc-finclude/petscsys.h>

  type(simulation_type) :: sim !! Simulation
  character(max_filename_length) :: filename !! filename
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
  call mpi%init(PETSC_COMM_WORLD)

  call output_program_info()

  call get_filename(filename)

  call sim%init(filename)
  call sim%run()
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
    integer :: num_args

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
