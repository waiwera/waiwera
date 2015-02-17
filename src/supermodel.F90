program supermodel
  !! Main supermodel driver program.

  use simulation_module

  implicit none

#include <petsc-finclude/petscsys.h>

  type(simulation_type) :: sim !! Simulation
  integer, parameter :: max_filename_length = 200 !! maximum filename length
  character(max_filename_length) :: filename !! filename
  PetscErrorCode :: ierr
  MPI_Comm :: comm

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
  comm = PETSC_COMM_WORLD

  call output_program_info()

  call get_filename(filename)

  call sim%init(filename, comm)
  call sim%run()
  call sim%destroy()

  call PetscFinalize(ierr); CHKERRQ(ierr)

contains

!------------------------------------------------------------------------

  subroutine output_program_info()
    !! Outputs program information.

    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_comm_rank(comm, rank, ierr); CHKERRQ(ierr)

    if (rank == 0) then
       write (*,*) 'Supermodel version 0.001'
    end if

  end subroutine output_program_info

!------------------------------------------------------------------------

  subroutine get_filename(filename)
    !! Gets filename from program argument or input.

    character(*), intent(out) :: filename
    ! Locals:
    PetscMPIInt :: rank
    PetscErrorCode :: ierr
    integer :: num_args

    call MPI_comm_rank(comm, rank, ierr); CHKERRQ(ierr)

    if (rank == 0) then

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
