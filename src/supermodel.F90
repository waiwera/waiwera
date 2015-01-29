program supermodel

  ! Main supermodel driver code.

  use simulation_module

  implicit none

#include <petsc-finclude/petscsys.h>

  type(simulation_type) :: sim
  PetscErrorCode :: ierr
  PetscInt, parameter :: max_filename_length = 200
  character(max_filename_length) :: filename

  call output_program_info()

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)

  call get_filename(filename)

  call sim%init(filename)
  call sim%run()
  call sim%destroy()

  call PetscFinalize(ierr); CHKERRQ(ierr)

contains

!------------------------------------------------------------------------

  subroutine output_program_info()

    ! Output program info.

    write (*,*) 'Supermodel version 0.001'

  end subroutine output_program_info

!------------------------------------------------------------------------

  subroutine get_filename(filename)

    ! Gets filename from program argument or input.

    character(*), intent(out) :: filename
    ! Locals:
    PetscInt :: num_args
    
    num_args = command_argument_count()
    if (num_args == 0) then
       write (*,*) 'Input file:'
       read (*,'(a120)') filename
    else 
       call get_command_argument(1, filename)
       filename = trim(filename)
    end if

  end subroutine get_filename

!------------------------------------------------------------------------

end program supermodel
