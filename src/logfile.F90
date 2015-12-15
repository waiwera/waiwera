module logfile_module
  !! Module for logging output messages to file.

  use mpi_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  PetscInt, parameter, public :: LOG_LEVEL_INFO  = 1, &
       LOG_LEVEL_WARN = 2, LOG_LEVEL_ERR = 3
  PetscInt, parameter :: max_log_level_name_length = 7
  character(max_log_level_name_length), parameter :: &
       log_level_name(3) = ['info   ', 'warning', 'error  ']
  PetscInt, parameter, public :: max_logfile_name_length = 120

  type logfile_type
     private
     PetscViewer :: viewer
     character(max_logfile_name_length), public :: filename
   contains
     private
     procedure, public :: init => logfile_init
     procedure, public :: write_string => logfile_write_string
     procedure, public :: write => logfile_write
     procedure, public :: destroy => logfile_destroy
  end type logfile_type

  public :: logfile_type

contains

!------------------------------------------------------------------------

  subroutine logfile_init(self, filename)
    !! Initialise logfile.

    class(logfile_type), intent(in out) :: self
    character(*), intent(in) :: filename !! Log file name
    ! Locals:
    PetscErrorCode :: ierr

    self%filename = filename
    call PetscViewerASCIIOpen(mpi%comm, filename, self%viewer, ierr)
    CHKERRQ(ierr)
    call PetscViewerASCIIPushSynchronized(self%viewer, ierr)
    CHKERRQ(ierr)
    call PetscViewerASCIIPushSynchronized(PETSC_VIEWER_STDOUT_WORLD, &
         ierr); CHKERRQ(ierr)

  end subroutine logfile_init

!------------------------------------------------------------------------

  subroutine logfile_write_string(self, string, echo)
    !! Write string to logfile, optionally echoing to console output.

    class(logfile_type), intent(in out) :: self
    character(*), intent(in) :: string
    PetscBool, intent(in) :: echo
    ! Locals:
    PetscErrorCode :: ierr

    call PetscViewerASCIISynchronizedPrintf(self%viewer, string, ierr)
    CHKERRQ(ierr)

    if (echo) then
       call PetscViewerASCIISynchronizedPrintf(PETSC_VIEWER_STDOUT_WORLD, &
            string, ierr); CHKERRQ(ierr)
    end if

  end subroutine logfile_write_string

!------------------------------------------------------------------------

  subroutine logfile_write(self, level, tag, content, echo)
    !! Write message to logfile, optionally echoing to console output.

    class(logfile_type), intent(in out) :: self
    PetscInt, intent(in) :: level
    character(*), intent(in) :: tag
    character(*), intent(in) :: content
    PetscBool, intent(in), optional :: echo
    ! Locals:
    PetscBool :: do_echo
    PetscBool, parameter :: default_echo = PETSC_FALSE
    character(:), allocatable ::  msg
    character, parameter :: lf = new_line('a')

    if (present(echo)) then
       do_echo = echo
    else
       do_echo = default_echo
    end if

    msg = trim(log_level_name(level)) // ' ' &
         // trim(tag) // ' ' // trim(content) // lf

    call self%write_string(msg, do_echo)

  end subroutine logfile_write

!------------------------------------------------------------------------

  subroutine logfile_destroy(self)
    !! Destroy logfile.

    class(logfile_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    self%filename = ""
    call PetscViewerASCIIPopSynchronized(self%viewer, ierr)
    CHKERRQ(ierr)
    call PetscViewerASCIIPopSynchronized(PETSC_VIEWER_STDOUT_WORLD, &
         ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(self%viewer, ierr); CHKERRQ(ierr)

  end subroutine logfile_destroy

!------------------------------------------------------------------------

end module logfile_module
