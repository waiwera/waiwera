module logfile_module
  !! Module for logging output messages to file on MPI output process.

  use mpi_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  PetscInt, parameter, public :: LOG_LEVEL_INFO  = 1, &
       LOG_LEVEL_WARN = 2, LOG_LEVEL_ERR = 3
  PetscInt, parameter :: max_log_level_name_length = 7
  character(max_log_level_name_length) :: log_level_name(3) = &
       ['info   ', 'warning', 'error  ']

  PetscInt, parameter, public :: max_logfile_name_length = 120
  PetscInt, parameter, public :: max_message_length = 80
  character(20), public :: message_format = '(a, 1x, a, 1x, a)'

  type logfile_type
     private
     PetscInt :: unit
     character(max_logfile_name_length), public :: filename
   contains
     private
     procedure, public :: init => logfile_init
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

    self%unit = 0
    self%filename = filename

    if (mpi%rank == mpi%output_rank) then
       open(newunit = self%unit, file = self%filename, &
            status = 'replace')
    end if

  end subroutine logfile_init

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
    character(max_message_length) :: msg
    PetscInt :: ierr

    if (present(echo)) then
       do_echo = echo
    else
       do_echo = default_echo
    end if

    write(msg, message_format) log_level_name(level), tag, content

    if (mpi%rank /= mpi%output_rank) then

       call MPI_send(msg, max_message_length, MPI_CHARACTER, mpi%output_rank, 0, &
            mpi%comm, ierr)
       call MPI_send(do_echo, 1, MPI_LOGICAL, mpi%output_rank, 0, &
            mpi%comm, ierr)

    else

       call MPI_recv(msg, max_message_length, MPI_CHARACTER, MPI_ANY_SOURCE, 0, &
            mpi%comm, ierr)
       call MPI_recv(do_echo, 1, MPI_LOGICAL, MPI_ANY_SOURCE, 0, &
            mpi%comm, ierr)

       write(self%unit, message_format) msg

       if (do_echo) then
          write(*, message_format) msg
       end if

    end if

  end subroutine logfile_write

!------------------------------------------------------------------------

  subroutine logfile_destroy(self)
    !! Destroy logfile.

    class(logfile_type), intent(in out) :: self

    if (mpi%rank == mpi%output_rank) then
       close(self%unit)
       self%unit = 0
    end if
    self%filename = ""

  end subroutine logfile_destroy

!------------------------------------------------------------------------

end module logfile_module
