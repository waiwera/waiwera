module logfile_module
  !! Module for logging output messages to file.

  use mpi_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  PetscInt, parameter, public :: LOG_LEVEL_INFO  = 1, &
       LOG_LEVEL_WARN = 2, LOG_LEVEL_ERR = 3
  PetscInt, parameter :: max_log_level_name_length = 4
  character(max_log_level_name_length), parameter :: &
       log_level_name(3) = ['info', 'warn', 'err ']
  PetscInt, parameter, public :: max_logfile_name_length = 120
  PetscInt, parameter :: max_log_key_length = 16

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
    !! Initialise logfile. If filename is empty, no disk file is
    !! created, but echoing messages to console is still possible.

    class(logfile_type), intent(in out) :: self
    character(*), intent(in) :: filename !! Log file name
    ! Locals:
    PetscErrorCode :: ierr

    self%filename = filename
    if (self%filename /= "") then
       call PetscViewerASCIIOpen(mpi%comm, filename, self%viewer, ierr)
       CHKERRQ(ierr)
       call PetscViewerASCIIPushSynchronized(self%viewer, ierr)
       CHKERRQ(ierr)
    end if
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

    if (self%filename /= "") then
       call PetscViewerASCIISynchronizedPrintf(self%viewer, string, ierr)
       CHKERRQ(ierr)
    end if

    if (echo) then
       call PetscViewerASCIISynchronizedPrintf(PETSC_VIEWER_STDOUT_WORLD, &
            string, ierr); CHKERRQ(ierr)
    end if

  end subroutine logfile_write_string

!------------------------------------------------------------------------

  subroutine logfile_write(self, level, source, event, int_keys, &
       int_values, real_keys, real_values, echo)
    !! Write message to logfile, optionally echoing to console output.

    class(logfile_type), intent(in out) :: self
    PetscInt, intent(in) :: level
    character(*), intent(in) :: source, event
    character(max_log_key_length), intent(in), optional :: int_keys(:), &
         real_keys(:)
    PetscInt, intent(in), optional :: int_values(:)
    PetscReal, intent(in), optional :: real_values(:)
    PetscBool, intent(in), optional :: echo
    ! Locals:
    PetscBool :: do_echo
    PetscBool, parameter :: default_echo = PETSC_FALSE
    character(:), allocatable ::  content, msg
    PetscInt, parameter :: max_str_len = 12, num_real_digits = 6
    character(max_str_len) :: int_str, real_str
    character(7) :: int_fmt, real_fmt
    PetscInt :: i, num_int, num_real
    character, parameter :: lf = new_line('a')

    if (present(echo)) then
       do_echo = echo
    else
       do_echo = default_echo
    end if

    content = ""
    if (present(int_keys) .and. present(int_values)) then
       num_int = size(int_keys)
       write(int_fmt, '(a,i2,a)') '(i', max_str_len, ')'
       do i = 1, num_int
          write(int_str, int_fmt) int_values(i)
          content = content  // trim(int_keys(i)) // &
               ': ' // trim(adjustl(int_str))
          if (i < num_int) then
             content = content // ', '
          end if
       end do
       if (present(real_keys) .and. present(real_values)) then
          content = content // ', '
       end if
    end if
    if (present(real_keys) .and. present(real_values)) then
       num_real = size(real_keys)
       write(real_fmt, '(a,i2,a,i1,a)') '(e', max_str_len, '.', &
            num_real_digits, ')'
       do i = 1, num_real
          write(real_str, real_fmt) real_values(i)
          content = content // trim(real_keys(i)) // &
               ': ' // trim(adjustl(real_str))
          if (i < num_real) then
             content = content // ', '
          end if
       end do
    end if

    msg = trim(log_level_name(level)) // ' ' &
         // trim(source) // ': ' &
         // trim(event) // ' [' &
         // trim(content) // ']' // lf

    call self%write_string(msg, do_echo)

  end subroutine logfile_write

!------------------------------------------------------------------------

  subroutine logfile_destroy(self)
    !! Destroy logfile.

    class(logfile_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    if (self%filename /= "") then
       call PetscViewerASCIIPopSynchronized(self%viewer, ierr)
       CHKERRQ(ierr)
       call PetscViewerDestroy(self%viewer, ierr); CHKERRQ(ierr)
    end if
    self%filename = ""
    call PetscViewerASCIIPopSynchronized(PETSC_VIEWER_STDOUT_WORLD, &
         ierr); CHKERRQ(ierr)

  end subroutine logfile_destroy

!------------------------------------------------------------------------

end module logfile_module
