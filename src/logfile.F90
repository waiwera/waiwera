module logfile_module
  !! Module for logging output messages to file.

  use mpi_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  PetscInt, parameter, public :: LOG_LEVEL_INFO  = 1, &
       LOG_LEVEL_WARN = 2, LOG_LEVEL_ERR = 3
  PetscInt, parameter, public :: max_log_message_string_length = 80
  PetscInt, parameter :: max_log_level_name_length = 4
  character(max_log_level_name_length), parameter :: &
       log_level_name(3) = ['info', 'warn', 'err ']
  PetscInt, parameter, public :: max_logfile_name_length = 120
  PetscInt, parameter :: max_log_key_length = 16
  PetscInt, parameter :: max_log_number_str_len = 12
  PetscInt, parameter :: num_log_real_digits = 6
  character, parameter :: lf = new_line('a')

  type logfile_type
     private
     PetscViewer :: viewer
     character(max_logfile_name_length), public :: filename
   contains
     private
     procedure, public :: init => logfile_init
     procedure, public :: write_string => logfile_write_string
     procedure, public :: append_int_data => logfile_append_int_data
     procedure, public :: append_real_data => logfile_append_real_data
     procedure, public :: append_real_array_data => &
          logfile_append_real_array_data
     procedure, public :: append_string_data => logfile_append_string_data
     procedure, public :: write => logfile_write
     procedure, public :: write_blank => logfile_write_blank
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

  subroutine logfile_append_int_data(self, int_keys, int_values, content)
    !! Append integer data to logfile content line.

    class(logfile_type), intent(in out) :: self
    character(*), intent(in) :: int_keys(:)
    PetscInt, intent(in) :: int_values(:)
    character(:), allocatable, intent(in out) ::  content
    ! Locals:
    PetscInt :: num_int, i
    character(7) :: int_fmt
    character(max_log_number_str_len) :: int_str

    if (len(content) > 0) then
       content = content // ', '
    end if

    num_int = size(int_keys)
    write(int_fmt, '(a,i2,a)') '(i', max_log_number_str_len, ')'
    do i = 1, num_int
       write(int_str, int_fmt) int_values(i)
       content = content  // trim(int_keys(i)) // &
            ': ' // trim(adjustl(int_str))
       if (i < num_int) then
          content = content // ', '
       end if
    end do

  end subroutine logfile_append_int_data

!------------------------------------------------------------------------

  subroutine logfile_append_real_data(self, real_keys, real_values, content)
    !! Append real data to logfile content line.

    class(logfile_type), intent(in out) :: self
    character(*), intent(in) :: real_keys(:)
    PetscReal, intent(in) :: real_values(:)
    character(:), allocatable, intent(in out) ::  content
    ! Locals:
    PetscInt :: num_real, i
    character(7) :: real_fmt
    character(max_log_number_str_len) :: real_str

    if (len(content) > 0) then
       content = content // ', '
    end if

    num_real = size(real_keys)
    write(real_fmt, '(a,i2,a,i1,a)') '(e', max_log_number_str_len, &
         '.', num_log_real_digits, ')'
    do i = 1, num_real
       write(real_str, real_fmt) real_values(i)
       content = content // trim(real_keys(i)) // &
            ': ' // trim(adjustl(real_str))
       if (i < num_real) then
          content = content // ', '
       end if
    end do

  end subroutine logfile_append_real_data

!------------------------------------------------------------------------

  subroutine logfile_append_real_array_data(self, real_array_key, &
       real_array_value, content)
    !! Append real array data to logfile content line. (Only one array
    !! can be appended.)

    class(logfile_type), intent(in out) :: self
    character(*), intent(in) :: real_array_key
    PetscReal, intent(in) :: real_array_value(:)
    character(:), allocatable, intent(in out) ::  content
    ! Locals:
    PetscInt :: num_real, i
    character(7) :: real_fmt
    character(max_log_number_str_len) :: real_str

    if (len(content) > 0) then
       content = content // ', '
    end if

    num_real = size(real_array_value)
    write(real_fmt, '(a,i2,a,i1,a)') '(e', &
         max_log_number_str_len, '.', num_log_real_digits, ')'
    content = content // trim(real_array_key) // ': ['
    do i = 1, num_real
       write(real_str, real_fmt) real_array_value(i)
       content = content // trim(adjustl(real_str))
       if (i < num_real) then
          content = content // ', '
       end if
    end do
    content = content // ']'


  end subroutine logfile_append_real_array_data

!------------------------------------------------------------------------

  subroutine logfile_append_string_data(self, str_key, str_value, content)
    !! Append string data to logfile content line. (Only one string
    !! can be appended.)

    class(logfile_type), intent(in out) :: self
    character(*), intent(in) :: str_key
    character(*), intent(in) :: str_value
    character(:), allocatable, intent(in out) ::  content

    if (len(content) > 0) then
       content = content // ', '
    end if

    content = content // trim(str_key) // ': ' // trim(str_value)

  end subroutine logfile_append_string_data

!------------------------------------------------------------------------

  subroutine logfile_write(self, level, source, event, int_keys, &
       int_values, real_keys, real_values, str_key, str_value, &
       real_array_key, real_array_value, echo)
    !! Write message to logfile, optionally echoing to console output.
    !! Output is in YAML format: each message is formatted as an
    !! inline list, with the first three elements being the level,
    !! source and event respectively. If there are additional integer,
    !! real, string or real array data, these are appended as an
    !! associative array.

    class(logfile_type), intent(in out) :: self
    PetscInt, intent(in) :: level
    character(*), intent(in) :: source, event
    character(*), intent(in), optional :: int_keys(:), &
         real_keys(:), str_key, real_array_key
    PetscInt, intent(in), optional :: int_values(:)
    PetscReal, intent(in), optional :: real_values(:)
    character(*), intent(in), optional :: str_value
    PetscReal, intent(in), optional :: real_array_value(:)
    PetscBool, intent(in), optional :: echo
    ! Locals:
    PetscBool :: do_echo, has_data
    PetscBool, parameter :: default_echo = PETSC_FALSE
    character(:), allocatable ::  content, msg

    if (present(echo)) then
       do_echo = echo
    else
       do_echo = default_echo
    end if

    content = ""
    has_data = PETSC_FALSE

    if (present(int_keys) .and. present(int_values)) then
       has_data = PETSC_TRUE
       call self%append_int_data(int_keys, int_values, content)
    end if

    if (present(real_keys) .and. present(real_values)) then
       has_data = PETSC_TRUE
       call self%append_real_data(real_keys, real_values, content)
    end if

    if (present(real_array_key) .and. present(real_array_value)) then
       has_data = PETSC_TRUE
       call self%append_real_array_data(real_array_key, &
            real_array_value, content)
    end if

    if (present(str_key) .and. present(str_value)) then
       has_data = PETSC_TRUE
       call self%append_string_data(str_key, str_value, content)
    end if

    msg = '- [' // trim(log_level_name(level)) // ', ' &
         // trim(source) // ', ' // trim(event)
    if (has_data) then
       msg = msg // ', {' // trim(content) // '}]' // lf
    else
       msg = msg // ']' // lf
    end if

    call self%write_string(msg, do_echo)

  end subroutine logfile_write

!------------------------------------------------------------------------

  subroutine logfile_write_blank(self)
    !! Writes blank line to logfile. 

    class(logfile_type), intent(in out) :: self

    call self%write_string(lf, PETSC_TRUE)

  end subroutine logfile_write_blank

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
