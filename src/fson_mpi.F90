module fson_mpi_module
  !! Subroutines for JSON input under MPI, with input on one rank
  !! broadcast to other ranks, and with the ability to specify
  !! default values if not present in the JSON input.

  use kinds_module
  use mpi_module
  use logfile_module
  use fson

  implicit none
  private

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdef.h>

  interface fson_get_default
     module procedure fson_get_default_integer
     module procedure fson_get_default_real
     module procedure fson_get_default_double
     module procedure fson_get_default_logical
     module procedure fson_get_default_character
     module procedure fson_get_default_array_1d_integer
     module procedure fson_get_default_array_1d_real
     module procedure fson_get_default_array_1d_double
     module procedure fson_get_default_array_1d_logical
     module procedure fson_get_default_array_2d_integer
     module procedure fson_get_default_array_2d_real
     module procedure fson_get_default_array_2d_double
     module procedure fson_get_default_array_2d_logical
  end interface fson_get_default

  interface fson_get_mpi
     module procedure fson_get_mpi_fson_value
     module procedure fson_get_mpi_integer
     module procedure fson_get_mpi_real
     module procedure fson_get_mpi_double
     module procedure fson_get_mpi_logical
     module procedure fson_get_mpi_character
     module procedure fson_get_mpi_array_1d_integer
     module procedure fson_get_mpi_array_1d_real
     module procedure fson_get_mpi_array_1d_double
     module procedure fson_get_mpi_array_1d_logical
     module procedure fson_get_mpi_array_2d_integer
     module procedure fson_get_mpi_array_2d_real
     module procedure fson_get_mpi_array_2d_double
     module procedure fson_get_mpi_array_2d_logical
  end interface fson_get_mpi

  public :: fson_get_default, fson_get_mpi, fson_has_mpi
  public :: fson_type_mpi, fson_value_count_mpi, fson_value_get_mpi
  public :: fson_parse_mpi, fson_destroy_mpi

contains

!------------------------------------------------------------------------

  PetscBool function assoc_non_null(self)
    !! Returns .true. if fson_value self is associated and has value
    !! type non-null.

    use fson_value_m, only : TYPE_NULL

    type(fson_value), pointer, intent(in) :: self

    assoc_non_null = .false.
    if (associated(self)) then
       if (self%value_type /= TYPE_NULL) then
          assoc_non_null = .true.
       end if
    end if

  end function assoc_non_null

!------------------------------------------------------------------------
! fson_get_default routines
!------------------------------------------------------------------------

  subroutine fson_get_default_integer(self, path, default, val, &
       logfile, log_key)
    !! Gets integer value with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscInt, intent(in) :: default
    PetscInt, intent(out) :: val
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(max_log_key_length) :: key

    call fson_get(self, path, p)
    if (assoc_non_null(p)) then
       call fson_get(p, ".", val)
    else
       val = default
       if (present(logfile)) then
          if (logfile%active) then
             if (present(log_key)) then
                key = log_key
             else
                key = path
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  int_keys = [key], &
                  int_values = [default])
          end if
       end if
    end if

  end subroutine fson_get_default_integer
  
!------------------------------------------------------------------------
  
  subroutine fson_get_default_real(self, path, default, val, logfile, &
       log_key)
    !! Gets real value with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    real, intent(in) :: default
    real, intent(out) :: val
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(max_log_key_length) :: key

    call fson_get(self, path, p)
    if (assoc_non_null(p)) then
       call fson_get(p, ".", val)
    else
       val = default
       if (present(logfile)) then
          if (logfile%active) then
             if (present(log_key)) then
                key = log_key
             else
                key = path
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  real_keys = [key], &
                  real_values = [dble(default)])
          end if
       end if
    end if

  end subroutine fson_get_default_real

!------------------------------------------------------------------------

  subroutine fson_get_default_double(self, path, default, val, logfile, &
       log_key)
    !! Gets double value with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscReal, intent(in) :: default
    PetscReal, intent(out) :: val
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(max_log_key_length) :: key

    call fson_get(self, path, p)
    if (assoc_non_null(p)) then
       call fson_get(p, ".", val)
    else
       val = default
       if (present(logfile)) then
          if (logfile%active) then
             if (present(log_key)) then
                key = log_key
             else
                key = path
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  real_keys = [key], &
                  real_values = [default])
          end if
       end if
    end if

  end subroutine fson_get_default_double

!------------------------------------------------------------------------

  subroutine fson_get_default_logical(self, path, default, val, &
       logfile, log_key)
    !! Gets logical value with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscBool, intent(in) :: default
    PetscBool, intent(out) :: val
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(max_log_key_length) :: key

    call fson_get(self, path, p)
    if (assoc_non_null(p)) then
       call fson_get(p, ".", val)
    else
       val = default
       if (present(logfile)) then
          if (logfile%active) then
             if (present(log_key)) then
                key = log_key
             else
                key = path
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  logical_keys = [key], &
                  logical_values = [default])
          end if
       end if
    end if

  end subroutine fson_get_default_logical

!------------------------------------------------------------------------

  subroutine fson_get_default_character(self, path, default, val, &
       logfile, log_key)
    !! Gets character value with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: default
    character(len=*), intent(out) :: val
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(max_log_key_length) :: key

    call fson_get(self, path, p)
    if (assoc_non_null(p)) then
       call fson_get(p, ".", val)
    else
       val = default
       if (present(logfile)) then
          if (logfile%active) then
             if (present(log_key)) then
                key = log_key
             else
                key = path
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  str_key = key, str_value = default)
          end if
       end if
    end if

  end subroutine fson_get_default_character

!------------------------------------------------------------------------

  subroutine fson_get_default_array_1d_integer(self, path, default, &
       val, logfile, log_key)
    !! Gets 1-D integer array with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscInt, intent(in) :: default(:)
    PetscInt, allocatable, intent(out) :: val(:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(:), allocatable :: intstr, str
    character(max_log_key_length) :: key

    call fson_get(self, path, p)
    if (assoc_non_null(p)) then
       call fson_get(p, ".", val)
    else
       val = default
       if (present(logfile)) then
          if (logfile%active) then
             if (present(log_key)) then
                key = log_key
             else
                key = path
             end if
             if (size(val) > 0) then
                write(intstr, logfile%int_format) val(1)
                str = '[' // intstr // ',...]'
             else
                str = '[]'
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  str_key = key, str_value = str)
          end if
       end if
    end if

  end subroutine fson_get_default_array_1d_integer

!------------------------------------------------------------------------

  subroutine fson_get_default_array_1d_real(self, path, default, val, &
       logfile, log_key)
    !! Gets 1-D real array with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    real, intent(in) :: default(:)
    real, allocatable, intent(out) :: val(:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(:), allocatable :: realstr, str
    character(max_log_key_length) :: key

    call fson_get(self, path, p)
    if (assoc_non_null(p)) then
       call fson_get(p, ".", val)
    else
       val = default
       if (present(logfile)) then
          if (logfile%active) then
             if (present(log_key)) then
                key = log_key
             else
                key = path
             end if
             if (size(val) > 0) then
                write(realstr, logfile%real_format) val(1)
                str = '[' // realstr // ',...]'
             else
                str = '[]'
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  str_key = key, str_value = str)
          end if
       end if
    end if

  end subroutine fson_get_default_array_1d_real

!------------------------------------------------------------------------

  subroutine fson_get_default_array_1d_double(self, path, default, val, &
       logfile, log_key)
    !! Gets 1-D double precision array with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscReal, intent(in) :: default(:)
    PetscReal, allocatable, intent(out) :: val(:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(:), allocatable :: realstr, str
    character(max_log_key_length) :: key

    call fson_get(self, path, p)
    if (assoc_non_null(p)) then
       call fson_get(p, ".", val)
    else
       val = default
       if (present(logfile)) then
          if (logfile%active) then
             if (present(log_key)) then
                key = log_key
             else
                key = path
             end if
             if (size(val) > 0) then
                write(realstr, logfile%real_format) val(1)
                str = '[' // realstr // ',...]'
             else
                str = '[]'
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  str_key = key, str_value = str)
          end if
       end if
    end if

  end subroutine fson_get_default_array_1d_double

!------------------------------------------------------------------------

  subroutine fson_get_default_array_1d_logical(self, path, default, &
       val, logfile, log_key)
    !! Gets 1-D logical array with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscBool, intent(in) :: default(:)
    PetscBool, allocatable, intent(out) :: val(:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(:), allocatable :: logstr, str
    character(max_log_key_length) :: key

    call fson_get(self, path, p)
    if (assoc_non_null(p)) then
       call fson_get(p, ".", val)
    else
       val = default
       if (present(logfile)) then
          if (logfile%active) then
             if (present(log_key)) then
                key = log_key
             else
                key = path
             end if
             if (size(val) > 0) then
                write(logstr, '(L)') val(1)
                str = '[' // logstr // ',...]'
             else
                str = '[]'
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  str_key = key, str_value = str)
          end if
       end if
    end if

  end subroutine fson_get_default_array_1d_logical

!------------------------------------------------------------------------

  subroutine fson_get_default_array_2d_integer(self, path, default, &
       val, logfile, log_key)
    !! Gets 2-D integer array with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscInt, intent(in) :: default(:,:)
    PetscInt, allocatable, intent(out) :: val(:,:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(:), allocatable :: intstr, str
    character(max_log_key_length) :: key

    call fson_get(self, path, p)
    if (assoc_non_null(p)) then
       call fson_get(p, ".", val)
    else
       val = default
       if (present(logfile)) then
          if (logfile%active) then
             if (present(log_key)) then
                key = log_key
             else
                key = path
             end if
             if (size(val) > 0) then
                write(intstr, logfile%int_format) val(1,1)
                str = '[[' // intstr // ',...]]'
             else
                str = '[]'
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  str_key = key, str_value = str)
          end if
       end if
    end if

  end subroutine fson_get_default_array_2d_integer

!------------------------------------------------------------------------

  subroutine fson_get_default_array_2d_real(self, path, default, val, &
       logfile, log_key)
    !! Gets 2-D real array with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    real, intent(in) :: default(:,:)
    real, allocatable, intent(out) :: val(:,:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(:), allocatable :: realstr, str
    character(max_log_key_length) :: key

    call fson_get(self, path, p)
    if (assoc_non_null(p)) then
       call fson_get(p, ".", val)
    else
       val = default
       if (present(logfile)) then
          if (logfile%active) then
             if (present(log_key)) then
                key = log_key
             else
                key = path
             end if
             if (size(val) > 0) then
                write(realstr, logfile%real_format) val(1,1)
                str = '[[' // realstr // ',...]]'
             else
                str = '[]'
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  str_key = key, str_value = str)
          end if
       end if
    end if

  end subroutine fson_get_default_array_2d_real

!------------------------------------------------------------------------

  subroutine fson_get_default_array_2d_double(self, path, default, val, &
       logfile, log_key)
    !! Gets 2-D double array with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscReal, intent(in) :: default(:,:)
    PetscReal, allocatable, intent(out) :: val(:,:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(:), allocatable :: realstr, str
    character(max_log_key_length) :: key

    call fson_get(self, path, p)
    if (assoc_non_null(p)) then
       call fson_get(p, ".", val)
    else
       val = default
       if (present(logfile)) then
          if (logfile%active) then
             if (present(log_key)) then
                key = log_key
             else
                key = path
             end if
             if (size(val) > 0) then
                write(realstr, logfile%real_format) val(1,1)
                str = '[[' // realstr // ',...]]'
             else
                str = '[]'
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  str_key = key, str_value = str)
          end if
       end if
    end if

  end subroutine fson_get_default_array_2d_double

!------------------------------------------------------------------------

  subroutine fson_get_default_array_2d_logical(self, path, default, &
       val, logfile, log_key)
    !! Gets 2-D logical array with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscBool, intent(in) :: default(:,:)
    PetscBool, allocatable, intent(out) :: val(:,:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(:), allocatable :: logstr, str
    character(max_log_key_length) :: key

    call fson_get(self, path, p)
    if (assoc_non_null(p)) then
       call fson_get(p, ".", val)
    else
       val = default
       if (present(logfile)) then
          if (logfile%active) then
             if (present(log_key)) then
                key = log_key
             else
                key = path
             end if
             if (size(val) > 0) then
                write(logstr, '(L)') val(1,1)
                str = '[[' // logstr // ',...]]'
             else
                str = '[]'
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  str_key = key, str_value = str)
          end if
       end if
    end if

  end subroutine fson_get_default_array_2d_logical

!------------------------------------------------------------------------
! fson_get_mpi routines
!------------------------------------------------------------------------

  subroutine fson_get_mpi_fson_value(self, path, val)
    !! Gets fson_value on MPI input rank, and null on all other ranks.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    type(fson_value), pointer, intent(out) :: val

    if (mpi%rank == mpi%input_rank) then
       call fson_get(self, path, val)
    else
       val => NULL()
    end if

  end subroutine fson_get_mpi_fson_value

!------------------------------------------------------------------------

  subroutine fson_get_mpi_integer(self, path, default, val, logfile, &
       log_key)
    !! Gets integer value on all ranks, with optional default.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscInt, intent(in), optional :: default
    PetscInt, intent(out) :: val
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    PetscInt :: ierr

    if (mpi%rank == mpi%input_rank) then
       if (present(default)) then
          call fson_get_default(self, path, default, val, logfile, &
               log_key)
       else
          call fson_get(self, path, val)
       end if
    end if
    call MPI_bcast(val, 1, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)

  end subroutine fson_get_mpi_integer

!------------------------------------------------------------------------

  subroutine fson_get_mpi_real(self, path, default, val, logfile, &
       log_key)
    !! Gets real value on all ranks, with optional default.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    real, intent(in), optional :: default
    real, intent(out) :: val
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    PetscInt :: ierr

    if (mpi%rank == mpi%input_rank) then
       if (present(default)) then
          call fson_get_default(self, path, default, val, logfile, &
               log_key)
       else
          call fson_get(self, path, val)
       end if
    end if
    call MPI_bcast(val, 1, MPI_REAL, mpi%input_rank, mpi%comm, ierr)

  end subroutine fson_get_mpi_real

!------------------------------------------------------------------------

  subroutine fson_get_mpi_double(self, path, default, val, logfile, &
       log_key)
    !! Gets double precision value on all ranks, with optional default.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscReal, intent(in), optional :: default
    PetscReal, intent(out) :: val
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    PetscInt :: ierr

    if (mpi%rank == mpi%input_rank) then
       if (present(default)) then
          call fson_get_default(self, path, default, val, logfile, &
               log_key)
       else
          call fson_get(self, path, val)
       end if
    end if
    call MPI_bcast(val, 1, MPI_DOUBLE_PRECISION, mpi%input_rank, &
         mpi%comm, ierr)

  end subroutine fson_get_mpi_double

!------------------------------------------------------------------------

  subroutine fson_get_mpi_logical(self, path, default, val, logfile, &
       log_key)
    !! Gets logical value on all ranks, with optional default.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscBool, intent(in), optional :: default
    PetscBool, intent(out) :: val
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    PetscInt :: ierr

    if (mpi%rank == mpi%input_rank) then
       if (present(default)) then
          call fson_get_default(self, path, default, val, logfile, &
               log_key)
       else
          call fson_get(self, path, val)
       end if
    end if
    call MPI_bcast(val, 1, MPI_LOGICAL, mpi%input_rank, mpi%comm, ierr)

  end subroutine fson_get_mpi_logical

!------------------------------------------------------------------------

  subroutine fson_get_mpi_character(self, path, default, val, logfile, &
       log_key)
    !! Gets character value on all ranks, with optional default.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    character(len=*), intent(in), optional :: default
    character(len=*), intent(out) :: val
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    PetscInt :: ierr, count

    if (mpi%rank == mpi%input_rank) then
       if (present(default)) then
          call fson_get_default(self, path, default, val, logfile, &
               log_key)
       else
          call fson_get(self, path, val)
       end if
    end if
    count = len(val)
    call MPI_bcast(val, count, MPI_CHARACTER, mpi%input_rank, &
         mpi%comm, ierr)

  end subroutine fson_get_mpi_character

!------------------------------------------------------------------------

  subroutine fson_get_mpi_array_1d_integer(self, path, default, val, &
       logfile, log_key)
    !! Gets 1-D integer array on all ranks, with optional default.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscInt, intent(in), optional :: default(:)
    PetscInt, allocatable, intent(out) :: val(:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    PetscInt :: ierr, count

    if (mpi%rank == mpi%input_rank) then
       if (present(default)) then
          call fson_get_default(self, path, default, val, logfile, &
               log_key)
       else
          call fson_get(self, path, val)
       end if
       count = size(val)
    end if
    call MPI_bcast(count, 1, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)
    if (mpi%rank /= mpi%input_rank) then
       allocate(val(count))
    end if
    call MPI_bcast(val, count, MPI_INTEGER, mpi%input_rank, mpi%comm, &
         ierr)

  end subroutine fson_get_mpi_array_1d_integer

!------------------------------------------------------------------------

  subroutine fson_get_mpi_array_1d_real(self, path, default, val, &
       logfile, log_key)
    !! Gets 1-D real array on all ranks, with optional default.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    real, intent(in), optional :: default(:)
    real, allocatable, intent(out) :: val(:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    PetscInt :: ierr, count

    if (mpi%rank == mpi%input_rank) then
       if (present(default)) then
          call fson_get_default(self, path, default, val, logfile, &
               log_key)
       else
          call fson_get(self, path, val)
       end if
       count = size(val)
    end if
    call MPI_bcast(count, 1, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)
    if (mpi%rank /= mpi%input_rank) then
       allocate(val(count))
    end if
    call MPI_bcast(val, count, MPI_REAL, mpi%input_rank, mpi%comm, ierr)

  end subroutine fson_get_mpi_array_1d_real

!------------------------------------------------------------------------

  subroutine fson_get_mpi_array_1d_double(self, path, default, val, &
       logfile, log_key)
    !! Gets 1-D double precision array on all ranks, with optional default.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscReal, intent(in), optional :: default(:)
    PetscReal, allocatable, intent(out) :: val(:)
    type(logfile_type), intent(in out), optional :: logfile 
    character(len=*), intent(in), optional :: log_key
   ! Locals:
    PetscInt :: ierr, count

    if (mpi%rank == mpi%input_rank) then
       if (present(default)) then
          call fson_get_default(self, path, default, val, &
               logfile, log_key)
       else
          call fson_get(self, path, val)
       end if
       count = size(val)
    end if
    call MPI_bcast(count, 1, MPI_INTEGER, mpi%input_rank, mpi%comm, &
         ierr)
    if (mpi%rank /= mpi%input_rank) then
       allocate(val(count))
    end if
    call MPI_bcast(val, count, MPI_DOUBLE_PRECISION, mpi%input_rank, &
         mpi%comm, ierr)

  end subroutine fson_get_mpi_array_1d_double

!------------------------------------------------------------------------

  subroutine fson_get_mpi_array_1d_logical(self, path, default, val, &
       logfile, log_key)
    !! Gets 1-D logical array on all ranks, with optional default.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscBool, intent(in), optional :: default(:)
    PetscBool, allocatable, intent(out) :: val(:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    PetscInt :: ierr, count

    if (mpi%rank == mpi%input_rank) then
       if (present(default)) then
          call fson_get_default(self, path, default, val, logfile, &
               log_key)
       else
          call fson_get(self, path, val)
       end if
       count = size(val)
    end if
    call MPI_bcast(count, 1, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)
    if (mpi%rank /= mpi%input_rank) then
       allocate(val(count))
    end if
    call MPI_bcast(val, count, MPI_LOGICAL, mpi%input_rank, mpi%comm, &
         ierr)

  end subroutine fson_get_mpi_array_1d_logical

!------------------------------------------------------------------------

  subroutine fson_get_mpi_array_2d_integer(self, path, default, val, &
       logfile, log_key)
    !! Gets 2-D integer array on all ranks, with optional default.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscInt, intent(in), optional :: default(:,:)
    PetscInt, allocatable, intent(out) :: val(:,:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    PetscInt :: ierr, count(2), total_count

    if (mpi%rank == mpi%input_rank) then
       if (present(default)) then
          call fson_get_default(self, path, default, val, logfile, &
               log_key)
       else
          call fson_get(self, path, val)
       end if
       count = shape(val)
    end if
    call MPI_bcast(count, 2, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)
    if (mpi%rank /= mpi%input_rank) then
       allocate(val(count(1), count(2)))
    end if
    total_count = count(1) * count(2)
    call MPI_bcast(val, total_count, MPI_INTEGER, mpi%input_rank, &
         mpi%comm, ierr)

  end subroutine fson_get_mpi_array_2d_integer

!------------------------------------------------------------------------

  subroutine fson_get_mpi_array_2d_real(self, path, default, val, &
       logfile, log_key)
    !! Gets 2-D real array on all ranks, with optional default.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    real, intent(in), optional :: default(:,:)
    real, allocatable, intent(out) :: val(:,:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    PetscInt :: ierr, count(2), total_count

    if (mpi%rank == mpi%input_rank) then
       if (present(default)) then
          call fson_get_default(self, path, default, val, logfile, &
               log_key)
       else
          call fson_get(self, path, val)
       end if
       count = shape(val)
    end if
    call MPI_bcast(count, 2, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)
    if (mpi%rank /= mpi%input_rank) then
       allocate(val(count(1), count(2)))
    end if
    total_count = count(1) * count(2)
    call MPI_bcast(val, total_count, MPI_REAL, mpi%input_rank, &
         mpi%comm, ierr)

  end subroutine fson_get_mpi_array_2d_real

!------------------------------------------------------------------------

  subroutine fson_get_mpi_array_2d_double(self, path, default, val, &
       logfile, log_key)
    !! Gets 2-D double array on all ranks, with optional default.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscReal, intent(in), optional :: default(:,:)
    PetscReal, allocatable, intent(out) :: val(:,:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    PetscInt :: ierr, count(2), total_count

    if (mpi%rank == mpi%input_rank) then
       if (present(default)) then
          call fson_get_default(self, path, default, val, logfile, &
               log_key)
       else
          call fson_get(self, path, val)
       end if
       count = shape(val)
    end if
    call MPI_bcast(count, 2, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)
    if (mpi%rank /= mpi%input_rank) then
       allocate(val(count(1), count(2)))
    end if
    total_count = count(1) * count(2)
    call MPI_bcast(val, total_count, MPI_DOUBLE_PRECISION, &
         mpi%input_rank, mpi%comm, ierr)

  end subroutine fson_get_mpi_array_2d_double

!------------------------------------------------------------------------

  subroutine fson_get_mpi_array_2d_logical(self, path, default, val, &
       logfile, log_key)
    !! Gets 2-D logical array on all ranks, with optional default.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscBool, intent(in), optional :: default(:,:)
    PetscBool, allocatable, intent(out) :: val(:,:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    integer :: ierr, count(2), total_count

    if (mpi%rank == mpi%input_rank) then
       if (present(default)) then
          call fson_get_default(self, path, default, val, logfile, &
               log_key)
       else
          call fson_get(self, path, val)
       end if
       count = shape(val)
    end if
    call MPI_bcast(count, 2, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)
    if (mpi%rank /= mpi%input_rank) then
       allocate(val(count(1), count(2)))
    end if
    total_count = count(1) * count(2)
    call MPI_bcast(val, total_count, MPI_LOGICAL, mpi%input_rank, &
         mpi%comm, ierr)

  end subroutine fson_get_mpi_array_2d_logical

!------------------------------------------------------------------------

  PetscBool function fson_has_mpi(self, path) result(has)
    !! Returns .true. on all ranks if fson object has the specified
    !! path, and .false. otherwise.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    ! Locals:
    type(fson_value), pointer :: p
    PetscInt :: ierr

    if (mpi%rank == mpi%input_rank) then
       call fson_get(self, path, p)
       has = (associated(p))
    end if

    call MPI_bcast(has, 1, MPI_LOGICAL, mpi%input_rank, mpi%comm, ierr)

  end function fson_has_mpi

!------------------------------------------------------------------------

  PetscInt function fson_type_mpi(self, path) result(t)
    !! Returns value type on all ranks of fson object with the specified
    !! path (TYPE_NULL if the path does not exist).

    use fson_value_m, only : TYPE_NULL

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    ! Locals:
    type(fson_value), pointer :: p
    PetscInt :: ierr

    if (mpi%rank == mpi%input_rank) then
       call fson_get(self, path, p)
       if (associated(p)) then
          t = p%value_type
       else
          t = TYPE_NULL
       end if
    end if

    call MPI_bcast(t, 1, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)

  end function fson_type_mpi

!------------------------------------------------------------------------

  PetscInt function fson_value_count_mpi(self, path) result(count)
    !! Returns value count on all ranks of fson object with the specified
    !! path (returns zero if the path does not exist).

    use fson_value_m, only : fson_value_count

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    ! Locals:
    type(fson_value), pointer :: p
    PetscInt :: ierr

    if (mpi%rank == mpi%input_rank) then
       call fson_get(self, path, p)
       if (assoc_non_null(p)) then
          count = fson_value_count(p)
       else
          count = 0
       end if
    end if

    call MPI_bcast(count, 1, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)

  end function fson_value_count_mpi

!------------------------------------------------------------------------

  function fson_value_get_mpi(self, i) result(p)
    !! Returns value i of fson object on MPI input rank, and null
    !! on all other ranks.

    use fson_value_m, only : fson_value_get

    type(fson_value), pointer, intent(in) :: self
    PetscInt, intent(in) :: i
    type(fson_value), pointer :: p

    if (mpi%rank == mpi%input_rank) then
       p => fson_value_get(self, i)
    else
       p => NULL()
    end if

  end function fson_value_get_mpi

!------------------------------------------------------------------------

  function fson_parse_mpi(file, unit, str) result(p)
    !! Returns parsed fson object corresponding to specified file
    !! on MPI input rank, and null on all other ranks.

    type(fson_value), pointer :: p
    character(len = *), intent(in), optional :: file
    PetscInt, intent(in out), optional :: unit
    character(len = *), intent(in), optional :: str

    if (mpi%rank == mpi%input_rank) then
       p => fson_parse(file, unit, str)
    else
       nullify(p)
    end if

  end function fson_parse_mpi

!------------------------------------------------------------------------

  subroutine fson_destroy_mpi(self)
    !! Destroys fson object on MPI input rank.
    type(fson_value), pointer, intent(in out) :: self

    if (mpi%rank == mpi%input_rank) then
       call fson_destroy(self)
       nullify(self)
    end if

  end subroutine fson_destroy_mpi

!------------------------------------------------------------------------

end module fson_mpi_module
