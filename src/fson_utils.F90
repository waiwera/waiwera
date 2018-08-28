!   Copyright 2018 University of Auckland.

!   This file is part of Waiwera.

!   Waiwera is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.

!   Waiwera is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.

!   You should have received a copy of the GNU Lesser General Public License
!   along with Waiwera.  If not, see <http://www.gnu.org/licenses/>.

module fson_utils_module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use logfile_module
  use fson

  implicit none
  private

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
     module procedure fson_get_default_array_1d_character
     module procedure fson_get_default_array_2d_integer
     module procedure fson_get_default_array_2d_real
     module procedure fson_get_default_array_2d_double
     module procedure fson_get_default_array_2d_logical
  end interface fson_get_default

  public :: assoc_non_null, fson_array_rank, fson_get_default

contains

!------------------------------------------------------------------------

  PetscBool function assoc_non_null(self)
    !! Returns .true. if fson_value self is associated and has value
    !! type non-null.

    use fson_value_m, only : TYPE_NULL

    type(fson_value), pointer, intent(in) :: self

    assoc_non_null = PETSC_FALSE
    if (associated(self)) then
       if (self%value_type /= TYPE_NULL) then
          assoc_non_null = PETSC_TRUE
       end if
    end if

  end function assoc_non_null

!------------------------------------------------------------------------

  PetscInt function fson_array_rank(self, path) result(r)
    !! Returns array rank of object: -1 if the object does not exist,
    !! zero if it is a scalar or zero-length array, 1 or 2 if it is a
    !! non-empty array. (Arrays of rank higher than 2 are not
    !! detected.)

    use fson_value_m, only : TYPE_ARRAY, fson_value_count, fson_value_get

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    ! Locals:
    type(fson_value), pointer :: p, p1
    PetscInt :: t, count

    call fson_get(self, path, p)
    if (associated(p)) then
       t = p%value_type
       if (t == TYPE_ARRAY) then
          count = fson_value_count(p)
          if (count == 0) then
             r = 0
          else
             p1 => fson_value_get(p, 1)
             if (p1%value_type == TYPE_ARRAY) then
                r = 2
             else
                r = 1
             end if
          end if
       else
          r = 0
       end if
    else
       r = -1
    end if

  end function fson_array_rank

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
                allocate(character(logfile%max_num_length) :: intstr)
                write(intstr, logfile%int_format) val(1)
                str = '[' // trim(intstr) // ',...]'
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
                allocate(character(logfile%max_num_length) :: realstr)
                write(realstr, logfile%real_format) val(1)
                str  = '[' // trim(realstr) // ',...]'
             else
                str  = '[]'
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
                allocate(character(logfile%max_num_length) :: realstr)
                write(realstr, logfile%real_format) val(1)
                str = '[' // trim(realstr) // ',...]'
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
    character :: logstr
    character(:), allocatable :: str
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

  subroutine fson_get_default_array_1d_character(self, path, default, &
       string_length, val, logfile, log_key)
    !! Gets 1-D character array with default if not present.

    type(fson_value), pointer, intent(in) :: self
    character(len=*), intent(in) :: path
    PetscInt, intent(in) :: string_length
    character(string_length), intent(in) :: default(:)
    character(string_length), allocatable, intent(out) :: val(:)
    type(logfile_type), intent(in out), optional :: logfile
    character(len=*), intent(in), optional :: log_key
    ! Locals:
    type(fson_value), pointer :: p
    character(string_length) :: logstr
    character(:), allocatable :: str
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
                write(logstr, '(a)') val(1)
                str = '[' // trim(logstr) // ',...]'
             else
                str = '[]'
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  str_key = key, str_value = str)
          end if
       end if
    end if

  end subroutine fson_get_default_array_1d_character

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
                allocate(character(logfile%max_num_length) :: intstr)
                write(intstr, logfile%int_format) val(1,1)
                str = '[[' // trim(intstr) // ',...]]'
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
                allocate(character(logfile%max_num_length) :: realstr)
                write(realstr, logfile%real_format) val(1,1)
                str  = '[[' // trim(realstr) // ',...]]'
             else
                str  = '[]'
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
                allocate(character(logfile%max_num_length) :: realstr)
                write(realstr, logfile%real_format) val(1,1)
                str = '[[' // trim(realstr) // ',...]]'
             else
                str  = '[]'
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
    character :: logstr
    character(:), allocatable :: str
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
                str ='[[' // logstr // ',...]]'
             else
                str  = '[]'
             end if
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  str_key = key, str_value = str)
          end if
       end if
    end if

  end subroutine fson_get_default_array_2d_logical

!------------------------------------------------------------------------


end module fson_utils_module
