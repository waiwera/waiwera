module fson_mpi_module
  !! Subroutines for JSON input under MPI, with input on one rank
  !! broadcast to other ranks, and with the ability to specify
  !! default values if not present in the JSON input.

  use kinds_module
  use mpi_module
  use fson

  implicit none
  private

#include <petsc-finclude/petscsys.h>

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

  contains
    
!------------------------------------------------------------------------
! fson_get_default routines
!------------------------------------------------------------------------

    subroutine fson_get_default_integer(self, path, default, val)
      !! Gets integer value with default if not present.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        integer, intent(in) :: default
        integer, intent(out) :: val
        ! Locals:
        type(fson_value), pointer :: p

        call fson_get(self, path, p)
        if (associated(p)) then
           call fson_get(p, ".", val)
        else
           val = default
        end if

    end subroutine fson_get_default_integer

!------------------------------------------------------------------------

    subroutine fson_get_default_real(self, path, default, val)
      !! Gets real value with default if not present.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        real, intent(in) :: default
        real, intent(out) :: val
        ! Locals:
        type(fson_value), pointer :: p

        call fson_get(self, path, p)
        if (associated(p)) then
           call fson_get(p, ".", val)
        else
           val = default
        end if

      end subroutine fson_get_default_real

!------------------------------------------------------------------------

    subroutine fson_get_default_double(self, path, default, val)
      !! Gets double value with default if not present.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        double precision, intent(in) :: default
        double precision, intent(out) :: val
        ! Locals:
        type(fson_value), pointer :: p

        call fson_get(self, path, p)
        if (associated(p)) then
           call fson_get(p, ".", val)
        else
           val = default
        end if

      end subroutine fson_get_default_double

!------------------------------------------------------------------------

    subroutine fson_get_default_logical(self, path, default, val)
      !! Gets logical value with default if not present.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        logical, intent(in) :: default
        logical, intent(out) :: val
        ! Locals:
        type(fson_value), pointer :: p

        call fson_get(self, path, p)
        if (associated(p)) then
           call fson_get(p, ".", val)
        else
           val = default
        end if

      end subroutine fson_get_default_logical

!------------------------------------------------------------------------

    subroutine fson_get_default_character(self, path, default, val)
      !! Gets character value with default if not present.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: default
        character(len=*), intent(out) :: val
        ! Locals:
        type(fson_value), pointer :: p

        call fson_get(self, path, p)
        if (associated(p)) then
           call fson_get(p, ".", val)
        else
           val = default
        end if

      end subroutine fson_get_default_character

!------------------------------------------------------------------------

      subroutine fson_get_default_array_1d_integer(self, path, default, val)
        !! Gets 1-D integer array with default if not present.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        integer, intent(in) :: default(:)
        integer, allocatable, intent(out) :: val(:)
        ! Locals:
        type(fson_value), pointer :: p

        call fson_get(self, path, p)
        if (associated(p)) then
           call fson_get(p, ".", val)
        else
           val = default
        end if

      end subroutine fson_get_default_array_1d_integer

!------------------------------------------------------------------------

      subroutine fson_get_default_array_1d_real(self, path, default, val)
        !! Gets 1-D real array with default if not present.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        real, intent(in) :: default(:)
        real, allocatable, intent(out) :: val(:)
        ! Locals:
        type(fson_value), pointer :: p

        call fson_get(self, path, p)
        if (associated(p)) then
           call fson_get(p, ".", val)
        else
           val = default
        end if

      end subroutine fson_get_default_array_1d_real

!------------------------------------------------------------------------

      subroutine fson_get_default_array_1d_double(self, path, default, val)
        !! Gets 1-D double precision array with default if not present.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: default(:)
        real(dp), allocatable, intent(out) :: val(:)
        ! Locals:
        type(fson_value), pointer :: p

        call fson_get(self, path, p)
        if (associated(p)) then
           call fson_get(p, ".", val)
        else
           val = default
        end if

      end subroutine fson_get_default_array_1d_double

!------------------------------------------------------------------------

      subroutine fson_get_default_array_1d_logical(self, path, default, val)
        !! Gets 1-D logical array with default if not present.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        logical, intent(in) :: default(:)
        logical, allocatable, intent(out) :: val(:)
        ! Locals:
        type(fson_value), pointer :: p

        call fson_get(self, path, p)
        if (associated(p)) then
           call fson_get(p, ".", val)
        else
           val = default
        end if

      end subroutine fson_get_default_array_1d_logical

!------------------------------------------------------------------------

      subroutine fson_get_default_array_2d_integer(self, path, default, val)
        !! Gets 2-D integer array with default if not present.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        integer, intent(in) :: default(:,:)
        integer, allocatable, intent(out) :: val(:,:)
        ! Locals:
        type(fson_value), pointer :: p

        call fson_get(self, path, p)
        if (associated(p)) then
           call fson_get(p, ".", val)
        else
           val = default
        end if

      end subroutine fson_get_default_array_2d_integer

!------------------------------------------------------------------------

      subroutine fson_get_default_array_2d_real(self, path, default, val)
        !! Gets 2-D real array with default if not present.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        real, intent(in) :: default(:,:)
        real, allocatable, intent(out) :: val(:,:)
        ! Locals:
        type(fson_value), pointer :: p

        call fson_get(self, path, p)
        if (associated(p)) then
           call fson_get(p, ".", val)
        else
           val = default
        end if

      end subroutine fson_get_default_array_2d_real

!------------------------------------------------------------------------

      subroutine fson_get_default_array_2d_double(self, path, default, val)
        !! Gets 2-D double array with default if not present.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: default(:,:)
        real(dp), allocatable, intent(out) :: val(:,:)
        ! Locals:
        type(fson_value), pointer :: p

        call fson_get(self, path, p)
        if (associated(p)) then
           call fson_get(p, ".", val)
        else
           val = default
        end if

      end subroutine fson_get_default_array_2d_double

!------------------------------------------------------------------------

      subroutine fson_get_default_array_2d_logical(self, path, default, val)
        !! Gets 2-D logical array with default if not present.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        logical, intent(in) :: default(:,:)
        logical, allocatable, intent(out) :: val(:,:)
        ! Locals:
        type(fson_value), pointer :: p

        call fson_get(self, path, p)
        if (associated(p)) then
           call fson_get(p, ".", val)
        else
           val = default
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

    subroutine fson_get_mpi_integer(self, path, default, val)
      !! Gets integer value on all ranks, with optional default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        integer, intent(in), optional :: default
        integer, intent(out) :: val
        ! Locals:
        integer :: ierr

        if (mpi%rank == mpi%input_rank) then
           if (present(default)) then
              call fson_get_default(self, path, default, val)
           else
              call fson_get(self, path, val)
           end if
        end if
        call MPI_bcast(val, 1, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)

    end subroutine fson_get_mpi_integer

!------------------------------------------------------------------------

    subroutine fson_get_mpi_real(self, path, default, val)
      !! Gets real value on all ranks, with optional default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        real, intent(in), optional :: default
        real, intent(out) :: val
        ! Locals:
        integer :: ierr

        if (mpi%rank == mpi%input_rank) then
           if (present(default)) then
              call fson_get_default(self, path, default, val)
           else
              call fson_get(self, path, val)
           end if
        end if
        call MPI_bcast(val, 1, MPI_REAL, mpi%input_rank, mpi%comm, ierr)

      end subroutine fson_get_mpi_real

!------------------------------------------------------------------------

    subroutine fson_get_mpi_double(self, path, default, val)
      !! Gets double precision value on all ranks, with optional default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        double precision, intent(in), optional :: default
        double precision, intent(out) :: val
        ! Locals:
        integer :: ierr

        if (mpi%rank == mpi%input_rank) then
           if (present(default)) then
              call fson_get_default(self, path, default, val)
           else
              call fson_get(self, path, val)
           end if
        end if
        call MPI_bcast(val, 1, MPI_DOUBLE_PRECISION, mpi%input_rank, &
             mpi%comm, ierr)

      end subroutine fson_get_mpi_double

!------------------------------------------------------------------------

    subroutine fson_get_mpi_logical(self, path, default, val)
      !! Gets logical value on all ranks, with optional default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        logical, intent(in), optional :: default
        logical, intent(out) :: val
        ! Locals:
        integer :: ierr

        if (mpi%rank == mpi%input_rank) then
           if (present(default)) then
              call fson_get_default(self, path, default, val)
           else
              call fson_get(self, path, val)
           end if
        end if
        call MPI_bcast(val, 1, MPI_LOGICAL, mpi%input_rank, mpi%comm, ierr)

      end subroutine fson_get_mpi_logical

!------------------------------------------------------------------------

    subroutine fson_get_mpi_character(self, path, default, val)
      !! Gets character value on all ranks, with optional default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        character(len=*), intent(in), optional :: default
        character(len=*), intent(out) :: val
        ! Locals:
        integer :: ierr, count

        if (mpi%rank == mpi%input_rank) then
           if (present(default)) then
              call fson_get_default(self, path, default, val)
           else
              call fson_get(self, path, val)
           end if
        end if
        count = len(val)
        call MPI_bcast(val, count, MPI_CHARACTER, mpi%input_rank, mpi%comm, ierr)

      end subroutine fson_get_mpi_character

!------------------------------------------------------------------------

    subroutine fson_get_mpi_array_1d_integer(self, path, default, val)
      !! Gets 1-D integer array on all ranks, with optional default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        integer, intent(in), optional :: default(:)
        integer, allocatable, intent(out) :: val(:)
        ! Locals:
        integer :: ierr, count

        if (mpi%rank == mpi%input_rank) then
           if (present(default)) then
              call fson_get_default(self, path, default, val)
           else
              call fson_get(self, path, val)
           end if
           count = size(val)
        end if
        call MPI_bcast(count, 1, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)
        if (mpi%rank /= mpi%input_rank) then
           allocate(val(count))
        end if
        call MPI_bcast(val, count, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)

      end subroutine fson_get_mpi_array_1d_integer

!------------------------------------------------------------------------

    subroutine fson_get_mpi_array_1d_real(self, path, default, val)
      !! Gets 1-D real array on all ranks, with optional default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        real, intent(in), optional :: default(:)
        real, allocatable, intent(out) :: val(:)
        ! Locals:
        integer :: ierr, count

        if (mpi%rank == mpi%input_rank) then
           if (present(default)) then
              call fson_get_default(self, path, default, val)
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

    subroutine fson_get_mpi_array_1d_double(self, path, default, val)
      !! Gets 1-D double precision array on all ranks, with optional default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        real(dp), intent(in), optional :: default(:)
        real(dp), allocatable, intent(out) :: val(:)
        ! Locals:
        integer :: ierr, count

        if (mpi%rank == mpi%input_rank) then
           if (present(default)) then
              call fson_get_default(self, path, default, val)
           else
              call fson_get(self, path, val)
           end if
           count = size(val)
        end if
        call MPI_bcast(count, 1, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)
        if (mpi%rank /= mpi%input_rank) then
           allocate(val(count))
        end if
        call MPI_bcast(val, count, MPI_DOUBLE_PRECISION, mpi%input_rank, &
             mpi%comm, ierr)

      end subroutine fson_get_mpi_array_1d_double

!------------------------------------------------------------------------

    subroutine fson_get_mpi_array_1d_logical(self, path, default, val)
      !! Gets 1-D logical array on all ranks, with optional default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        logical, intent(in), optional :: default(:)
        logical, allocatable, intent(out) :: val(:)
        ! Locals:
        integer :: ierr, count

        if (mpi%rank == mpi%input_rank) then
           if (present(default)) then
              call fson_get_default(self, path, default, val)
           else
              call fson_get(self, path, val)
           end if
           count = size(val)
        end if
        call MPI_bcast(count, 1, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)
        if (mpi%rank /= mpi%input_rank) then
           allocate(val(count))
        end if
        call MPI_bcast(val, count, MPI_LOGICAL, mpi%input_rank, mpi%comm, ierr)

      end subroutine fson_get_mpi_array_1d_logical

!------------------------------------------------------------------------

    subroutine fson_get_mpi_array_2d_integer(self, path, default, val)
      !! Gets 2-D integer array on all ranks, with optional default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        integer, intent(in), optional :: default(:,:)
        integer, allocatable, intent(out) :: val(:,:)
        ! Locals:
        integer :: ierr, count(2), total_count

        if (mpi%rank == mpi%input_rank) then
           if (present(default)) then
              call fson_get_default(self, path, default, val)
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
        call MPI_bcast(val, total_count, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)

      end subroutine fson_get_mpi_array_2d_integer

!------------------------------------------------------------------------

    subroutine fson_get_mpi_array_2d_real(self, path, default, val)
      !! Gets 2-D real array on all ranks, with optional default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        real, intent(in), optional :: default(:,:)
        real, allocatable, intent(out) :: val(:,:)
        ! Locals:
        integer :: ierr, count(2), total_count

        if (mpi%rank == mpi%input_rank) then
           if (present(default)) then
              call fson_get_default(self, path, default, val)
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
        call MPI_bcast(val, total_count, MPI_REAL, mpi%input_rank, mpi%comm, ierr)

      end subroutine fson_get_mpi_array_2d_real

!------------------------------------------------------------------------

    subroutine fson_get_mpi_array_2d_double(self, path, default, val)
      !! Gets 2-D double array on all ranks, with optional default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        real(dp), intent(in), optional :: default(:,:)
        real(dp), allocatable, intent(out) :: val(:,:)
        ! Locals:
        integer :: ierr, count(2), total_count

        if (mpi%rank == mpi%input_rank) then
           if (present(default)) then
              call fson_get_default(self, path, default, val)
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

    subroutine fson_get_mpi_array_2d_logical(self, path, default, val)
      !! Gets 2-D logical array on all ranks, with optional default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        logical, intent(in), optional :: default(:,:)
        logical, allocatable, intent(out) :: val(:,:)
        ! Locals:
        integer :: ierr, count(2), total_count

        if (mpi%rank == mpi%input_rank) then
           if (present(default)) then
              call fson_get_default(self, path, default, val)
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
        call MPI_bcast(val, total_count, MPI_LOGICAL, mpi%input_rank, mpi%comm, ierr)

      end subroutine fson_get_mpi_array_2d_logical

!------------------------------------------------------------------------

      logical function fson_has_mpi(self, path) result(has)
        !! Returns .true. on all ranks if fson object has the specified
        !! path, and .false. otherwise.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        ! Locals:
        type(fson_value), pointer :: p
        integer :: ierr

        if (mpi%rank == mpi%input_rank) then
           call fson_get(self, path, p)
           has = (associated(p))
        end if

        call MPI_bcast(has, 1, MPI_LOGICAL, mpi%input_rank, mpi%comm, ierr)

      end function fson_has_mpi

!------------------------------------------------------------------------

      integer function fson_type_mpi(self, path) result(t)
        !! Returns value type on all ranks of fson object with the specified
        !! path (TYPE_NULL if the path does not exist).

        use fson_value_m, only : TYPE_NULL

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        ! Locals:
        type(fson_value), pointer :: p
        integer :: ierr

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

      integer function fson_value_count_mpi(self, path) result(count)
        !! Returns value count on all ranks of fson object with the specified
        !! path (returns zero if the path does not exist).

        use fson_value_m, only : fson_value_count

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        ! Locals:
        type(fson_value), pointer :: p
        integer :: ierr

        if (mpi%rank == mpi%input_rank) then
           call fson_get(self, path, p)
           if (associated(p)) then
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
        integer, intent(in) :: i
        type(fson_value), pointer :: p

        if (mpi%rank == mpi%input_rank) then
           p => fson_value_get(self, i)
        else
           p => NULL()
        end if

      end function fson_value_get_mpi

!------------------------------------------------------------------------

end module fson_mpi_module
