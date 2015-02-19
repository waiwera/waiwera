module fson_mpi_module
  !! Subroutines for JSON input under MPI, with input on one rank
  !! broadcast to other ranks, and with the ability to specify
  !! default values if not present in the JSON input.

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
  end interface fson_get_default

  interface fson_get_mpi
     module procedure fson_get_mpi_integer
     module procedure fson_get_mpi_real
     module procedure fson_get_mpi_double
     module procedure fson_get_mpi_logical
     module procedure fson_get_mpi_character
  end interface fson_get_mpi

  public :: fson_get_default, fson_get_mpi

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
! fson_get_mpi routines
!------------------------------------------------------------------------

    subroutine fson_get_mpi_integer(self, path, default, val)
      !! Gets integer value on all ranks, with default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        integer, intent(in) :: default
        integer, intent(out) :: val
        ! Locals:
        integer :: ierr

        if (mpi%rank == mpi%input_rank) then
           call fson_get_default(self, path, default, val)
        end if
        call MPI_bcast(val, 1, MPI_INTEGER, mpi%input_rank, mpi%comm, ierr)

    end subroutine fson_get_mpi_integer

!------------------------------------------------------------------------

    subroutine fson_get_mpi_real(self, path, default, val)
      !! Gets real value on all ranks, with default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        real, intent(in) :: default
        real, intent(out) :: val
        ! Locals:
        integer :: ierr

        if (mpi%rank == mpi%input_rank) then
           call fson_get_default(self, path, default, val)
        end if
        call MPI_bcast(val, 1, MPI_REAL, mpi%input_rank, mpi%comm, ierr)

      end subroutine fson_get_mpi_real

!------------------------------------------------------------------------

    subroutine fson_get_mpi_double(self, path, default, val)
      !! Gets double precision value on all ranks, with default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        double precision, intent(in) :: default
        double precision, intent(out) :: val
        ! Locals:
        integer :: ierr

        if (mpi%rank == mpi%input_rank) then
           call fson_get_default(self, path, default, val)
        end if
        call MPI_bcast(val, 1, MPI_DOUBLE_PRECISION, mpi%input_rank, &
             mpi%comm, ierr)

      end subroutine fson_get_mpi_double

!------------------------------------------------------------------------

    subroutine fson_get_mpi_logical(self, path, default, val)
      !! Gets logical value on all ranks, with default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        logical, intent(in) :: default
        logical, intent(out) :: val
        ! Locals:
        integer :: ierr

        if (mpi%rank == mpi%input_rank) then
           call fson_get_default(self, path, default, val)
        end if
        call MPI_bcast(val, 1, MPI_LOGICAL, mpi%input_rank, mpi%comm, ierr)

      end subroutine fson_get_mpi_logical

!------------------------------------------------------------------------

    subroutine fson_get_mpi_character(self, path, default, val)
      !! Gets character value on all ranks, with default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: default
        character(len=*), intent(out) :: val
        ! Locals:
        integer :: ierr, count

        if (mpi%rank == mpi%input_rank) then
           call fson_get_default(self, path, default, val)
        end if
        count = len(val)
        call MPI_bcast(val, count, MPI_CHARACTER, mpi%input_rank, mpi%comm, ierr)

      end subroutine fson_get_mpi_character

!------------------------------------------------------------------------

end module fson_mpi_module
