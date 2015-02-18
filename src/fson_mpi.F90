module fson_mpi_module
  !! Subroutines for JSON input under MPI, with input on one rank
  !! broadcast to other ranks, and with the ability to specify
  !! default values if not present in the JSON input.

  use mpi_module
  use fson

  private

#include <petsc-finclude/petscsys.h>

  public :: fson_get_default, fson_get_mpi

  interface fson_get_default
     module procedure fson_get_default_integer
     module procedure fson_get_default_real
     module procedure fson_get_default_character
  end interface fson_get_default

  interface fson_get_mpi
     module procedure fson_get_mpi_integer
     module procedure fson_get_mpi_real
     module procedure fson_get_mpi_character
  end interface fson_get_mpi

  contains
    
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

    subroutine fson_get_mpi_integer(self, path, default, val)
      !! Gets integer value on all ranks, with default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        integer, intent(in) :: default
        integer, intent(out) :: val
        ! Locals:
        integer :: ierr

        if (rank == input_rank) then
           call fson_get_default(self, path, default, val)
        end if
        call MPI_bcast(val, 1, MPI_INTEGER, input_rank, comm, ierr)

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

        if (rank == input_rank) then
           call fson_get_default(self, path, default, val)
        end if
        call MPI_bcast(val, 1, MPI_REAL, input_rank, comm, ierr)

      end subroutine fson_get_mpi_real

!------------------------------------------------------------------------

    subroutine fson_get_mpi_character(self, path, default, val)
      !! Gets character value on all ranks, with default.

        type(fson_value), pointer, intent(in) :: self
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: default
        character(len=*), intent(out) :: val
        ! Locals:
        integer :: ierr, count

        if (rank == input_rank) then
           call fson_get_default(self, path, default, val)
        end if
        count = len(val)
        call MPI_bcast(val, count, MPI_CHARACTER, input_rank, comm, ierr)

      end subroutine fson_get_mpi_character

!------------------------------------------------------------------------

end module fson_mpi_module
