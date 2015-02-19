module mpi_module
  !! MPI variables.

  implicit none
  private

#include <petsc-finclude/petscdef.h>

  type, public :: mpi_type
     MPI_Comm, public :: comm !! MPI communicator
     PetscMPIInt, public :: size !! number of processors
     PetscMPIInt, public :: rank !! processor rank
     PetscMPIInt, public :: input_rank = 0 !! Processor rank for input handling
     PetscMPIInt, public :: output_rank = 0 !! Processor rank for output handling
   contains
     private
     procedure, public :: setup => mpi_setup
  end type mpi_type

  type(mpi_type), public :: mpi

contains

!------------------------------------------------------------------------

  subroutine mpi_setup(self, comm)
    !! Sets up MPI object.

    class(mpi_type), intent(in out) :: self
    MPI_Comm, intent(in) :: comm !! MPI communicator
    ! Locals:
    integer :: ierr

    self%comm = comm
    call MPI_COMM_RANK(self%comm, self%rank, ierr)
    call MPI_COMM_SIZE(self%comm, self%size, ierr)
    
  end subroutine mpi_setup

!------------------------------------------------------------------------

end module mpi_module
