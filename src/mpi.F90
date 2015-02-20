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
     procedure, public :: init => mpi_initialize
  end type mpi_type

  type(mpi_type), public :: mpi

contains

!------------------------------------------------------------------------

  subroutine mpi_initialize(self, comm)
    !! Initializes MPI object. (Note this does not call MPI_INIT- that is
    !! already done by PETSC_INITIALIZE.)

    class(mpi_type), intent(in out) :: self
    MPI_Comm, intent(in) :: comm !! MPI communicator
    ! Locals:
    integer :: ierr

    self%comm = comm
    call MPI_COMM_RANK(self%comm, self%rank, ierr)
    call MPI_COMM_SIZE(self%comm, self%size, ierr)
    
  end subroutine mpi_initialize

!------------------------------------------------------------------------

end module mpi_module
