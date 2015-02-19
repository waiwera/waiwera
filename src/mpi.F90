module mpi_module
  !! MPI variables.

  implicit none

  private

#include <petsc-finclude/petscdef.h>

  MPI_Comm, public :: comm !! MPI communicator
  PetscMPIInt, public :: rank !! processor rank
  PetscMPIInt, public, parameter :: input_rank = 0 !! Processor rank for input handling
  PetscMPIInt, public, parameter :: output_rank = input_rank !! Processor rank for output handling

end module mpi_module
