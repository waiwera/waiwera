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
  end type mpi_type

  type(mpi_type), public :: mpi

end module mpi_module
