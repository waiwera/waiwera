module version_module
  !! Module for software version information.

  implicit none
  private

#include <petsc/finclude/petscdef.h>
  
  PetscInt, parameter :: max_version_string_length = 64
  character(len = max_version_string_length), public :: &
       supermodel_version = "0.0.1"

end module version_module
