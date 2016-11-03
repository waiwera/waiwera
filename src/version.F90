module version_module
  !! Module for software version information.

  implicit none
  private

#include <petsc/finclude/petscdef.h>
  
  PetscInt, parameter :: max_version_string_length = 64
  character(len = max_version_string_length), public :: &
       waiwera_version = "0.1.0"

end module version_module
