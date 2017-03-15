module rock_parameters_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  type rock_type
     PetscReal, pointer :: density

   contains
     private

  end type rock_type

  PetscInt, parameter :: num_rock_variables = 3
  PetscInt, parameter :: max_rock_variable_name_length = 32

end module rock_parameters_module
