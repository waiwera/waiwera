module unit_test_utils_module
  !! Unit testing utilities

#include <petsc/finclude/petscsys.h>

  use petscsys
  use zofu

  implicit none
  private

  public :: transition_compare

contains

!------------------------------------------------------------------------

  subroutine transition_compare(test, expected_primary, expected_region, &
       expected_transition, primary, fluid, transition, err, message)

    use fluid_module

    ! Runs assertions to test EOS transition

    class(unit_test_type), intent(in out) :: test
    PetscReal, intent(in) :: expected_primary(:), primary(:)
    PetscInt, intent(in) :: expected_region
    type(fluid_type), intent(in) :: fluid
    PetscBool, intent(in) :: expected_transition, transition
    PetscErrorCode, intent(in) :: err
    character(60), intent(in) :: message

    call test%assert(expected_primary, primary, &
         trim(message) // "primary")
    call test%assert(expected_region, nint(fluid%region), &
         trim(message) // " region")
    call test%assert(expected_transition, transition, &
         trim(message) // " transition")
    call test%assert(0, err, trim(message) // " error")

  end subroutine transition_compare

!------------------------------------------------------------------------

end module unit_test_utils_module
