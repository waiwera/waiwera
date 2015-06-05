module eos_setup_module
  !! Module for setting up EOS from input data.

  use eos_w_module
  use eos_module
  use fson
  use fson_mpi_module
  use thermodynamics_module

  implicit none
  private

#include <petsc-finclude/petsc.h90>

  public :: setup_eos

  contains

!------------------------------------------------------------------------

  subroutine setup_eos(json, thermo, eos)
    !! Reads equation of state from JSON input file.
    !! If not present, a default value is assigned.

    use utils_module, only : str_to_lower

    type(fson_value), pointer, intent(in) :: json
    class(thermodynamics_type), intent(in) :: thermo
    class(eos_type), pointer, intent(in out) :: eos
    ! Locals:
    character(max_eos_name_length), parameter :: &
         default_eos_name = "w"
    character(max_eos_name_length) :: eos_name

    call fson_get_mpi(json, "eos.name", default_eos_name, eos_name)
    eos_name = str_to_lower(eos_name)

    select case (eos_name)
    case ("ew")
       allocate(eos_w_type :: eos)  ! change to eos_ew when it's ready
    case default
       allocate(eos_w_type :: eos)
    end select

    call eos%init(json, thermo)

  end subroutine setup_eos

!------------------------------------------------------------------------

end module eos_setup_module
