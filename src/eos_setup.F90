module eos_setup_module
  !! Module for setting up EOS from input data.

  use eos_w_module
  use eos_module
  use fson
  use fson_mpi_module
  use thermodynamics_module

  implicit none
  private

  integer, parameter :: max_eos_ID_length = 12
  character(max_eos_ID_length), parameter :: &
       default_eos_ID = "W"

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
    character(max_eos_ID_length) :: eos_ID

    call fson_get_mpi(json, "eos", default_eos_ID, eos_ID)
    eos_ID = str_to_lower(eos_ID)

    select case (eos_ID)
    case ("ew")
       allocate(eos_w_type :: eos)  ! change to eos_ew when it's ready
    case default
       allocate(eos_w_type :: eos)
    end select

    call eos%init(thermo)

  end subroutine setup_eos

!------------------------------------------------------------------------

end module eos_setup_module
