module eos_setup_module

  use fson
  use utils_module, only : str_to_lower
  use fson_mpi_module, only: fson_get_mpi
  use logfile_module
  use thermodynamics_module
  use eos_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>
#include <petsc/finclude/petscsys.h>

  public :: setup_eos

contains

!------------------------------------------------------------------------
  
  subroutine setup_eos(json, thermo, eos, logfile)
    !! Reads equation of state from JSON input file.
    !! If not present, a default value is assigned.

    type(fson_value), pointer, intent(in) :: json
    class(thermodynamics_type), intent(in) :: thermo
    class(eos_type), allocatable, intent(in out) :: eos
    type(logfile_type), intent(in out) :: logfile
    ! Locals:
    character(max_eos_name_length), parameter :: &
         default_eos_name = "we"
    character(max_eos_name_length) :: eos_name

    call fson_get_mpi(json, "eos.name", default_eos_name, &
         eos_name, logfile)
    eos_name = str_to_lower(eos_name)

    select case (eos_name)
    case ("w")
       allocate(eos_w_type :: eos)
    case ("we")
       allocate(eos_we_type :: eos)
    case default
       allocate(eos_we_type :: eos)
    end select

    call eos%init(json, thermo, logfile)

  end subroutine setup_eos

!------------------------------------------------------------------------

end module eos_setup_module
