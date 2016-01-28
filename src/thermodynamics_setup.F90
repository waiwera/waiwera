module thermodynamics_setup_module
  !! Module for setting up thermodynamics from input.

  use thermodynamics_module
  use IAPWS_module
  use IFC67_module
  use fson
  use fson_mpi_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>

  PetscInt, parameter :: max_thermo_ID_length = 8
  character(max_thermo_ID_length), parameter :: &
       default_thermo_ID = "IAPWS"

  public :: setup_thermodynamics

contains

!------------------------------------------------------------------------

  subroutine setup_thermodynamics(json, thermo, logfile)
    !! Reads thermodynamic formulation from JSON input file.
    !! If not present, a default value is assigned.

    use utils_module, only : str_to_lower
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    class(thermodynamics_type), allocatable, intent(in out) :: thermo
    type(logfile_type), intent(in out) :: logfile
    ! Locals:
    character(max_thermo_ID_length) :: thermo_ID

    call fson_get_mpi(json, "thermodynamics", default_thermo_ID, &
         thermo_ID, logfile)
    thermo_ID = str_to_lower(thermo_ID)

    select case (thermo_ID)
    case ("ifc67")
       allocate(IFC67_type :: thermo)
    case default
       allocate(IAPWS_type :: thermo)
    end select

    call thermo%init()

  end subroutine setup_thermodynamics

!------------------------------------------------------------------------

end module thermodynamics_setup_module
