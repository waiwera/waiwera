!   Copyright 2016 University of Auckland.

!   This file is part of Waiwera.

!   Waiwera is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.

!   Waiwera is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.

!   You should have received a copy of the GNU Lesser General Public License
!   along with Waiwera.  If not, see <http://www.gnu.org/licenses/>.

module thermodynamics_setup_module
  !! Module for setting up thermodynamics from input.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use thermodynamics_module
  use IAPWS_module
  use IFC67_module
  use fson
  use fson_mpi_module

  implicit none
  private

  PetscInt, parameter :: max_thermo_ID_length = 8
  character(max_thermo_ID_length), parameter :: &
       default_thermo_ID = "IAPWS"

  public :: setup_thermodynamics

contains

!------------------------------------------------------------------------

  subroutine setup_thermodynamics(json, thermo, logfile)
    !! Reads thermodynamic formulation from JSON input file.
    !! If not present, a default value is assigned.

    use fson_value_m, only : TYPE_STRING
    use utils_module, only : str_to_lower
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    class(thermodynamics_type), allocatable, intent(in out) :: thermo
    type(logfile_type), intent(in out) :: logfile
    ! Locals:
    character(max_thermo_ID_length) :: thermo_ID
    PetscInt :: thermo_type
    PetscBool :: extrapolate
    PetscBool, parameter :: default_extrapolate = PETSC_FALSE

    if (fson_has_mpi(json, "thermodynamics")) then
       thermo_type = fson_type_mpi(json, "thermodynamics")
       if (thermo_type == TYPE_STRING) then
          call fson_get_mpi(json, "thermodynamics", val = thermo_ID)
          extrapolate = default_extrapolate
       else
          call fson_get_mpi(json, "thermodynamics.name", &
               default_thermo_ID, thermo_ID, logfile)
          call fson_get_mpi(json, "thermodynamics.extrapolate", &
               default_extrapolate, extrapolate, logfile)
       end if
    else
       thermo_ID = default_thermo_ID
       call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
            str_key = 'thermodynamics.name', str_value = default_thermo_ID)
       extrapolate = default_extrapolate
       call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
            logical_keys = ['thermodynamics.extrapolate'], &
            logical_values = [default_extrapolate])
    end if

    thermo_ID = str_to_lower(thermo_ID)
    select case (thermo_ID)
    case ("ifc67")
       allocate(IFC67_type :: thermo)
    case default
       allocate(IAPWS_type :: thermo)
    end select

    call thermo%init(extrapolate)

  end subroutine setup_thermodynamics

!------------------------------------------------------------------------

end module thermodynamics_setup_module
