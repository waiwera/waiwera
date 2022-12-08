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

module eos_setup_module
  !! Module for selecting and initialising equation of state (EOS)
  !! modules.

#include <petsc/finclude/petsc.h>

  use petsc

  implicit none
  private

  public :: setup_eos

contains

!------------------------------------------------------------------------
  
  subroutine setup_eos(json, thermo, eos, logfile)
    !! Reads equation of state from JSON input file.  If not present,
    !! a default value is assigned.

    use fson
    use fson_value_m, only : TYPE_STRING, TYPE_OBJECT, TYPE_NULL
    use fson_mpi_module, only: fson_get_mpi, fson_type_mpi
    use logfile_module
    use thermodynamics_module
    use utils_module, only : str_to_lower
    use eos_module
    use eos_w_module
    use eos_we_module
    use eos_wce_module
    use eos_wae_module
    use eos_wse_module
    use eos_wsae_module
    use eos_wsce_module

    type(fson_value), pointer, intent(in) :: json
    class(thermodynamics_type), intent(in) :: thermo
    class(eos_type), allocatable, intent(in out) :: eos
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: eos_json_type
    character(max_eos_name_length), parameter :: &
         default_eos_name = "we"
    character(max_eos_name_length) :: eos_name

    eos_json_type = fson_type_mpi(json, "eos")
    select case (eos_json_type)
    case (TYPE_STRING, TYPE_NULL)
       call fson_get_mpi(json, "eos", default_eos_name, &
            eos_name, logfile)
    case (TYPE_OBJECT)
       call fson_get_mpi(json, "eos.name", default_eos_name, &
            eos_name, logfile)
    end select
    eos_name = str_to_lower(eos_name)

    select case (eos_name)
    case ("w")
       allocate(eos_w_type :: eos)
    case ("we")
       allocate(eos_we_type :: eos)
    case ("wce")
       allocate(eos_wce_type :: eos)
    case ("wae")
       allocate(eos_wae_type :: eos)
    case ("wse")
       allocate(eos_wse_type :: eos)
    case ("wsae")
       allocate(eos_wsae_type :: eos)
    case ("wsce")
       allocate(eos_wsce_type :: eos)
    case default
       allocate(eos_we_type :: eos)
    end select

    call eos%init(json, thermo, logfile)

  end subroutine setup_eos

!------------------------------------------------------------------------

end module eos_setup_module
