!   Copyright 2016 University of Auckland.

!   This file is part of Waiwera.

!   Waiwera is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.

!   Waiwera is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.

!   You should have received a copy of the GNU General Public License
!   along with Waiwera.  If not, see <http://www.gnu.org/licenses/>.

module eos_setup_module
  !! Module for selecting and initialising equation of state (EOS)
  !! modules.

  implicit none
  private

  public :: setup_eos

contains

!------------------------------------------------------------------------
  
  subroutine setup_eos(json, thermo, eos, logfile)
    !! Reads equation of state from JSON input file.  If not present,
    !! a default value is assigned.

    use fson
    use fson_mpi_module, only: fson_get_mpi
    use logfile_module
    use thermodynamics_module
    use utils_module, only : str_to_lower
    use eos_module
    use eos_w_module
    use eos_we_module

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
