!   Copyright 2017 University of Auckland.

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

module zone_module
  !! Module for a zone in the mesh and identifying which cells are inside it.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use logfile_module
  use fson_mpi_module
  use fson

  implicit none
  private

  PetscInt, parameter, public :: max_zone_name_length = 80
  PetscInt, parameter, public :: ZONE_TYPE_CELL_ARRAY = 1

  type, public :: zone_type
     !! 3-D zone in the mesh, used to identify cells.
     private
     character(max_zone_name_length), public :: name
     PetscInt, allocatable, public :: cells(:)
   contains
     procedure, public :: init => zone_init
     procedure, public :: find_cells => zone_find_cells
     procedure, public :: destroy => zone_destroy
  end type zone_type

  type, public, extends(zone_type) :: zone_cell_array_type
     !! Zone defined by an array of cells. (Note: it is not necessary
     !! to override the find_cells() method, because the cells are
     !! explicitly specified in the input.)
   contains
     procedure, public :: init => zone_cell_array_init
  end type zone_cell_array_type

  public :: get_zone_type

contains

!------------------------------------------------------------------------

  PetscInt function get_zone_type(json) result(zone_type)
    !! Determines zone type from JSON input.

    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscInt, parameter :: max_type_str_len = 16
    character(max_type_str_len) :: type_str

    zone_type = -1

    if (fson_has_mpi(json, "type")) then

       call fson_get_mpi(json, "type", val = type_str)

       select case (type_str)
       case ('array')
          zone_type = ZONE_TYPE_CELL_ARRAY
       end select

    else ! determine type from object keys:

       if (fson_has_mpi(json, "cells")) then

          zone_type = ZONE_TYPE_CELL_ARRAY

       end if

    end if

  end function get_zone_type

!------------------------------------------------------------------------

  subroutine zone_init(self, json, logfile, err)
    !! Initialise zone. This routine just gets the zone name - to be
    !! overridden by derived types.

    class(zone_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err

    if (fson_has_mpi(json, "name")) then
       call fson_get_mpi(json, "name", val = self%name)
       err = 0
    else
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_ERR, 'zone', 'init', &
               str_key = 'stop', &
               str_value = 'unnamed zone found.')
       end if
       err = 1
    end if

  end subroutine zone_init

!------------------------------------------------------------------------

  subroutine zone_destroy(self)
    !! Destroys a zone.

    class(zone_type), intent(in out) :: self

    if (allocated(self%cells)) then
       deallocate(self%cells)
    end if

  end subroutine zone_destroy

!------------------------------------------------------------------------

  subroutine zone_find_cells(self, err)
    !! Find cells in a zone (dummy routine to be overridden by derived
    !! types).

    class(zone_type), intent(in out) :: self
    PetscErrorCode, intent(out) :: err

    err = 0
    
  end subroutine zone_find_cells

!------------------------------------------------------------------------
! zone_cell_array_type
!------------------------------------------------------------------------

  subroutine zone_cell_array_init(self, json, logfile, err)
    !! Initialise cell array zone.

    use fson_value_m, only: TYPE_ARRAY, TYPE_OBJECT

    class(zone_cell_array_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: json_type
    PetscInt, allocatable :: default_cells(:)

    call self%zone_type%init(json, logfile, err)

    if (err == 0) then

       json_type = fson_type_mpi(json, ".")

       select case (json_type)
       case (TYPE_ARRAY)
          call fson_get_mpi(json, ".", val = self%cells)
       case (TYPE_OBJECT)
          default_cells = [PetscInt::]
          call fson_get_mpi(json, "cells", default_cells, &
               self%cells, logfile)
       end select

    end if

  end subroutine zone_cell_array_init

!------------------------------------------------------------------------  

end module zone_module
