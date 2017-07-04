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

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use logfile_module
  use fson_mpi_module
  use fson

  implicit none
  private

  PetscInt, parameter, public :: max_zone_name_length = 80
  PetscInt, parameter, public :: ZONE_TYPE_CELL_ARRAY = 1
  character(len = 8), public :: zone_label_name = "zone"

  type, public :: zone_type
     !! 3-D zone in the mesh, used to identify cells.
     private
     PetscInt :: index !! Zone index
   contains
     procedure, public :: init => zone_init
     procedure, public :: find_cells => zone_find_cells
     procedure, public :: destroy => zone_destroy
  end type zone_type

  type, public, extends(zone_type) :: zone_cell_array_type
     !! Zone defined by an array of global cell indices.
     PetscInt, allocatable, public :: cells(:)
   contains
     procedure, public :: init => zone_cell_array_init
     procedure, public :: destroy => zone_cell_array_destroy
     procedure, public :: find_cells => zone_cell_array_find_cells
  end type zone_cell_array_type

  public :: get_zone_type

contains

!------------------------------------------------------------------------

  PetscInt function get_zone_type(json) result(ztype)
    !! Determines zone type from JSON input.

    use fson_value_m, only: TYPE_ARRAY, TYPE_OBJECT

    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscInt, parameter :: max_type_str_len = 16
    character(max_type_str_len) :: type_str
    PetscInt :: json_type

    zone_type = -1

    json_type = fson_type_mpi(json, ".")

    select case (json_type)

    case (TYPE_ARRAY)

       zone_type = ZONE_TYPE_CELL_ARRAY

    case (TYPE_OBJECT)

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

    end select

  end function get_zone_type

!------------------------------------------------------------------------
! zone_type
!------------------------------------------------------------------------

  subroutine zone_init(self, index, json)
    !! Initialise zone.

    class(zone_type), intent(in out) :: self
    PetscInt, intent(in) :: index
    type(fson_value), pointer, intent(in) :: json

    self%index = index

  end subroutine zone_init

!------------------------------------------------------------------------

  subroutine zone_destroy(self)
    !! Destroys a zone- dummy routine to be overridden.

    class(zone_type), intent(in out) :: self

    continue

  end subroutine zone_destroy

!------------------------------------------------------------------------

  subroutine zone_find_cells(self, dm, cell_geometry, err)
    !! Find cells in a zone- dummy routine to be overridden.

    class(zone_type), intent(in out) :: self
    DM, intent(in out) :: dm
    Vec, intent(in) :: cell_geometry
    PetscErrorCode, intent(out) :: err

    err = 0

  end subroutine zone_find_cells

!------------------------------------------------------------------------
! zone_cell_array_type
!------------------------------------------------------------------------

  subroutine zone_cell_array_init(self, index, json)
    !! Initialise cell array zone.

    use fson_value_m, only: TYPE_ARRAY, TYPE_OBJECT

    class(zone_cell_array_type), intent(in out) :: self
    PetscInt, intent(in) :: index
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscInt :: json_type
    PetscInt, allocatable :: default_cells(:)

    call self%zone_type%init(index, json)

    json_type = fson_type_mpi(json, ".")

    select case (json_type)
    case (TYPE_ARRAY)
       call fson_get_mpi(json, ".", val = self%cells)
    case (TYPE_OBJECT)
       default_cells = [PetscInt::]
       call fson_get_mpi(json, "cells", default_cells, &
            self%cells)
    end select

  end subroutine zone_cell_array_init

!------------------------------------------------------------------------

  subroutine zone_cell_array_destroy(self)
    !! Destroys a cell array zone.

    class(zone_cell_array_type), intent(in out) :: self

    if (allocated(self%cells)) then
       deallocate(self%cells)
    end if

  end subroutine zone_cell_array_destroy

!------------------------------------------------------------------------

  subroutine zone_cell_array_find_cells(self, dm, cell_geometry, err)
    !! Find cells in a zone from array of global node indices.

    use cell_order_module, only: cell_order_label_name

    class(zone_cell_array_type), intent(in out) :: self
    DM, intent(in out) :: dm
    Vec, intent(in) :: cell_geometry
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: num_matching
    IS :: cell_IS
    PetscInt, pointer :: cells(:)
    PetscInt :: ghost, c, ic, global_cell_index
    DMLabel :: ghost_label
    PetscErrorCode :: ierr

    err = 0

    call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)

    do ic = 1, size(self%cells)
       global_cell_index = self%cells(ic)
       call DMGetStratumSize(dm, cell_order_label_name, &
            global_cell_index, num_matching, ierr); CHKERRQ(ierr)
       if (num_matching > 0) then
          call DMGetStratumIS(dm, cell_order_label_name, &
               global_cell_index, cell_IS, ierr); CHKERRQ(ierr)
          call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
          do c = 1, num_matching
             call DMLabelGetValue(ghost_label, cells(c), ghost, ierr)
             if (ghost < 0) then
                call DMSetLabelValue(dm, zone_label_name, cells(c), &
                     self%index, ierr); CHKERRQ(ierr)
             end if
          end do
          call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
          call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
       end if
    end do

  end subroutine zone_cell_array_find_cells

!------------------------------------------------------------------------  

end module zone_module
