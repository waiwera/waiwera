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
  use list_module
  use fson

  implicit none
  private

  PetscInt, parameter, public :: max_zone_name_length = 80
  PetscInt, parameter, public :: ZONE_TYPE_CELL_ARRAY = 1, ZONE_TYPE_BOX = 2, &
       ZONE_TYPE_COMBINE = 3

  type, public :: zone_type
     !! 3-D zone in the mesh, used to identify cells.
     private
     PetscInt :: index !! Zone index
     character(max_zone_name_length), public :: name
   contains
     procedure, public :: init => zone_init
     procedure, public :: label_cells => zone_label_cells
     procedure, public :: destroy => zone_destroy
     procedure, public :: label_name => zone_label_name
  end type zone_type

  type, public, extends(zone_type) :: zone_cell_array_type
     !! Zone defined by an array of global cell indices.
     private
     PetscInt, allocatable, public :: cells(:)
   contains
     procedure, public :: init => zone_cell_array_init
     procedure, public :: destroy => zone_cell_array_destroy
     procedure, public :: label_cells => zone_cell_array_label_cells
  end type zone_cell_array_type

  type, public, extends(zone_type) :: zone_box_type
     !! Zone defined by coordinate ranges. Bounds for any coordinate
     !! may be omitted, in which case they will not be tested. Hence,
     !! a box with no bounds corresponds to all cells in the mesh.
     private
     PetscReal, public :: coord_range(2, 3)
     PetscBool, public :: coord_specified(3)
   contains
     procedure, public :: init => zone_box_init
     procedure, public :: label_cells => zone_box_label_cells
  end type zone_box_type

  type, public, extends(zone_type) :: zone_combine_type
     !! Zone defined by combining other zones, adding zones together
     !! and/or subtracting (excluding) them.
     private
     character(max_zone_name_length), allocatable :: plus(:), minus(:)
   contains
     procedure, public :: init => zone_combine_init
     procedure, public :: destroy => zone_combine_destroy
     procedure, public :: label_cells => zone_combine_label_cells
  end type zone_combine_type

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

    ztype = -1

    json_type = fson_type_mpi(json, ".")

    select case (json_type)

    case (TYPE_ARRAY)

       ztype = ZONE_TYPE_CELL_ARRAY

    case (TYPE_OBJECT)

       if (fson_has_mpi(json, "type")) then

          call fson_get_mpi(json, "type", val = type_str)

          select case (type_str)
          case ('array')
             ztype = ZONE_TYPE_CELL_ARRAY
          case ('box')
             ztype = ZONE_TYPE_BOX
          case ('combine')
             ztype = ZONE_TYPE_COMBINE
          end select

       else ! determine type from object keys:

          if (fson_has_mpi(json, "cells")) then

             ztype = ZONE_TYPE_CELL_ARRAY

          else if (fson_has_mpi(json, "x") .or. &
               fson_has_mpi(json, "r") .or. &
               fson_has_mpi(json, "y") .or. &
               fson_has_mpi(json, "z")) then

             ztype = ZONE_TYPE_BOX

          else if (fson_has_mpi(json, "+") .or. &
               fson_has_mpi(json, "-")) then

             ztype = ZONE_TYPE_COMBINE

          end if

       end if

    end select

  end function get_zone_type

!------------------------------------------------------------------------
! zone_type
!------------------------------------------------------------------------

  subroutine zone_init(self, index, name, json)
    !! Initialise zone.

    class(zone_type), intent(in out) :: self
    PetscInt, intent(in) :: index
    character(*), intent(in) :: name
    type(fson_value), pointer, intent(in) :: json

    self%index = index
    self%name = name

  end subroutine zone_init

!------------------------------------------------------------------------

  subroutine zone_destroy(self)
    !! Destroys a zone- dummy routine to be overridden.

    class(zone_type), intent(in out) :: self

    continue

  end subroutine zone_destroy

!------------------------------------------------------------------------

  subroutine zone_label_cells(self, dm, cell_geometry, err)
    !! Label cells in a zone on DM- dummy routine to be overridden.

    class(zone_type), intent(in out) :: self
    DM, intent(in out) :: dm
    Vec, intent(in) :: cell_geometry
    PetscErrorCode, intent(out) :: err

    err = 0

  end subroutine zone_label_cells

!------------------------------------------------------------------------

  function zone_label_name(self) result(name)
    !! Returns name of zone label associated with a zone.

    class(zone_type), intent(in) :: self
    character(:), allocatable :: name

    allocate(name, source = 'zone_' // self%name)

  end function zone_label_name

!------------------------------------------------------------------------
! zone_cell_array_type
!------------------------------------------------------------------------

  subroutine zone_cell_array_init(self, index, name, json)
    !! Initialise cell array zone.

    use fson_value_m, only: TYPE_ARRAY, TYPE_OBJECT

    class(zone_cell_array_type), intent(in out) :: self
    PetscInt, intent(in) :: index
    character(*), intent(in) :: name
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscInt :: json_type
    PetscInt, allocatable :: default_cells(:)

    call self%zone_type%init(index, name, json)

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

  subroutine zone_cell_array_label_cells(self, dm, cell_geometry, err)
    !! Label cells on a DM in a zone defined by an array of global
    !! node indices.

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
    character(:), allocatable :: label_name
    PetscErrorCode :: ierr

    err = 0

    label_name = self%label_name()
    call DMCreateLabel(dm, label_name, ierr); CHKERRQ(ierr)
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
                call DMSetLabelValue(dm, label_name, cells(c), &
                     1, ierr); CHKERRQ(ierr)
             end if
          end do
          call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
          call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
       end if
    end do

  end subroutine zone_cell_array_label_cells

!------------------------------------------------------------------------
! zone_box_type
!------------------------------------------------------------------------

  subroutine zone_box_init(self, index, name, json)
    !! Initialise box zone.

    class(zone_box_type), intent(in out) :: self
    PetscInt, intent(in) :: index
    character(*), intent(in) :: name
    type(fson_value), pointer, intent(in) :: json

    call self%zone_type%init(index, name, json)

    self%coord_specified = PETSC_FALSE
    call set_coords(1, "x")
    call set_coords(1, "r")
    call set_coords(2, "y")
    call set_coords(3, "z")

  contains

    subroutine set_coords(index, coord)

      PetscInt, intent(in) :: index
      character(*), intent(in) :: coord
      ! Locals:
      PetscReal, allocatable :: r(:)

      if (fson_has_mpi(json, coord)) then
         self%coord_specified(index) = PETSC_TRUE
         call fson_get_mpi(json, coord, val = r)
         self%coord_range(:, index) = r
      end if

    end subroutine set_coords

  end subroutine zone_box_init

!------------------------------------------------------------------------

  subroutine zone_box_label_cells(self, dm, cell_geometry, err)
    !! Label cells on a DM in a zone with specified coordinate ranges.

    use dm_utils_module, only: local_vec_section, section_offset
    use cell_module

    class(zone_box_type), intent(in out) :: self
    DM, intent(in out) :: dm
    Vec, intent(in) :: cell_geometry
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: dim, ghost, i, c, offset
    DMLabel :: ghost_label
    PetscInt :: start_cell, end_cell, end_interior_cell
    PetscInt :: cmax, fmax, emax, vmax
    PetscSection :: section
    PetscReal, contiguous, pointer :: cell_geom_array(:)
    type(cell_type) :: cell
    PetscBool :: found(3)
    character(:), allocatable :: label_name
    PetscErrorCode :: ierr

    err = 0

    label_name = self%label_name()
    call DMCreateLabel(dm, label_name, ierr); CHKERRQ(ierr)

    call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call DMPlexGetHybridBounds(dm, cmax, fmax, emax, vmax, ierr)
    CHKERRQ(ierr)
    end_interior_cell = cmax
    call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call local_vec_section(cell_geometry, section)
    call VecGetArrayReadF90(cell_geometry, cell_geom_array, ierr); CHKERRQ(ierr)

    do c = start_cell, end_interior_cell - 1

       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then

          call section_offset(section, c, offset, ierr); CHKERRQ(ierr)
          call cell%assign_geometry(cell_geom_array, offset)
          found = PETSC_TRUE

          do i = 1, dim
             if (self%coord_specified(i)) then
                found(i) = &
                     (self%coord_range(1, i) <= cell%centroid(i)) .and. &
                     (cell%centroid(i) <= self%coord_range(2, i))
             else
                found(i) = PETSC_TRUE
             end if
          end do

          if (all(found)) then
             call DMSetLabelValue(dm, label_name, c, 1, &
                  ierr); CHKERRQ(ierr)
          end if

       end if
    end do

    call VecRestoreArrayReadF90(cell_geometry, cell_geom_array, ierr)
    CHKERRQ(ierr)

  end subroutine zone_box_label_cells

!------------------------------------------------------------------------
! zone_combine_type
!------------------------------------------------------------------------

  subroutine zone_combine_init(self, index, name, json)
    !! Initialise a combined zone.

    class(zone_combine_type), intent(in out) :: self
    PetscInt, intent(in) :: index
    character(*), intent(in) :: name
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    character(max_zone_name_length), allocatable :: default_names(:)

    call self%zone_type%init(index, name, json)
    default_names = [character(max_zone_name_length)::]

    call fson_get_mpi(json, "+", default_names, max_zone_name_length, &
         self%plus)

    call fson_get_mpi(json, "-", default_names, max_zone_name_length, &
         self%minus)

  end subroutine zone_combine_init

!------------------------------------------------------------------------

  subroutine zone_combine_destroy(self)
    !! Destroys a combined zone.

    class(zone_combine_type), intent(in out) :: self

    deallocate(self%plus, self%minus)

  end subroutine zone_combine_destroy

!------------------------------------------------------------------------

  subroutine zone_combine_label_cells(self, dm, cell_geometry, err)
    !! Label cells in a combined zone on a DM.

    class(zone_combine_type), intent(in out) :: self
    DM, intent(in out) :: dm
    Vec, intent(in) :: cell_geometry
    PetscErrorCode, intent(out) :: err
    ! Locals:
    DMLabel :: ghost_label
    PetscInt :: start_cell, end_cell, end_interior_cell
    PetscInt :: cmax, fmax, emax, vmax
    PetscInt :: i, c, p, num_matching, ghost
    IS :: zone_IS
    PetscInt, pointer :: cells(:)
    character(:), allocatable :: label_name
    PetscErrorCode :: ierr

    err = 0
    label_name = self%label_name()
    call DMCreateLabel(dm, label_name, ierr); CHKERRQ(ierr)

    call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call DMPlexGetHybridBounds(dm, cmax, fmax, emax, vmax, ierr)
    CHKERRQ(ierr)
    end_interior_cell = cmax
    
    associate(num_plus => size(self%plus), num_minus => size(self%minus))

      if ((num_plus == 0) .and. (num_minus > 0)) then
         ! Label all cells:
         do c = start_cell, end_interior_cell - 1
            call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
            if (ghost < 0) then
               call DMSetLabelValue(dm, label_name, c, 1, ierr)
               CHKERRQ(ierr)
            end if
         end do
      else ! Label all zones in plus array:
         do i = 1, num_plus
            call DMGetStratumSize(dm, self%plus(i), &
                 1, num_matching, ierr); CHKERRQ(ierr)
            if (num_matching > 0) then
               call DMGetStratumIS(dm, self%plus(i), &
                    1, zone_IS, ierr); CHKERRQ(ierr)
               call ISGetIndicesF90(zone_IS, cells, ierr); CHKERRQ(ierr)
               do c = 1, num_matching
                  p = cells(c)
                  call DMLabelGetValue(ghost_label, p, ghost, ierr)
                  if (ghost < 0) then
                     call DMSetLabelValue(dm, label_name, p, &
                          1, ierr); CHKERRQ(ierr)
                  end if
               end do
               call ISRestoreIndicesF90(zone_IS, cells, ierr); CHKERRQ(ierr)
            end if
         end do
      end if

      ! Unlabel all zones in minus array:
      do i = 1, num_minus
         call DMGetStratumSize(dm, self%minus(i), &
              1, num_matching, ierr); CHKERRQ(ierr)
         if (num_matching > 0) then
            call DMGetStratumIS(dm, self%minus(i), &
                 1, zone_IS, ierr); CHKERRQ(ierr)
            call ISGetIndicesF90(zone_IS, cells, ierr); CHKERRQ(ierr)
            do c = 1, num_matching
               p = cells(c)
               call DMLabelGetValue(ghost_label, p, ghost, ierr)
               if (ghost < 0) then
                  call DMSetLabelValue(dm, label_name, p, &
                       -1, ierr); CHKERRQ(ierr)
               end if
            end do
            call ISRestoreIndicesF90(zone_IS, cells, ierr); CHKERRQ(ierr)
         end if
      end do

    end associate

  end subroutine zone_combine_label_cells

!------------------------------------------------------------------------  

end module zone_module
