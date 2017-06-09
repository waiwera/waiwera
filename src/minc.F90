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

module minc_module
  !! Module for MINC (Multiple INteracting Continua) treatment of fractured media .

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module

  implicit none
  private

  PetscInt, parameter :: max_minc_label_name_length = 4
  character(max_minc_label_name_length), parameter, public :: minc_label_name = "minc"

  type, public :: minc_type
     !! MINC parameters for a particular zone.
     private
     PetscInt, public :: num_levels !! Number of MINC levels (additional MINC cells per fracture cell)
     PetscReal, allocatable, public :: volume_fractions(:) !! Fraction of cell volume occupied by fractures
     PetscInt, public :: num_fracture_planes !! Number of fracture planes (1.. 3)
     PetscReal, allocatable, public :: fracture_spacing(:) !! Fracture spacings for each set of fracture planes
   contains
     private
     procedure, public :: init => minc_init
     procedure, public :: destroy => minc_destroy
  end type minc_type

contains

!------------------------------------------------------------------------

  subroutine minc_init(self, json, dm, iminc, str, logfile)
    !! Initialises MINC object from JSON input, and sets minc label on
    !! DM at cells where these MINC parameters are to be applied.

    use fson
    use fson_mpi_module
    use logfile_module
    use fson_value_m, only : TYPE_ARRAY, TYPE_REAL
    
    class(minc_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    DM, intent(in out) :: dm !! DM to set labels on
    PetscInt, intent(in) :: iminc !! Index of MINC zone (1-based)
    character(*), intent(in) :: str !! Logfile string for current MINC object
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    DMLabel :: minc_label
    PetscInt :: start_cell, end_cell
    PetscInt :: num_cells, ic, c
    PetscInt, allocatable :: cells(:)
    PetscInt, allocatable :: default_cells(:)
    type(fson_value), pointer :: spacing_json
    PetscInt :: spacing_type, num_spacings
    PetscReal :: fracture_spacing
    PetscReal, allocatable :: fracture_spacing_array(:)
    PetscErrorCode :: ierr
    PetscReal, parameter :: default_volume_fractions(2) = [0.1_dp, 0.9_dp]
    PetscInt, parameter :: default_num_fracture_planes = 1
    PetscReal, parameter :: default_fracture_spacing = 50._dp

    call fson_get_mpi(json, "volume_fractions", default_volume_fractions, &
         self%volume_fractions, logfile, trim(str) // "volume_fractions")
    self%volume_fractions = self%volume_fractions / sum(self%volume_fractions)
    self%num_levels = size(self%volume_fractions) - 1

    call fson_get_mpi(json, "fracture.planes", default_num_fracture_planes, &
         self%num_fracture_planes, logfile, trim(str) // "fracture.planes")

    allocate(self%fracture_spacing(self%num_fracture_planes))
    if (fson_has_mpi(json, "fracture.spacing")) then
       call fson_get_mpi(json, "fracture.spacing", spacing_json)
       spacing_type = fson_type_mpi(spacing_json, ".")
       select case (spacing_type)
       case (TYPE_REAL)
          call fson_get_mpi(json, "fracture.spacing", default_fracture_spacing, &
               fracture_spacing, logfile, trim(str) // "fracture.spacing")
          self%fracture_spacing = fracture_spacing
       case (TYPE_ARRAY)
          call fson_get_mpi(json, "fracture.spacing", [default_fracture_spacing], &
               fracture_spacing_array, logfile, trim(str) // "fracture.spacing")
          num_spacings = size(fracture_spacing_array)
          self%fracture_spacing(1: num_spacings) = fracture_spacing_array
          if (num_spacings < self%num_fracture_planes) then
             self%fracture_spacing(num_spacings + 1:) = fracture_spacing_array(1)
          end if
       end select
    else
       call fson_get_mpi(json, "fracture.spacing", default_fracture_spacing, &
            fracture_spacing, logfile, trim(str) // "fracture.spacing")
       self%fracture_spacing = fracture_spacing
    end if

    call DMGetLabel(dm, minc_label_name, minc_label, ierr); CHKERRQ(ierr)
    default_cells = [PetscInt::]
    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call fson_get_mpi(json, "cells", default_cells, cells, logfile, &
         trim(str) // "cells")
    if (allocated(cells)) then
       num_cells = size(cells)
       do ic = 1, num_cells
          c = cells(ic)
          if ((c >= start_cell) .and. (c < end_cell)) then
             call DMSetLabelValue(dm, minc_label_name, &
                  c, iminc, ierr); CHKERRQ(ierr)
          end if
       end do
       deallocate(cells)
    end if

    ! TODO: rock properties
    
  end subroutine minc_init

!------------------------------------------------------------------------

  subroutine minc_destroy(self)
    !! Destroys MINC object.

    class(minc_type), intent(in out) :: self

    deallocate(self%volume_fractions)
    deallocate(self%fracture_spacing)

  end subroutine minc_destroy

!------------------------------------------------------------------------

end module minc_module
