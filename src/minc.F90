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

  character(9), parameter, public :: minc_zone_label_name = "minc_zone"
  character(10), parameter, public :: minc_level_label_name = "minc_level"

  type, public :: minc_type
     !! MINC parameters for a particular zone.
     private
     PetscInt, public :: num_levels !! Number of MINC levels (additional MINC cells per fracture cell)
     PetscReal, allocatable, public :: volume(:) !! Volume occupied by fracture and matrix levels (scaled by original cell volume)
     PetscReal, public :: fracture_connection_distance !! Connection distance from fracture to matrix
     PetscReal, allocatable, public :: connection_distance(:) !! Connection distances for matrix cells
     PetscReal, allocatable, public :: connection_area(:) !! Areas for matrix cell connections (scaled by original cell volume)
     PetscInt, public :: num_fracture_planes !! Number of fracture planes (1.. 3)
     PetscReal, allocatable, public :: fracture_spacing(:) !! Fracture spacings for each set of fracture planes
   contains
     private
     procedure, public :: init => minc_init
     procedure, public :: destroy => minc_destroy
     procedure, public :: proximity => minc_proximity
     procedure, public :: proximity_derivative => minc_proximity_derivative
     procedure, public :: inner_connection_distance => minc_inner_connection_distance
     procedure, public :: setup_geometry => minc_setup_geometry
  end type minc_type

contains

!------------------------------------------------------------------------

  subroutine minc_init(self, json, dm, iminc, str, logfile, err)
    !! Initialises MINC object from JSON input, and sets minc label on
    !! DM at cells where these MINC parameters are to be applied.

    use fson
    use fson_mpi_module
    use logfile_module
    use fson_value_m, only : TYPE_ARRAY, TYPE_REAL, TYPE_INTEGER
    use zone_label_module
    
    class(minc_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    DM, intent(in out) :: dm !! DM to set labels on
    PetscInt, intent(in) :: iminc !! Index of MINC zone (1-based)
    character(*), intent(in) :: str !! Logfile string for current MINC object
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: start_cell, end_cell
    PetscInt :: num_cells, ic, c, iz
    PetscInt, allocatable :: cells(:)
    character(max_zone_name_length), allocatable :: zones(:)
    character(:), allocatable :: label_name
    IS :: cell_IS
    PetscInt, pointer :: zone_cells(:)
    PetscBool :: has_label
    type(fson_value), pointer :: spacing_json
    PetscInt :: spacing_type, num_spacings
    PetscReal :: fracture_spacing
    PetscReal, allocatable :: fracture_spacing_array(:)
    PetscErrorCode :: ierr
    PetscReal, parameter :: default_volume_fraction(2) = [0.1_dp, 0.9_dp]
    PetscInt, parameter :: default_num_fracture_planes = 1
    PetscReal, parameter :: default_fracture_spacing = 50._dp
    PetscReal, parameter :: default_fracture_connection_distance = 0._dp

    err = 0

    call fson_get_mpi(json, "volume_fractions", default_volume_fraction, &
         self%volume, logfile, trim(str) // "volume_fractions")
    self%volume = self%volume / sum(self%volume)
    self%num_levels = size(self%volume) - 1

    call fson_get_mpi(json, "fracture.planes", default_num_fracture_planes, &
         self%num_fracture_planes, logfile, trim(str) // "fracture.planes")

    allocate(self%fracture_spacing(self%num_fracture_planes))
    if (fson_has_mpi(json, "fracture.spacing")) then
       call fson_get_mpi(json, "fracture.spacing", spacing_json)
       spacing_type = fson_type_mpi(spacing_json, ".")
       select case (spacing_type)
       case (TYPE_REAL, TYPE_INTEGER)
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

    call fson_get_mpi(json, "fracture.connection_distance", &
         default_fracture_connection_distance, &
         self%fracture_connection_distance, logfile, &
         trim(str) // "fracture.connection_distance")

    call self%setup_geometry(err)

    if (err == 0) then

       if (fson_has_mpi(json, "cells")) then
          call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
          call fson_get_mpi(json, "cells", val = cells)
          num_cells = size(cells)
          if (num_cells > 0) then
             call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
             CHKERRQ(ierr)
             do ic = 1, num_cells
                c = cells(ic)
                if ((c >= start_cell) .and. (c < end_cell)) then
                   call DMSetLabelValue(dm, minc_zone_label_name, &
                        c, iminc, ierr); CHKERRQ(ierr)
                end if
             end do
          end if
          deallocate(cells)
       end if

       if (fson_has_mpi(json, "zones")) then
          call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
          call fson_get_mpi(json, "zones", string_length = max_zone_name_length, &
               val = zones)
          associate(num_zones => size(zones))
            do iz = 1, num_zones
               label_name = zone_label_name(zones(iz))
               call DMHasLabel(dm, label_name, has_label, ierr); CHKERRQ(ierr)
               if (has_label) then
                  call DMGetStratumSize(dm, label_name, 1, num_cells, &
                       ierr); CHKERRQ(ierr)
                  if (num_cells > 0) then
                     call DMGetStratumIS(dm, label_name, 1, cell_IS, &
                          ierr); CHKERRQ(ierr)
                     call ISGetIndicesF90(cell_IS, zone_cells, ierr)
                     CHKERRQ(ierr)
                     do ic = 1, num_cells
                        c = zone_cells(ic)
                        if ((c >= start_cell) .and. (c < end_cell)) then
                           call DMSetLabelValue(dm, minc_zone_label_name, &
                                c, iminc, ierr); CHKERRQ(ierr)
                        end if
                     end do
                     call ISRestoreIndicesF90(cell_IS, zone_cells, ierr)
                     CHKERRQ(ierr)
                     call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
                  end if
               else
                  err = 1
                  if (present(logfile)) then
                     call logfile%write(LOG_LEVEL_ERR, "input", &
                          "unrecognised zone", &
                          str_key = "name", str_value = zones(iz))
                  end if
                  exit
               end if
            end do
          end associate
          deallocate(zones)
       end if

    else
       call logfile%write(LOG_LEVEL_ERR, 'minc', 'init', &
            str_key = 'stop', &
            str_value = 'Could not compute MINC geometry.', &
            rank = 0)
    end if

    ! TODO: rock properties
    
  end subroutine minc_init

!------------------------------------------------------------------------

  subroutine minc_destroy(self)
    !! Destroys MINC object.

    class(minc_type), intent(in out) :: self

    deallocate(self%volume)
    deallocate(self%fracture_spacing)
    deallocate(self%connection_distance)
    deallocate(self%connection_area)

  end subroutine minc_destroy

!------------------------------------------------------------------------

  PetscReal function minc_proximity(self, d) result(p)
    !! MINC proximity function: total fraction of matrix rock within
    !! distance d of a fracture, for 'nested cube' geometry.

    class(minc_type), intent(in) :: self
    PetscReal, intent(in) :: d

    ! fout is an array of the fraction not within d of a fracture, for
    ! each set of fracture planes:
    associate(fout => 1._dp - 2._dp * d / self%fracture_spacing)
      if (any(fout < 0._dp)) then
         p = 1._dp
      else
         p = 1._dp - product(fout)
      end if
    end associate

  end function minc_proximity

!------------------------------------------------------------------------

  PetscReal function minc_proximity_derivative(self, d) result(pd)
    !! Derivative of MINC proximity function with respect to distance
    !! d, for 'nested cube' geometry.

    use utils_module, only: array_exclusive_products

    class(minc_type), intent(in) :: self
    PetscReal, intent(in) :: d

    associate(fout => 1._dp - 2._dp * d / self%fracture_spacing)
      if (any(fout < 0._dp)) then
         pd = 0._dp
      else
         pd = 2._dp * sum(array_exclusive_products(fout) / &
              self%fracture_spacing)
      end if
    end associate

  end function minc_proximity_derivative

!------------------------------------------------------------------------

  PetscReal function minc_inner_connection_distance(self, x) result(c)
    !! Connection distance for innermost cell, at distance x, for
    !! 'nested cube' geometry. From Pruess (1983), "GMINC - a mesh
    !! generator for flow simulations in fractured reservoirs",
    !! report LBL-15227.

    use utils_module, only: array_pair_sum

    class(minc_type), intent(in) :: self
    PetscReal, intent(in) :: x

    associate(u => self%fracture_spacing - 2._dp * x)
      select case (self%num_fracture_planes)
      case (1)
         c = u(1) / 6._dp
      case (2)
         c = 0.25_dp * product(u) / sum(u)
      case (3)
         c = 0.3_dp * product(u) / array_pair_sum(u)
      end select
    end associate

  end function minc_inner_connection_distance

!------------------------------------------------------------------------

  subroutine minc_setup_geometry(self, err)
    !! Calculates distances and areas for connections between matrix
    !! cells.

    use utils_module, only: array_cumulative_sum
    use root_finder_module

    class(minc_type), intent(in out) :: self
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal, target :: volsum(self%num_levels)
    PetscReal :: x, xl, xr
    PetscInt :: i
    type(root_finder_type) :: root_finder
    procedure(root_finder_function), pointer :: f
    class(*), pointer :: v

    err = 0
    allocate(self%connection_distance(self%num_levels + 1))
    allocate(self%connection_area(self%num_levels))

    f => volume_difference
    call root_finder%init(f)

    associate(vmatrix => 1._dp - self%volume(1))

      volsum = array_cumulative_sum(self%volume(2:)) / vmatrix

      x = 0._dp
      self%connection_distance(1) = self%fracture_connection_distance
      self%connection_area(1) = vmatrix * self%proximity_derivative(x)
      xr = self%volume(2) / self%connection_area(1)

      do i = 1, self%num_levels - 1

         xl = x
         v => volsum(i)
         do while (f(xr, v) < 0._dp)
            xr = xr * 2._dp
         end do
         root_finder%interval = [xl, xr]
         root_finder%context => v
         call root_finder%find()
         if (root_finder%err == 0) then
            x = root_finder%root
         else
            err = 1
            exit
         end if

         self%connection_distance(i + 1) = 0.5_dp * (x - xl)
         self%connection_area(i + 1) = vmatrix * self%proximity_derivative(x)

      end do

      if (err == 0) then
         self%connection_distance(self%num_levels + 1) = &
              self%inner_connection_distance(x)
      end if

    end associate

    call root_finder%destroy()
    f => null()
    v => null()

  contains

    PetscReal function volume_difference(x, context) result(y)

      PetscReal, intent(in) :: x
      class(*), pointer, intent(in out) :: context

      select type (context)
      type is (PetscReal)
         associate(v => context)
           y = self%proximity(x) - v
         end associate
      end select

    end function volume_difference

  end subroutine minc_setup_geometry

!------------------------------------------------------------------------

end module minc_module
