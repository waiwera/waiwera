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
  !! Module for MINC (Multiple INteracting Continua) treatment of fractured media.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module

  implicit none
  private

  character(9), parameter, public :: minc_zone_label_name = "minc_zone"
  character(18), parameter, public :: minc_rocktype_zone_label_name = "minc_rocktype_zone"
  character(10), parameter, public :: minc_level_label_name = "minc_level"

  type, public :: minc_type
     !! MINC parameters for a particular zone.
     private
     PetscInt, public, allocatable :: rocktype_zones(:)
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

  subroutine minc_init(self, json, dm, iminc, str, rock_dict, &
       minc_rocktype_zone_index, logfile, err)
    !! Initialises MINC object from JSON input, and sets minc label on
    !! DM at cells where these MINC parameters are to be applied.

    use fson
    use fson_mpi_module
    use logfile_module
    use fson_value_m, only: TYPE_ARRAY, TYPE_REAL, TYPE_INTEGER, &
         TYPE_STRING, TYPE_OBJECT
    use zone_label_module
    use rock_module, only: max_rockname_length, rock_type_label_name, &
         rock_dict_item_type
    use list_module
    use dictionary_module
    
    class(minc_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    DM, intent(in out) :: dm !! DM to set labels on
    PetscInt, intent(in) :: iminc !! Index of MINC zone (1-based)
    character(*), intent(in) :: str !! Logfile string for current MINC object
    type(dictionary_type), intent(in out) :: rock_dict
    PetscInt, intent(in out) :: minc_rocktype_zone_index !! Index to assign next MINC rocktype zone
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: start_cell, end_cell
    PetscBool :: has_label
    type(fson_value), pointer :: spacing_json, rock_json, rocki_json
    PetscInt :: spacing_type, num_spacings, rock_type, num_rocks, irock
    PetscReal :: fracture_spacing
    PetscReal, allocatable :: fracture_spacing_array(:)
    PetscErrorCode :: ierr
    PetscReal :: fracture_volume
    PetscReal, allocatable :: matrix_volume(:)
    PetscReal, parameter :: default_matrix_volume = 0.9_dp
    PetscInt, parameter :: default_num_fracture_planes = 1
    PetscReal, parameter :: default_fracture_spacing = 50._dp
    PetscReal, parameter :: default_fracture_connection_distance = 0._dp

    err = 0

    if (fson_has_mpi(json, "geometry.fracture.volume")) then
       call fson_get_mpi(json, "geometry.fracture.volume", val = fracture_volume)
       call get_matrix_volumes(json, matrix_volume, 1._dp - fracture_volume)
    else
       call get_matrix_volumes(json, matrix_volume, default_matrix_volume)
       fracture_volume = 1._dp - sum(matrix_volume)
    end if
    self%volume = [fracture_volume, matrix_volume]
    self%volume = self%volume / sum(self%volume)
    self%num_levels = size(matrix_volume)

    call fson_get_mpi(json, "geometry.fracture.planes", default_num_fracture_planes, &
         self%num_fracture_planes, logfile, trim(str) // "geometry.fracture.planes")

    allocate(self%fracture_spacing(self%num_fracture_planes))
    if (fson_has_mpi(json, "geometry.fracture.spacing")) then
       call fson_get_mpi(json, "geometry.fracture.spacing", spacing_json)
       spacing_type = fson_type_mpi(spacing_json, ".")
       select case (spacing_type)
       case (TYPE_REAL, TYPE_INTEGER)
          call fson_get_mpi(json, "geometry.fracture.spacing", default_fracture_spacing, &
               fracture_spacing, logfile, trim(str) // "geometry.fracture.spacing")
          self%fracture_spacing = fracture_spacing
       case (TYPE_ARRAY)
          call fson_get_mpi(json, "geometry.fracture.spacing", [default_fracture_spacing], &
               fracture_spacing_array, logfile, trim(str) // "geometry.fracture.spacing")
          num_spacings = size(fracture_spacing_array)
          self%fracture_spacing(1: num_spacings) = fracture_spacing_array
          if (num_spacings < self%num_fracture_planes) then
             self%fracture_spacing(num_spacings + 1:) = fracture_spacing_array(1)
          end if
       end select
    else
       call fson_get_mpi(json, "geometry.fracture.spacing", default_fracture_spacing, &
            fracture_spacing, logfile, trim(str) // "geometry.fracture.spacing")
       self%fracture_spacing = fracture_spacing
    end if

    call fson_get_mpi(json, "geometry.fracture.connection", &
         default_fracture_connection_distance, &
         self%fracture_connection_distance, logfile, &
         trim(str) // "geometry.fracture.connection")

    call self%setup_geometry(err)

    if (err == 0) then

       if (fson_has_mpi(json, "rock")) then

          call fson_get_mpi(json, "rock", rock_json)
          rock_type = fson_type_mpi(rock_json, ".")
          select case (rock_type)
          case (TYPE_OBJECT)
             num_rocks = 1
             rocki_json => rock_json
          case (TYPE_ARRAY)
             num_rocks = fson_value_count_mpi(rock_json, ".")
             rocki_json => fson_value_children_mpi(rock_json)
          end select
          allocate(self%rocktype_zones(num_rocks))

          call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr); CHKERRQ(ierr)

          do irock = 1, num_rocks
             self%rocktype_zones(irock) = minc_rocktype_zone_index
             call label_rock_zones(iminc, minc_rocktype_zone_index, err)
             if (err == 0) then
                call label_rock_types(iminc, minc_rocktype_zone_index, err)
             else
                exit
             end if
             if (rock_type == TYPE_ARRAY) then
                rocki_json => fson_value_next_mpi(rocki_json)
             end if
             minc_rocktype_zone_index = minc_rocktype_zone_index + 1
          end do

       else
          if (present(logfile)) then
             call logfile%write(LOG_LEVEL_WARN, 'minc', &
                  'init', str_key = trim(str) // 'rock', &
                  str_value = 'not found')
          end if
       end if

    else
       call logfile%write(LOG_LEVEL_ERR, 'minc', 'init', &
            str_key = 'stop', &
            str_value = 'Could not compute MINC geometry.', &
            rank = 0)
    end if

  contains

!........................................................................

    subroutine get_matrix_volumes(json, matrix_volume, default_matrix_volume)

      type(fson_value), pointer, intent(in) :: json
      PetscReal, allocatable, intent(out) :: matrix_volume(:)
      PetscReal, intent(in) :: default_matrix_volume
      ! Locals:
      PetscInt :: matrix_type

      if (fson_has_mpi(json, "geometry.matrix.volume")) then
         matrix_type = fson_type_mpi(json, "geometry.matrix.volume")
         select case (matrix_type)
         case (TYPE_REAL)
            allocate(matrix_volume(1))
            call fson_get_mpi(json, "geometry.matrix.volume", val = matrix_volume(1))
         case (TYPE_ARRAY)
            call fson_get_mpi(json, "geometry.matrix.volume", val = matrix_volume)
         end select
      else
         allocate(matrix_volume(1))
         call fson_get_mpi(json, "geometry.matrix.volume", default_matrix_volume, &
              matrix_volume(1), logfile, trim(str) // "geometry.matrix.volume")
      end if

    end subroutine get_matrix_volumes

!........................................................................

    subroutine label_rock_zones(iminc, minc_rocktype_zone_index, err)
      !! Sets MINC zone and rocktype labels on specified zones.

      PetscInt, intent(in) :: iminc, minc_rocktype_zone_index
      PetscErrorCode, intent(out) :: err
      ! Locals:
      PetscInt :: zones_type, iz, c, ic, num_cells
      PetscInt, pointer :: zone_cells(:)
      character(max_zone_name_length), allocatable :: zones(:)
      character(:), allocatable :: label_name
      IS :: cell_IS
      PetscErrorCode :: ierr

      err = 0

      if (fson_has_mpi(rocki_json, "zones")) then
         zones_type = fson_type_mpi(rocki_json, "zones")
         select case (zones_type)
         case (TYPE_STRING)
            allocate(zones(1))
            call fson_get_mpi(rocki_json, "zones", val = zones(1))
         case (TYPE_ARRAY)
            call fson_get_mpi(rocki_json, "zones", &
                 string_length = max_zone_name_length, val = zones)
         end select
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
                          call DMSetLabelValue(dm, minc_rocktype_zone_label_name, &
                               c, minc_rocktype_zone_index, ierr); CHKERRQ(ierr)
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

    end subroutine label_rock_zones

!........................................................................

    subroutine label_rock_types(iminc, minc_rocktype_zone_index, err)
      !! Sets MINC zone and rocktype labels on specified original rock
      !! types.

      PetscInt, intent(in) :: iminc, minc_rocktype_zone_index
      PetscErrorCode, intent(out) :: err
      ! Locals:
      PetscInt :: itype, ic, c, num_cells
      PetscInt :: types_type
      type(list_node_type), pointer :: node
      character(max_rockname_length), allocatable :: types(:)
      character(max_rockname_length) :: rock_name
      IS :: cell_IS
      PetscInt, pointer :: type_cells(:)
      PetscErrorCode :: ierr

      err = 0

      if (fson_has_mpi(rocki_json, "types")) then
         call DMHasLabel(dm, rock_type_label_name, has_label, ierr)
         CHKERRQ(ierr)
         if (has_label) then
            types_type = fson_type_mpi(rocki_json, "types")
            select case (types_type)
            case (TYPE_STRING)
               allocate(types(1))
               call fson_get_mpi(rocki_json, "types", val = types(1))
            case (TYPE_ARRAY)
               call fson_get_mpi(rocki_json, "types", &
                    string_length = max_rockname_length, val = types)
            end select
            associate(num_types => size(types))
              do itype = 1, num_types
                 rock_name = types(itype)
                 node => rock_dict%get(rock_name)
                 if (associated(node)) then
                    select type (item => node%data)
                    type is (rock_dict_item_type)
                       call DMGetStratumSize(dm, rock_type_label_name, &
                            item%label_value, num_cells, ierr); CHKERRQ(ierr)
                       if (num_cells > 0) then
                          call DMGetStratumIS(dm, rock_type_label_name, &
                               item%label_value, cell_IS, ierr); CHKERRQ(ierr)
                          call ISGetIndicesF90(cell_IS, type_cells, ierr)
                          CHKERRQ(ierr)
                          do ic = 1, num_cells
                             c = type_cells(ic)
                             if ((c >= start_cell) .and. (c < end_cell)) then
                                call DMSetLabelValue(dm, minc_zone_label_name, &
                                     c, iminc, ierr); CHKERRQ(ierr)
                                call DMSetLabelValue(dm, minc_rocktype_zone_label_name, &
                                     c, minc_rocktype_zone_index, ierr); CHKERRQ(ierr)
                             end if
                          end do
                          call ISRestoreIndicesF90(cell_IS, type_cells, ierr)
                          CHKERRQ(ierr)
                          call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
                       end if
                    end select
                 else
                    err = 1
                    if (present(logfile)) then
                       call logfile%write(LOG_LEVEL_ERR, "input", &
                            "unrecognised rock type", &
                            str_key = "name", str_value = rock_name)
                    end if
                    exit
                 end if
              end do
            end associate
         else
            err = 1
            if (present(logfile)) then
               call logfile%write(LOG_LEVEL_ERR, "input", &
                    "no rock types defined for MINC")
            end if
         end if
      end if

    end subroutine label_rock_types

!........................................................................

  end subroutine minc_init

!------------------------------------------------------------------------

  subroutine minc_destroy(self)
    !! Destroys MINC object.

    class(minc_type), intent(in out) :: self

    deallocate(self%volume)
    deallocate(self%fracture_spacing)
    deallocate(self%connection_distance)
    deallocate(self%connection_area)
    if (allocated(self%rocktype_zones)) deallocate(self%rocktype_zones)

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
