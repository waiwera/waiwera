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

module rock_module
  !! Defines type for accessing local rock properties on cells and faces.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use relative_permeability_module
  use capillary_pressure_module
  use fson

  implicit none
  private

  type rock_dict_item_type
     !! Item in rock dictionary.
     private
     type(fson_value), pointer, public :: rock !! JSON input data value
     PetscInt, public :: label_value !! DM rock label value
  end type rock_dict_item_type

  type rock_type
     !! Local rock properties.
     private
     PetscReal, pointer, contiguous, public :: permeability(:)   !! Permeability
     PetscReal, pointer, public :: wet_conductivity, dry_conductivity !! Heat conductivities
     PetscReal, pointer, public :: porosity          !! Porosity
     PetscReal, pointer, public :: density           !! Grain density
     PetscReal, pointer, public :: specific_heat     !! Specific heat
     class(relative_permeability_type), pointer, &
          public :: relative_permeability !! Relative permeability functions
     class(capillary_pressure_type), pointer, &
          public :: capillary_pressure !! Capillary pressure function
     PetscInt, public :: dof !! Number of degrees of freedom
   contains
     private
     procedure, public :: init => rock_init
     procedure, public :: assign => rock_assign
     procedure, public :: assign_relative_permeability => &
          rock_assign_relative_permeability
     procedure, public :: assign_capillary_pressure => &
          rock_assign_capillary_pressure
     procedure, public :: destroy => rock_destroy
     procedure, public :: energy => rock_energy
  end type rock_type

  PetscInt, parameter :: num_rock_variables = 6
  PetscInt, parameter :: max_rock_variable_name_length = 32
  character(max_rock_variable_name_length), parameter, public :: &
       rock_variable_names(num_rock_variables) = &
       [character(max_rock_variable_name_length):: &
       "permeability", "wet_conductivity", "dry_conductivity", "porosity", &
       "density", "specific_heat"]
  PetscInt, parameter, public :: &
       rock_variable_num_components(num_rock_variables) = &
       [3, 1, 1, 1, 1, 1]
  PetscInt, parameter, public :: max_rockname_length = 24

  ! Default rock properties:
  PetscReal, parameter, public :: default_permeability_scalar = 1.e-13_dp
  PetscReal, parameter, public :: default_permeability(3) = [ &
       default_permeability_scalar, default_permeability_scalar, &
       default_permeability_scalar]
  PetscReal, parameter, public :: default_porosity = 0.1_dp
  PetscReal, parameter, public :: default_density = 2200.0_dp
  PetscReal, parameter, public :: default_specific_heat = 1000._dp
  PetscReal, parameter, public :: default_heat_conductivity = 2.5_dp

  character(len = 9), public :: rock_type_label_name = "rock_type"

  public :: rock_dict_item_type, rock_type, setup_rock_types, setup_rock_vector, &
       create_rock_vector, label_rock_cell

contains

!------------------------------------------------------------------------

  subroutine rock_init(self)
    !! Initialises rock object.

    class(rock_type), intent(in out) :: self

    self%dof = sum(rock_variable_num_components)

  end subroutine rock_init

!------------------------------------------------------------------------

  subroutine rock_assign(self, data, offset)
    !! Assigns pointers in a rock object to elements of the specified
    !! data array, starting from the given offset.

    class(rock_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! rock data array
    PetscInt, intent(in) :: offset !! rock array offset

    self%permeability => data(offset: offset + 2)
    self%wet_conductivity => data(offset + 3)
    self%dry_conductivity => data(offset + 4)
    self%porosity => data(offset + 5)
    self%density => data(offset + 6)
    self%specific_heat => data(offset + 7)

  end subroutine rock_assign
    
!------------------------------------------------------------------------

  subroutine rock_assign_relative_permeability(self, relative_permeability)
    !! Assigns relative permeability pointer for a rock object.

    class(rock_type), intent(in out) :: self
    class(relative_permeability_type), intent(in), &
         target :: relative_permeability

    self%relative_permeability => relative_permeability

  end subroutine rock_assign_relative_permeability

!------------------------------------------------------------------------

  subroutine rock_assign_capillary_pressure(self, capillary_pressure)
    !! Assigns capillary pressure pointer for a rock object.

    class(rock_type), intent(in out) :: self
    class(capillary_pressure_type), intent(in), &
         target :: capillary_pressure

    self%capillary_pressure => capillary_pressure

  end subroutine rock_assign_capillary_pressure

!------------------------------------------------------------------------

  subroutine rock_destroy(self)
    !! Destroys a rock object (nullifies all pointer components).

    class(rock_type), intent(in out) :: self

    self%permeability => null()
    self%wet_conductivity => null()
    self%dry_conductivity => null()
    self%porosity => null()
    self%density => null()
    self%specific_heat => null()
    self%relative_permeability => null()
    self%capillary_pressure => null()

  end subroutine rock_destroy

!------------------------------------------------------------------------

  PetscReal function rock_energy(self, temperature)
    !! Returns rock energy density at a given temperature.

    class(rock_type), intent(in) :: self
    PetscReal, intent(in) :: temperature !! Temperature (deg C)

    rock_energy = self%density * self%specific_heat * temperature

  end function rock_energy

!------------------------------------------------------------------------

  subroutine set_default_rock_properties(dm, rock_vector, range_start, &
       ghost_cell)
    !! Assigns default rock properties to all cells.

    use dm_utils_module, only: global_section_offset, global_vec_section

    DM, intent(in) :: dm
    Vec, intent(in out) :: rock_vector
    PetscInt, intent(in) :: range_start
    PetscInt, allocatable, intent(in) :: ghost_cell(:)
    ! Locals:
    PetscSection :: rock_section
    PetscReal, pointer, contiguous :: rock_array(:)
    PetscInt :: start_cell, end_cell, c, offset
    type(rock_type) :: rock
    PetscErrorCode :: ierr

    call global_vec_section(rock_vector, rock_section)
    call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call rock%init()

    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)

    do c = start_cell, end_cell - 1
       if (ghost_cell(c) < 0) then
          offset = global_section_offset(rock_section, c, range_start)
          call rock%assign(rock_array, offset)
          rock%permeability = default_permeability
          rock%wet_conductivity = default_heat_conductivity
          rock%dry_conductivity = default_heat_conductivity
          rock%porosity = default_porosity
          rock%density = default_density
          rock%specific_heat = default_specific_heat
       end if
    end do

    call rock%destroy()
    call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)

  end subroutine set_default_rock_properties

!------------------------------------------------------------------------

  subroutine label_rock_cell(dm, start_cell, end_cell, p, rock_type_index)
    !! Sets DM rocktype label on single cell.

    DM, intent(in out) :: dm
    PetscInt, intent(in) :: start_cell, end_cell
    PetscInt, intent(in) :: p, rock_type_index
    ! Locals:
    PetscErrorCode :: ierr

    if ((p >= start_cell) .and. (p < end_cell)) then
       call DMSetLabelValue(dm, rock_type_label_name, &
            p, rock_type_index, ierr); CHKERRQ(ierr)
    end if

  end subroutine label_rock_cell

!------------------------------------------------------------------------

  subroutine setup_rock_types(json, dm, rock_dict, logfile, err)
    !! Sets up rock type dictionary and rock type labels on serial DM, from
    !! JSON input.

    use dictionary_module
    use fson
    use fson_mpi_module
    use fson_value_m, only : TYPE_STRING, TYPE_ARRAY
    use logfile_module
    use zone_label_module

    type(fson_value), pointer, intent(in) :: json
    DM, intent(in out) :: dm
    type(dictionary_type), intent(in out) :: rock_dict
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: start_cell, end_cell, num_rocktypes, ir
    DMLabel :: ghost_label
    type(fson_value), pointer :: rocktypes, r
    character(len=64) :: rockstr
    character(len=12) :: irstr
    character(max_rockname_length) :: name
    type(rock_dict_item_type), pointer :: item
    PetscErrorCode :: ierr

    err = 0
    call rock_dict%init(owner = PETSC_TRUE)

    if (fson_has_mpi(json, "rock.types")) then

       call DMCreateLabel(dm, rock_type_label_name, ierr); CHKERRQ(ierr)

       call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
       CHKERRQ(ierr)
       call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)

       call fson_get_mpi(json, "rock.types", rocktypes)
       num_rocktypes = fson_value_count_mpi(rocktypes, ".")
       r => fson_value_children_mpi(rocktypes)

       do ir = 1, num_rocktypes
          write(irstr, '(i0)') ir - 1
          rockstr = 'rock.types[' // trim(irstr) // '].'
          call fson_get_mpi(r, "name", "", name, logfile, trim(rockstr) // "name")
          if (name /= "") then
             allocate(item)
             item%label_value = ir
             item%rock => r
             call rock_dict%add(name, item)
          end if
          call label_rock_zones(ir, err)
          if (err > 0) exit
          r => fson_value_next_mpi(r)
       end do

    else
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_WARN, "input", "no rocktypes")
       end if
    end if

  contains

    subroutine label_rock_zones(ir, err)
      !! Sets DM rocktype label in specified zones.

      PetscInt, intent(in) :: ir
      PetscErrorCode, intent(out) :: err
      ! Locals:
      PetscInt :: iz, num_zone_cells, c
      PetscBool :: has_label
      character(max_zone_name_length), allocatable :: zones(:)
      PetscInt :: zones_type
      character(:), allocatable :: label_name
      IS :: cell_IS
      PetscInt, pointer :: cells(:)
      PetscErrorCode :: ierr

      err = 0

      if (fson_has_mpi(r, "zones")) then
         zones_type = fson_type_mpi(r, "zones")
         select case (zones_type)
         case (TYPE_STRING)
            allocate(zones(1))
            call fson_get_mpi(r, "zones", val = zones(1))
         case (TYPE_ARRAY)
            call fson_get_mpi(r, "zones", string_length = max_zone_name_length, &
                 val = zones)
         end select
         associate(num_zones => size(zones))
           do iz = 1, num_zones
              label_name = zone_label_name(zones(iz))
              call DMHasLabel(dm, label_name, has_label, ierr); CHKERRQ(ierr)
              if (has_label) then
                 call DMGetStratumSize(dm, label_name, 1, num_zone_cells, &
                      ierr); CHKERRQ(ierr)
                 if (num_zone_cells > 0) then
                    call DMGetStratumIS(dm, label_name, 1, cell_IS, &
                         ierr); CHKERRQ(ierr)
                    call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
                    do c = 1, num_zone_cells
                       call label_rock_cell(dm, start_cell, end_cell, cells(c), ir)
                    end do
                    call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
                    call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
                 end if
              else
                 err = 1
                 if (present(logfile)) then
                    call logfile%write(LOG_LEVEL_ERR, "input", "unrecognised zone", &
                         str_key = "name", str_value = zones(iz))
                 end if
                 exit
              end if
           end do
         end associate
         deallocate(zones)
      end if

    end subroutine label_rock_zones

  end subroutine setup_rock_types

!------------------------------------------------------------------------

  subroutine setup_rock_vector_types(json, dm, rock_vector, rock_dict, &
       range_start, logfile)
    !! Sets up rock vector on DM from rock types in JSON input.

    use dictionary_module
    use dm_utils_module, only: global_section_offset, global_vec_section
    use zone_label_module
    use fson
    use fson_mpi_module
    use logfile_module
    use fson_value_m, only : TYPE_ARRAY

    type(fson_value), pointer, intent(in) :: json
    DM, intent(in out) :: dm
    Vec, intent(in out) :: rock_vector
    type(dictionary_type), intent(in out) :: rock_dict
    PetscInt, intent(in) :: range_start
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: num_rocktypes, ic, num_cells, offset, ghost
    PetscInt :: ir, permeability_type, dim
    PetscInt :: perm_size
    DMLabel :: ghost_label
    type(fson_value), pointer :: rocktypes, r
    PetscInt, pointer :: cells(:)
    type(rock_type) :: rock
    character(max_rockname_length) :: name
    IS :: cell_IS
    PetscReal :: porosity, density, specific_heat
    PetscReal :: wet_conductivity, dry_conductivity
    PetscReal :: permeability_scalar
    PetscReal, allocatable :: permeability(:)
    PetscReal, pointer, contiguous :: rock_array(:)
    PetscSection :: section
    PetscErrorCode :: ierr
    character(len=64) :: rockstr
    character(len=12) :: irstr

    if (fson_has_mpi(json, "rock.types")) then

       call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
       call rock%init()

       call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
       call global_vec_section(rock_vector, section)

       call fson_get_mpi(json, "rock.types", rocktypes)
       num_rocktypes = fson_value_count_mpi(rocktypes, ".")
       r => fson_value_children_mpi(rocktypes)
       call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)

       do ir = 1, num_rocktypes

          write(irstr, '(i0)') ir - 1
          rockstr = 'rock.types[' // trim(irstr) // '].'
          call fson_get_mpi(r, "name", "", name, logfile, trim(rockstr) // "name")
          permeability_type = fson_type_mpi(r, "permeability")
          select case (permeability_type)
          case (TYPE_ARRAY)
             call fson_get_mpi(r, "permeability", default_permeability, &
                  permeability, logfile, trim(rockstr) // "permeability")
          case default ! real or null
             call fson_get_mpi(r, "permeability", default_permeability_scalar, &
                  permeability_scalar, logfile, trim(rockstr) // "permeability")
             allocate(permeability(dim))
             permeability = permeability_scalar
          end select
          call fson_get_mpi(r, "wet_conductivity", default_heat_conductivity, &
               wet_conductivity, logfile, trim(rockstr) // "wet_conductivity")
          call fson_get_mpi(r, "dry_conductivity", wet_conductivity, &
               dry_conductivity, logfile, trim(rockstr) // "dry_conductivity")
          call fson_get_mpi(r, "porosity", default_porosity, porosity, logfile, &
               trim(rockstr) // "porosity")
          call fson_get_mpi(r, "density", default_density, density, logfile, &
               trim(rockstr) // "density")
          call fson_get_mpi(r, "specific_heat", default_specific_heat, &
               specific_heat, logfile, trim(rockstr) // "specific_heat")
          perm_size = size(permeability)

          call DMGetStratumSize(dm, rock_type_label_name, ir, num_cells, &
               ierr); CHKERRQ(ierr)
          if (num_cells > 0) then
             call DMGetStratumIS(dm, rock_type_label_name, ir, cell_IS, &
                  ierr); CHKERRQ(ierr)
             call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
             do ic = 1, num_cells
                associate(c => cells(ic))
                  call DMLabelGetValue(ghost_label, c, ghost, ierr)
                  if (ghost < 0) then
                     offset = global_section_offset(section, c, range_start)
                     call rock%assign(rock_array, offset)
                     rock%permeability = 0._dp
                     rock%permeability(1: perm_size) = permeability
                     rock%wet_conductivity = wet_conductivity
                     rock%dry_conductivity = dry_conductivity
                     rock%porosity = porosity
                     rock%density = density
                     rock%specific_heat = specific_heat
                  end if
                end associate
             end do
             call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
             call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
          end if
          deallocate(permeability)

          r => fson_value_next_mpi(r)

       end do

       call rock%destroy()
       call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)

    end if

  end subroutine setup_rock_vector_types

!------------------------------------------------------------------------

  subroutine create_rock_vector(dm, rock_vector, range_start)
    !! Creates and returns rock vector and corresponding range start.

    use dm_utils_module, only: dm_set_data_layout, global_vec_range_start

    DM, intent(in) :: dm !! DM on which to create rock vector
    Vec, intent(out) :: rock_vector !! Output rock vector
    PetscInt, intent(out) :: range_start !! Range start, used for computing offsets
    ! Locals:
    DM :: dm_rock
    PetscInt :: dim, rock_variable_dim(num_rock_variables)
    PetscErrorCode :: ierr

    call DMClone(dm, dm_rock, ierr); CHKERRQ(ierr)

    call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)
    rock_variable_dim = dim
    call dm_set_data_layout(dm_rock, rock_variable_num_components, &
         rock_variable_dim, rock_variable_names)

    call DMCreateGlobalVector(dm_rock, rock_vector, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(rock_vector, "rock", ierr); CHKERRQ(ierr)
    call global_vec_range_start(rock_vector, range_start)
    call DMDestroy(dm_rock, ierr); CHKERRQ(ierr)

  end subroutine create_rock_vector

!------------------------------------------------------------------------

  subroutine setup_rock_vector(json, dm, rock_vector, rock_dict, &
       range_start, ghost_cell, logfile, err)

    !! Sets up rock vector on specified DM from JSON input. If
    !! initialising rock properties using rock types, this routine
    !! also sets up a dictionary of rock types for efficient access by
    !! rock type name.

    use dictionary_module
    use fson
    use fson_mpi_module
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    DM, intent(in out) :: dm
    Vec, intent(in out) :: rock_vector
    type(dictionary_type), intent(in out) :: rock_dict
    PetscInt, intent(out) :: range_start
    PetscInt, allocatable, intent(in) :: ghost_cell(:)
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err

    err = 0

    call create_rock_vector(dm, rock_vector, range_start)
    call set_default_rock_properties(dm, rock_vector, range_start, &
         ghost_cell)

    if (fson_has_mpi(json, "rock")) then

       if (fson_has_mpi(json, "rock.types")) then

          call setup_rock_vector_types(json, dm, rock_vector, rock_dict, &
               range_start, logfile)

       else
          ! other types of rock initialization here- TODO
          ! e.g. read from a rock section in initial conditions HDF5 file
       end if

    end if

  end subroutine setup_rock_vector

!------------------------------------------------------------------------

end module rock_module
