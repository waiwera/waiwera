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

  implicit none
  private

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
  PetscReal, parameter, public :: default_permeability(3) = [1.e-13_dp, 1.e-13_dp, 1.e-13_dp]
  PetscReal, parameter, public :: default_porosity = 0.1_dp
  PetscReal, parameter, public :: default_density = 2200.0_dp
  PetscReal, parameter, public :: default_specific_heat = 1000._dp
  PetscReal, parameter, public :: default_heat_conductivity = 2.5_dp

  public :: rock_type, setup_rock_vector

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
          call global_section_offset(rock_section, c, range_start, &
               offset, ierr); CHKERRQ(ierr)
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

  subroutine setup_rock_vector_types(json, dm, ao, rock_vector, rock_dict, &
       range_start, logfile, err)
    !! Sets up rock vector on DM from rock types in JSON input.

    use dictionary_module
    use cell_order_module, only: cell_order_label_name
    use dm_utils_module, only: global_section_offset, global_vec_section, &
         natural_to_local_cell_index
    use zone_label_module
    use fson
    use fson_mpi_module
    use logfile_module
    use fson_value_m, only : TYPE_STRING, TYPE_ARRAY

    type(fson_value), pointer, intent(in) :: json
    DM, intent(in) :: dm
    AO, intent(in) :: ao
    Vec, intent(out) :: rock_vector
    type(dictionary_type), intent(in out) :: rock_dict
    PetscInt, intent(in) :: range_start
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: num_rocktypes, ic, iz, c
    PetscInt, target :: ir
    PetscInt :: perm_size, num_zone_cells
    type(fson_value), pointer :: rocktypes, r
    PetscInt :: start_cell, end_cell
    PetscInt, allocatable :: natural_cell_indices(:), local_cell_indices(:)
    character(max_zone_name_length), allocatable :: zones(:)
    character(:), allocatable :: label_name
    PetscInt, pointer :: cells(:)
    DMLabel :: ghost_label
    type(rock_type) :: rock
    character(max_rockname_length) :: name
    PetscBool :: has_label
    IS :: cell_IS
    PetscReal :: porosity, density, specific_heat
    PetscReal :: wet_conductivity, dry_conductivity
    PetscReal, allocatable :: permeability(:)
    PetscReal, pointer, contiguous :: rock_array(:)
    PetscSection :: section
    ISLocalToGlobalMapping :: l2g
    PetscErrorCode :: ierr
    character(len=64) :: rockstr
    character(len=12) :: irstr
    PetscInt :: zones_type

    err = 0

    if (fson_has_mpi(json, "rock.types")) then

       call rock%init()

       call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
       CHKERRQ(ierr)

       call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
       call global_vec_section(rock_vector, section)

       call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
       call DMGetLocalToGlobalMapping(dm, l2g, ierr); CHKERRQ(ierr)

       call fson_get_mpi(json, "rock.types", rocktypes)
       num_rocktypes = fson_value_count_mpi(rocktypes, ".")
       r => fson_value_children_mpi(rocktypes)

       do ir = 1, num_rocktypes

          write(irstr, '(i0)') ir - 1
          rockstr = 'rock.types[' // trim(irstr) // '].'
          call fson_get_mpi(r, "name", "", name, logfile, trim(rockstr) // "name")
          if (name /= "") then
             call rock_dict%add(name, r)
          end if
          call fson_get_mpi(r, "permeability", default_permeability, &
               permeability, logfile, trim(rockstr) // "permeability")
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

          if (fson_has_mpi(r, "cells")) then
             call fson_get_mpi(r, "cells", val = natural_cell_indices)
             if (allocated(natural_cell_indices)) then
                associate(num_cells => size(natural_cell_indices))
                  allocate(local_cell_indices(num_cells))
                  local_cell_indices = natural_to_local_cell_index(ao, &
                       l2g, natural_cell_indices)
                  do ic = 1, num_cells
                     associate(c => local_cell_indices(ic))
                       if (c >= 0) call assign_rock_parameters(c)
                     end associate
                  end do
                end associate
                deallocate(natural_cell_indices, local_cell_indices)
             end if
          end if

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
                           call assign_rock_parameters(cells(c))
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

          r => fson_value_next_mpi(r)

       end do

       call rock%destroy()
       call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)

    else
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_WARN, "input", "no rocktypes")
       end if
    end if

  contains

    subroutine assign_rock_parameters(p)
      !! Assigns rock parameters for cell at mesh point p.

      PetscInt, intent(in) :: p
      ! Locals:
      PetscInt :: offset, ghost

      if ((p >= start_cell) .and. (p < end_cell)) then
         call DMLabelGetValue(ghost_label, p, ghost, ierr); CHKERRQ(ierr)
         if (ghost < 0) then
            call global_section_offset(section, p, range_start, &
                 offset, ierr); CHKERRQ(ierr)
            call rock%assign(rock_array, offset)
            rock%permeability = 0._dp
            rock%permeability(1: perm_size) = permeability
            rock%wet_conductivity = wet_conductivity
            rock%dry_conductivity = dry_conductivity
            rock%porosity = porosity
            rock%density = density
            rock%specific_heat = specific_heat
         end if
      end if

    end subroutine assign_rock_parameters

  end subroutine setup_rock_vector_types

!------------------------------------------------------------------------

  subroutine setup_rock_vector(json, dm, ao, rock_vector, rock_dict, &
       range_start, ghost_cell, logfile, err)

    !! Sets up rock vector on specified DM from JSON input. If
    !! initialising rock properties using rock types, this routine
    !! also sets up a dictionary of rock types for efficient access by
    !! rock type name.

    use dictionary_module
    use dm_utils_module, only: set_dm_data_layout, global_vec_range_start
    use fson
    use fson_mpi_module
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    DM, intent(in) :: dm
    AO, intent(in) :: ao
    Vec, intent(out) :: rock_vector
    type(dictionary_type), intent(in out) :: rock_dict
    PetscInt, intent(out) :: range_start
    PetscInt, allocatable, intent(in) :: ghost_cell(:)
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    DM :: dm_rock
    PetscInt :: dim, rock_variable_dim(num_rock_variables)
    PetscErrorCode :: ierr

    err = 0

    call DMClone(dm, dm_rock, ierr); CHKERRQ(ierr)

    call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)
    rock_variable_dim = dim
    call set_dm_data_layout(dm_rock, rock_variable_num_components, &
         rock_variable_dim, rock_variable_names)

    call DMCreateGlobalVector(dm_rock, rock_vector, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(rock_vector, "rock", ierr); CHKERRQ(ierr)
    call global_vec_range_start(rock_vector, range_start)

    call set_default_rock_properties(dm, rock_vector, range_start, &
         ghost_cell)

    if (fson_has_mpi(json, "rock")) then

       if (fson_has_mpi(json, "rock.types")) then

          call setup_rock_vector_types(json, dm, ao, rock_vector, rock_dict, &
               range_start, logfile, err)

       else
          ! other types of rock initialization here- TODO
          ! e.g. read from a rock section in initial conditions HDF5 file
       end if

    end if

    call DMDestroy(dm_rock, ierr); CHKERRQ(ierr)

  end subroutine setup_rock_vector

!------------------------------------------------------------------------

end module rock_module
