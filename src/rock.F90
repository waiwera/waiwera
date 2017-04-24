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

  use kinds_module
  use relative_permeability_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  type rock_type
     !! Local rock properties.
     private
     PetscReal, pointer, contiguous, public :: permeability(:)   !! Permeability
     PetscReal, pointer, public :: wet_conductivity, dry_conductivity !! Heat conductivities
     PetscReal, pointer, public :: porosity          !! Porosity
     PetscReal, pointer, public :: density           !! Grain density
     PetscReal, pointer, public :: specific_heat     !! Specific heat
     PetscReal, pointer, public :: youngs_modulus    !! Young's modulus for rock mechanics
     PetscReal, pointer, public :: poissons_ratio    !! Poisson's ratio for rock mechanics
     class(relative_permeability_type), pointer, &
          public :: relative_permeability !! Relative permeability functions
     PetscInt, public :: dof !! Number of degrees of freedom
   contains
     private
     procedure, public :: init => rock_init
     procedure, public :: assign => rock_assign
     procedure, public :: assign_relative_permeability => &
          rock_assign_relative_permeability
     procedure, public :: destroy => rock_destroy
     procedure, public :: energy => rock_energy
  end type rock_type

  PetscInt, parameter :: num_rock_variables = 8
  PetscInt, parameter :: max_rock_variable_name_length = 32
  character(max_rock_variable_name_length), parameter, public :: &
       rock_variable_names(num_rock_variables) = &
       [character(max_rock_variable_name_length):: &
       "permeability", "wet_conductivity", "dry_conductivity", "porosity", &
       "density", "specific_heat", "youngs_modulus","poissons_ratio"]
  PetscInt, parameter, public :: &
       rock_variable_num_components(num_rock_variables) = &
       [3, 1, 1, 1, 1, 1, 1, 1]
  PetscInt, parameter, public :: &
       rock_variable_dim(num_rock_variables) = &
       [3, 3, 3, 3, 3, 3, 3, 3]
  PetscInt, parameter, public :: max_rockname_length = 24

  PetscInt, parameter, public :: max_rocktype_label_length = 9
  character(max_rocktype_label_length), parameter, public :: &
       rocktype_label_name = "rock_type"

  ! Default rock properties:
  PetscReal, parameter, public :: default_permeability(3) = [1.e-13_dp, 1.e-13_dp, 1.e-13_dp]
  PetscReal, parameter, public :: default_porosity = 0.1_dp
  PetscReal, parameter, public :: default_density = 2200.0_dp
  PetscReal, parameter, public :: default_specific_heat = 1000._dp
  PetscReal, parameter, public :: default_heat_conductivity = 2.5_dp
  PetscReal, parameter, public :: default_youngs_modulus = 1.0e10_dp
  PetscReal, parameter, public :: default_poissons_ratio = 0.25_dp

  public :: rock_type, setup_rock_vector, setup_rocktype_labels

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
    self%youngs_modulus => data(offset + 8)
    self%poissons_ratio => data(offset + 9)

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

  subroutine rock_destroy(self)
    !! Destroys a rock object (nullifies all pointer components).

    class(rock_type), intent(in out) :: self

    nullify(self%permeability)
    nullify(self%wet_conductivity)
    nullify(self%dry_conductivity)
    nullify(self%porosity)
    nullify(self%density)
    nullify(self%specific_heat)
    nullify(self%relative_permeability)
    nullify(self%youngs_modulus)
    nullify(self%poissons_ratio)

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
          rock%youngs_modulus = default_youngs_modulus
          rock%poissons_ratio = default_poissons_ratio
       end if
    end do

    call rock%destroy()
    call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)

  end subroutine set_default_rock_properties

!------------------------------------------------------------------------

  subroutine setup_rock_vector_types(json, dm, rock_vector, range_start, &
       logfile)
    !! Sets up rock vector on DM from rock types in JSON input.

    use dm_utils_module, only: global_section_offset, global_vec_section
    use fson
    use fson_mpi_module
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    DM, intent(in) :: dm
    Vec, intent(out) :: rock_vector
    PetscInt, intent(in) :: range_start
    type(logfile_type), intent(in out) :: logfile
    ! Locals:
    PetscInt :: num_rocktypes, ir, ic, c, num_cells, offset, ghost
    type(fson_value), pointer :: rocktypes, r
    IS :: rock_IS
    PetscInt, pointer :: rock_cells(:)
    DMLabel :: ghost_label
    type(rock_type) :: rock
    character(max_rockname_length) :: name
    PetscReal :: porosity, density, specific_heat
    PetscReal :: youngs_modulus, poissons_ratio
    PetscReal :: wet_conductivity, dry_conductivity
    PetscReal, allocatable :: permeability(:)
    PetscReal, pointer, contiguous :: rock_array(:)
    PetscSection :: section
    PetscErrorCode :: ierr
    character(len=64) :: rockstr
    character(len=12) :: irstr

    call rock%init()

    call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call global_vec_section(rock_vector, section)

    call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    
    call fson_get_mpi(json, "rock.types", rocktypes)
    num_rocktypes = fson_value_count_mpi(rocktypes, ".")
    do ir = 1, num_rocktypes
       write(irstr, '(i0)') ir - 1
       rockstr = 'rock.types[' // trim(irstr) // '].'
       r => fson_value_get_mpi(rocktypes, ir)
       call fson_get_mpi(r, "name", "", name, logfile, trim(rockstr) // "name")
       call fson_get_mpi(r, "permeability", default_permeability, &
            permeability, logfile, trim(rockstr) // "permeability")
       call fson_get_mpi(r, "wet conductivity", default_heat_conductivity, &
            wet_conductivity, logfile, trim(rockstr) // "wet conductivity")
       call fson_get_mpi(r, "dry conductivity", wet_conductivity, &
            dry_conductivity, logfile, trim(rockstr) // "dry conductivity")
       call fson_get_mpi(r, "porosity", default_porosity, porosity, logfile, &
            trim(rockstr) // "porosity")
       call fson_get_mpi(r, "density", default_density, density, logfile, &
            trim(rockstr) // "density")
       call fson_get_mpi(r, "specific heat", default_specific_heat, &
            specific_heat, logfile, trim(rockstr) // "specific heat")
       call fson_get_mpi(r, "youngs modulus", default_youngs_modulus, youngs_modulus, logfile, &
            trim(rockstr) // "youngs modulus")
       call fson_get_mpi(r, "poissons ratio", default_poissons_ratio, poissons_ratio, logfile, &
            trim(rockstr) // "poissons ratio")
       call DMGetStratumIS(dm, rocktype_label_name, ir, rock_IS, &
            ierr); CHKERRQ(ierr)
       if (rock_IS /= 0) then
          call ISGetIndicesF90(rock_IS, rock_cells, ierr); CHKERRQ(ierr)
          num_cells = size(rock_cells)
          do ic = 1, num_cells
             c = rock_cells(ic)
             call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
             if (ghost < 0) then
                call global_section_offset(section, c, range_start, &
                     offset, ierr); CHKERRQ(ierr)
                call rock%assign(rock_array, offset)
                rock%permeability = permeability
                rock%wet_conductivity = wet_conductivity
                rock%dry_conductivity = dry_conductivity
                rock%porosity = porosity
                rock%density = density
                rock%specific_heat = specific_heat
                rock%youngs_modulus = youngs_modulus
                rock%poissons_ratio = poissons_ratio
             end if
          end do
          call ISRestoreIndicesF90(rock_IS, rock_cells, ierr); CHKERRQ(ierr)
       end if
       call ISDestroy(rock_IS, ierr); CHKERRQ(ierr)
    end do
    call rock%destroy()
    call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)

  end subroutine setup_rock_vector_types

!------------------------------------------------------------------------

  subroutine setup_rock_vector(json, dm, rock_vector, range_start, &
       ghost_cell, logfile)
    !! Sets up rock vector on specified DM from JSON input.

    use dm_utils_module, only: set_dm_data_layout, global_vec_range_start
    use fson
    use fson_mpi_module
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    DM, intent(in) :: dm
    Vec, intent(out) :: rock_vector
    PetscInt, intent(out) :: range_start
    PetscInt, allocatable, intent(in) :: ghost_cell(:)
    type(logfile_type), intent(in out) :: logfile
    ! Locals:
    DM :: dm_rock
    PetscErrorCode :: ierr

    call DMClone(dm, dm_rock, ierr); CHKERRQ(ierr)

    call set_dm_data_layout(dm_rock, rock_variable_num_components, &
         rock_variable_dim, rock_variable_names)

    call DMCreateGlobalVector(dm_rock, rock_vector, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(rock_vector, "rock", ierr); CHKERRQ(ierr)
    call global_vec_range_start(rock_vector, range_start)

    call set_default_rock_properties(dm, rock_vector, range_start, &
         ghost_cell)

    if (fson_has_mpi(json, "rock")) then

       if (fson_has_mpi(json, "rock.types")) then

          call setup_rock_vector_types(json, dm, rock_vector, range_start, &
               logfile)

       else
          ! other types of rock initialization here- TODO
          ! e.g. read from a rock section in initial conditions HDF5 file
       end if

    end if

    call DMDestroy(dm_rock, ierr); CHKERRQ(ierr)

  end subroutine setup_rock_vector

!------------------------------------------------------------------------

  subroutine setup_rocktype_labels(json, dm, logfile)
    !! Sets up rocktype label on a DM. The values of the rock type
    !! label are the indices (1-based) of the rocktypes specified in the 
    !! JSON input file.

    use fson
    use fson_mpi_module
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    DM, intent(in out) :: dm
    type(logfile_type), intent(in out) :: logfile
    ! Locals:
    PetscErrorCode :: ierr
    type(fson_value), pointer :: rocktypes, r
    PetscInt :: num_rocktypes, num_cells, ir, ic, c
    PetscInt, allocatable :: cells(:)
    PetscInt, allocatable :: default_cells(:)
    character(len=64) :: rockstr
    character(len=12) :: irstr

    default_cells = [PetscInt::] ! empty integer array

    if (fson_has_mpi(json, "rock.types")) then
       call fson_get_mpi(json, "rock.types", rocktypes)
       call DMCreateLabel(dm, rocktype_label_name, ierr); CHKERRQ(ierr)
       num_rocktypes = fson_value_count_mpi(rocktypes, ".")
       do ir = 1, num_rocktypes
          write(irstr, '(i0)') ir - 1
          rockstr = 'rock.types[' // trim(irstr) // '].'
          r => fson_value_get_mpi(rocktypes, ir)
          call fson_get_mpi(r, "cells", default_cells, cells, &
               logfile, trim(rockstr) // "cells")
          if (allocated(cells)) then
             num_cells = size(cells)
             do ic = 1, num_cells
                c = cells(ic)
                call DMSetLabelValue(dm, rocktype_label_name, &
                     c, ir, ierr); CHKERRQ(ierr)
             end do
             deallocate(cells)
          end if
       end do
    else
       call logfile%write(LOG_LEVEL_WARN, "input", "no rocktypes")
    end if

  end subroutine setup_rocktype_labels

!------------------------------------------------------------------------

end module rock_module
