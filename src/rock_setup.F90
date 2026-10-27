!   Copyright 2021 University of Auckland.

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

module rock_setup_module
  !! Module for setting up rock properties and controls.

#include <petsc/finclude/petsc.h>

  use petsc
  use fson
  use kinds_module
  use rock_module
  use rock_control_module
  use list_module

  implicit none
  private

  type rock_dict_item_type
     !! Item in rock dictionary.
     private
     type(fson_value), pointer, public :: rock !! JSON input data value
     PetscInt, public :: label_value !! DM rock label value
  end type rock_dict_item_type
  
  public :: rock_dict_item_type
  public :: setup_rock_types, setup_rocks, create_rock_vector, label_rock_cell

contains

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

  subroutine setup_rocks_from_types(json, dm, rock_vector, &
       rock_dict, rock_controls, range_start, eos, logfile, err)
    !! Sets up rock vector on DM and rock controls from rock types in
    !! JSON input.

    use dictionary_module
    use dm_utils_module, only: global_section_offset, global_vec_section
    use zone_label_module
    use fson
    use fson_mpi_module
    use eos_module, only: eos_type
    use logfile_module
    use fson_value_m, only : TYPE_ARRAY, TYPE_OBJECT, TYPE_STRING
    use interpolation_module

    type(fson_value), pointer, intent(in) :: json
    DM, intent(in out) :: dm
    Vec, intent(in out) :: rock_vector
    type(dictionary_type), intent(in out) :: rock_dict
    type(list_type), intent(in out) :: rock_controls
    PetscInt, intent(in) :: range_start
    class(eos_type), intent(in) :: eos
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: num_rocktypes, ir, ic, offset, perm_size
    DMLabel :: ghost_label
    type(fson_value), pointer :: rocktypes, r
    PetscInt, allocatable :: rock_cell(:)
    type(rock_type) :: rock
    character(max_rockname_length) :: name
    PetscReal :: porosity, density, specific_heat
    PetscReal :: wet_conductivity, dry_conductivity
    PetscReal, allocatable :: permeability(:)
    PetscReal, pointer, contiguous :: rock_array(:)
    PetscSection :: section
    PetscBool :: permeability_init, porosity_init, specific_heat_init
    PetscErrorCode :: ierr
    character(len=64) :: rockstr
    character(len=12) :: irstr

    err = 0

    if (fson_has_mpi(json, "rock.types")) then

       call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
       call rock%init()

       call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
       call global_vec_section(rock_vector, section)

       call fson_get_mpi(json, "rock.types", rocktypes)
       num_rocktypes = fson_value_count_mpi(rocktypes, ".")
       r => fson_value_children_mpi(rocktypes)

       do ir = 1, num_rocktypes

          write(irstr, '(i0)') ir - 1
          rockstr = 'rock.types[' // trim(irstr) // '].'
          call fson_get_mpi(r, "name", "", name, logfile, trim(rockstr) // "name")

          call get_rock_cells(dm, rock_type_label_name, ir, rock_cell)

          call setup_permeability(r, rock_cell, permeability, &
               permeability_init, err)
          if (err == 0) then

             call setup_porosity(r, rock_cell, porosity, porosity_init)

             call fson_get_mpi(r, "wet_conductivity", default_heat_conductivity, &
                  wet_conductivity, logfile, trim(rockstr) // "wet_conductivity")
             call fson_get_mpi(r, "dry_conductivity", wet_conductivity, &
                  dry_conductivity, logfile, trim(rockstr) // "dry_conductivity")
             call fson_get_mpi(r, "density", default_density, density, logfile, &
                  trim(rockstr) // "density")
             call setup_specific_heat(r, rock_cell, specific_heat, specific_heat_init)
             if (permeability_init) then
                perm_size = size(permeability)
             end if

             do ic = 1, size(rock_cell)
                associate(c => rock_cell(ic))
                  offset = global_section_offset(section, c, range_start)
                  call rock%assign(rock_array, offset)
                  if (permeability_init) then
                     rock%permeability = 0._dp
                     rock%permeability(1: perm_size) = permeability
                  end if
                  rock%wet_conductivity = wet_conductivity
                  rock%dry_conductivity = dry_conductivity
                  if (porosity_init) then
                     rock%porosity = porosity
                  end if
                  rock%density = density
                  if (specific_heat_init) then
                     rock%specific_heat = specific_heat
                  end if
                end associate
             end do

          end if

          if (allocated(permeability)) deallocate(permeability)
          r => fson_value_next_mpi(r)

       end do

       call rock%destroy()
       call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)

    end if

  contains

    subroutine get_rock_cells(dm, rock_type_label_name, ir, rock_cell)
      !! Returns array of cell indices on local process with specified
      !! rock type name.

      DM, intent(in out) :: dm
      character(*), intent(in) :: rock_type_label_name
      PetscInt, intent(in) :: ir
      PetscInt, allocatable, intent(out) :: rock_cell(:)
      ! Locals:
      IS :: cell_IS
      PetscInt, pointer :: cells(:), ghost_cell(:)
      PetscInt :: num_cells, i
      PetscErrorCode :: ierr

      call DMGetStratumSize(dm, rock_type_label_name, ir, num_cells, &
           ierr); CHKERRQ(ierr)
      if (num_cells > 0) then
         call DMGetStratumIS(dm, rock_type_label_name, ir, cell_IS, &
              ierr); CHKERRQ(ierr)
         call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
         allocate(ghost_cell(num_cells))
         do i = 1, num_cells
            call DMLabelGetValue(ghost_label, cells(i), ghost_cell(i), &
                 ierr); CHKERRQ(ierr)
         end do
         rock_cell = pack(cells, ghost_cell < 0)
         call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
         call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
      else
         allocate(rock_cell(0))
      end if

    end subroutine get_rock_cells

!........................................................................

    subroutine setup_permeability(rock_json, rock_cell, permeability, &
         permeability_init, err)
      !! Sets up rock permeability and any associated controls.

      type(fson_value), pointer, intent(in) :: rock_json
      PetscInt, intent(in) :: rock_cell(:)
      PetscReal, allocatable, intent(out) :: permeability(:)
      PetscBool, intent(out) :: permeability_init
      PetscErrorCode, intent(in out) :: err
      ! Locals:
      PetscInt :: permeability_type, permeability_rank
      PetscInt :: permeability_value_type, dim
      PetscReal :: permeability_scalar
      type(fson_value), pointer :: perm

      permeability_init = PETSC_TRUE
       call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)

       permeability_type = fson_type_mpi(rock_json, "permeability")
      select case (permeability_type)
      case (TYPE_ARRAY)
         permeability_rank = fson_mpi_array_rank(rock_json, "permeability")
         select case (permeability_rank)
         case (1)
            call fson_get_mpi(rock_json, "permeability", default_permeability, &
                 permeability, logfile, trim(rockstr) // "permeability")
         case (2) ! table of permeabilities vs. time
            call setup_permeability_time_rock_control(rock_json, "permeability", &
                 rock_cell)
            permeability_init = PETSC_FALSE
         end select
      case (TYPE_OBJECT)
         if (fson_has_mpi(rock_json, "permeability.time")) then
            call fson_get_mpi(rock_json, "permeability", perm)
            call setup_permeability_time_rock_control(perm, "time", rock_cell)
            permeability_init = PETSC_FALSE
         else if (fson_has_mpi(rock_json, "permeability.temperature")) then
            permeability_value_type = fson_type_mpi(rock_json, "permeability.value")
            select case (permeability_value_type)
            case (TYPE_ARRAY)
               call fson_get_mpi(rock_json, "permeability.value", default_permeability, &
                    permeability, logfile, trim(rockstr) // "permeability.value")
            case default ! real or null
               call fson_get_mpi(rock_json, "permeability.value", &
                    default_permeability_scalar, permeability_scalar, &
                    logfile, trim(rockstr) // "permeability.value")
               allocate(permeability(dim))
               permeability = permeability_scalar
            end select
            call setup_permeability_temperature_rock_control(rock_json, rock_cell, &
                 permeability, eos, err)
            permeability_init = PETSC_FALSE
         else if (fson_has_mpi(rock_json, "permeability.value")) then
            permeability_type = fson_type_mpi(rock_json, "permeability.value")
            select case (permeability_type)
            case (TYPE_ARRAY)
               call fson_get_mpi(rock_json, "permeability.value", default_permeability, &
                    permeability, logfile, trim(rockstr) // "permeability.value")
            case default
               call fson_get_mpi(rock_json, "permeability.value", &
                    default_permeability_scalar, permeability_scalar, &
                    logfile, trim(rockstr) // "permeability.value")
               allocate(permeability(dim))
               permeability = permeability_scalar
            end select
         end if
      case default ! real or null
         call fson_get_mpi(rock_json, "permeability", default_permeability_scalar, &
              permeability_scalar, logfile, trim(rockstr) // "permeability")
         allocate(permeability(dim))
         permeability = permeability_scalar
      end select

    end subroutine setup_permeability

!........................................................................

    subroutine setup_porosity(rock_json, rock_cell, porosity, porosity_init)
      !! Sets up rock porosity and any associated controls.

      type(fson_value), pointer, intent(in) :: rock_json
      PetscInt, intent(in) :: rock_cell(:)
      PetscReal, intent(out) :: porosity
      PetscBool, intent(out) :: porosity_init
      ! Locals:
      PetscInt :: porosity_type

      porosity_init = PETSC_TRUE

      porosity_type = fson_type_mpi(rock_json, "porosity")
      select case (porosity_type)
      case (TYPE_ARRAY)
         call setup_porosity_time_rock_control(rock_json, rock_cell)
         porosity_init = PETSC_FALSE
      case default
         call fson_get_mpi(rock_json, "porosity", default_porosity, &
              porosity, logfile, trim(rockstr) // "porosity")
      end select

    end subroutine setup_porosity

!........................................................................

    subroutine setup_specific_heat(rock_json, rock_cell, specific_heat, &
         specific_heat_init)
      !! Sets up rock specific heat and any associated controls.

      type(fson_value), pointer, intent(in) :: rock_json
      PetscInt, intent(in) :: rock_cell(:)
      PetscReal, intent(out) :: specific_heat
      PetscBool, intent(out) :: specific_heat_init
      ! Locals:
      PetscInt :: specific_heat_type, interpolation_type
      PetscReal, allocatable :: specific_heat_table(:,:)
      type(fson_value), pointer :: sph
      type(temperature_dependent_specific_heat_rock_control_type), pointer :: control
      character(max_interpolation_str_length) :: interpolation_str
      character(4), parameter :: default_specific_heat_interpolation_str = &
           "step"

      specific_heat_init = PETSC_TRUE

      specific_heat_type = fson_type_mpi(rock_json, "specific_heat")
      select case (specific_heat_type)
      case (TYPE_OBJECT)
         call fson_get_mpi(rock_json, "specific_heat", sph)
         call fson_get_mpi(sph, "temperature", val = specific_heat_table)
         call fson_get_mpi(sph, "interpolation", &
              default_specific_heat_interpolation_str, interpolation_str)
         interpolation_type = interpolation_type_from_str(interpolation_str)
         specific_heat_init = PETSC_FALSE

         if (size(rock_cell) > 0) then
            allocate(control)
            call control%init(specific_heat_table, rock_cell, interpolation_type, eos)
            call rock_controls%append(control)
         end if
         deallocate(specific_heat_table)

      case default ! real or null
         call fson_get_mpi(rock_json, "specific_heat", default_specific_heat, &
              specific_heat, logfile, trim(rockstr) // "specific_heat")
      end select

    end subroutine setup_specific_heat

!........................................................................

      subroutine setup_permeability_time_rock_control(rock_json, key, rock_cell)
      !! Sets up time-dependent permeability rock control from JSON input.

      type(fson_value), pointer, intent(in) :: rock_json
      character(*), intent(in) :: key
      PetscInt, intent(in) :: rock_cell(:)
      ! Locals:
      PetscReal, allocatable :: permeability_table(:,:)
      type(permeability_table_rock_control_type), pointer :: control
      character(max_interpolation_str_length) :: interpolation_str
      PetscInt :: interpolation_type

      call fson_get_mpi(rock_json, key, val = permeability_table)
      call fson_get_mpi(rock_json, "interpolation", &
           default_interpolation_str, interpolation_str)
      interpolation_type = interpolation_type_from_str(interpolation_str)

      if (size(rock_cell) > 0) then
         allocate(control)
         call control%init(permeability_table, rock_cell, interpolation_type)
         call rock_controls%append(control)
      end if
      deallocate(permeability_table)

    end subroutine setup_permeability_time_rock_control

!........................................................................

    subroutine setup_permeability_temperature_rock_control(rock_json, rock_cell, &
         permeability, eos, err)
      !! Sets up temperature-dependent permeability rock control from JSON input.

      use utils_module, only: str_to_lower

      type(fson_value), pointer, intent(in) :: rock_json
      PetscInt, intent(in) :: rock_cell(:)
      PetscReal, intent(in) :: permeability(:)
      class(eos_type), intent(in) :: eos
      PetscErrorCode, intent(in out) :: err

      ! Locals:
      PetscReal, allocatable :: permeability_table(:,:)
      class(temperature_dependent_permeability_rock_control_type), pointer :: control
      character(max_interpolation_str_length) :: interpolation_str
      character(32) :: permeability_control_name
      PetscInt :: control_type, interpolation_type

      control_type = fson_type_mpi(rock_json, "permeability.temperature")
      select case (control_type)
      case (TYPE_ARRAY)
         call fson_get_mpi(rock_json, "permeability.temperature", &
              val = permeability_table)
         call fson_get_mpi(rock_json, "interpolation", &
              default_interpolation_str, interpolation_str)
         interpolation_type = interpolation_type_from_str(interpolation_str)
         allocate(temperature_dependent_permeability_rock_control_type :: control)
      case (TYPE_STRING)
         call fson_get_mpi(rock_json, "permeability.temperature", &
              val = permeability_control_name)
         select case (str_to_lower(permeability_control_name))
         case ('hayba-ingebritsen')
            associate(nk => size(permeability))
              allocate(permeability_table(1, nk))
            end associate
            permeability_table(1, :) = permeability
            allocate(hayba_ingebritsen_permeability_rock_control_type :: control)
         case default
            if (present(logfile)) then
               if (logfile%active) then
                  call logfile%write(LOG_LEVEL_WARN, 'input', 'unrecognised', &
                       str_key = trim(rockstr) // ".permeability.temperature", &
                       str_value = permeability_control_name)
               end if
            end if
            err = 1
         end select
      end select

      if (size(rock_cell) > 0) then
         call control%init(permeability_table, rock_cell, interpolation_type, eos)
         call rock_controls%append(control)
      else
         deallocate(control)
      end if
      deallocate(permeability_table)

    end subroutine setup_permeability_temperature_rock_control

!........................................................................

    subroutine setup_porosity_time_rock_control(rock_json, rock_cell)
      !! Sets up time-dependent porosity rock control from JSON input.

      type(fson_value), pointer, intent(in) :: rock_json
      PetscInt, intent(in) :: rock_cell(:)
      ! Locals:
      PetscReal, allocatable :: porosity_table(:,:)
      type(porosity_table_rock_control_type), pointer :: control
      character(max_interpolation_str_length) :: interpolation_str
      PetscInt :: interpolation_type

      call fson_get_mpi(rock_json, "porosity", val = porosity_table)
      call fson_get_mpi(rock_json, "interpolation", &
           default_interpolation_str, interpolation_str)
      interpolation_type = interpolation_type_from_str(interpolation_str)
      if (size(rock_cell) > 0) then
         allocate(control)
         call control%init(porosity_table, rock_cell, interpolation_type)
         call rock_controls%append(control)
      end if
      deallocate(porosity_table)

    end subroutine setup_porosity_time_rock_control

  end subroutine setup_rocks_from_types

!------------------------------------------------------------------------

  subroutine create_rock_vector(dm, rock_vector, range_start)
    !! Creates and returns rock vector and corresponding range start.

    use rock_module, only: num_rock_variables
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

  subroutine setup_rocks(json, dm, rock_vector, rock_dict, &
       rock_controls, range_start, ghost_cell, eos, logfile, err)

    !! Sets up rock vector on specified DM and rock controls from JSON
    !! input. If initialising rock properties using rock types, this
    !! routine also sets up a dictionary of rock types for efficient
    !! access by rock type name.

    use dictionary_module
    use fson
    use fson_mpi_module
    use eos_module, only: eos_type
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    DM, intent(in out) :: dm
    Vec, intent(in out) :: rock_vector
    type(dictionary_type), intent(in out) :: rock_dict
    type(list_type), intent(in out) :: rock_controls
    PetscInt, intent(out) :: range_start
    PetscInt, allocatable, intent(in) :: ghost_cell(:)
    class(eos_type), intent(in) :: eos
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err

    err = 0

    call create_rock_vector(dm, rock_vector, range_start)
    call set_default_rock_properties(dm, rock_vector, range_start, &
         ghost_cell)
    call rock_controls%init(owner = PETSC_TRUE)

    if (fson_has_mpi(json, "rock")) then

       if (fson_has_mpi(json, "rock.types")) then

          call setup_rocks_from_types(json, dm, rock_vector, &
               rock_dict, rock_controls, range_start, eos, logfile, err)

       else
          ! other types of rock initialization here- TODO
          ! e.g. read from a rock section in initial conditions HDF5 file
       end if

    end if

  end subroutine setup_rocks

!------------------------------------------------------------------------
  
end module rock_setup_module
  
