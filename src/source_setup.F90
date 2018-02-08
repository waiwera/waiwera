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

module source_setup_module
  !! Module for setting up sources and sinks.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fson
  use fson_value_m, only: TYPE_INTEGER, TYPE_REAL, TYPE_ARRAY, &
       TYPE_STRING, TYPE_OBJECT
  use fson_mpi_module
  use list_module
  use logfile_module
  use eos_module
  use thermodynamics_module, only: thermodynamics_type
  use source_module
  use source_control_module

  implicit none
  private

  public :: setup_sources

contains

!------------------------------------------------------------------------

  subroutine setup_sources(json, dm, ao, eos, thermo, start_time, &
       fluid_vector, fluid_range_start, source_vector, source_range_start, &
       source_controls, logfile, err)
    !! Sets up sinks / sources and source controls.

    use dm_utils_module, only: global_vec_section, global_vec_range_start, &
         global_section_offset, create_path_dm, set_dm_data_layout

    type source_spec_type
       !! Specification for a source.
       PetscInt :: spec_index
       type(fson_value), pointer :: json
       PetscInt, allocatable :: cell_natural_index(:), cell_local_index(:)
    end type source_spec_type

    type(fson_value), pointer, intent(in) :: json !! JSON file object
    DM, intent(in) :: dm !! Mesh DM
    AO, intent(in) :: ao !! Application ordering for natural to global cell indexing
    class(eos_type), intent(in) :: eos !! Equation of state
    class(thermodynamics_type), intent(in) :: thermo !! Thermodynamics formulation
    PetscReal, intent(in) :: start_time
    Vec, intent(in) :: fluid_vector !! Fluid vector
    PetscInt, intent(in) :: fluid_range_start !! Range start for global fluid vector
    Vec, intent(out) :: source_vector !! Source vector
    PetscInt, intent(out) :: source_range_start !! Range start for global source vector
    type(list_type), intent(in out) :: source_controls !! List of source controls
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    type(list_type) :: source_list
    type(source_spec_type), pointer :: spec
    DM :: dm_source
    type(source_type) :: source
    PetscInt :: i, j, spec_index, local_source_index
    PetscInt :: total_num_cell_sources
    type(fson_value), pointer :: sources_json, source_json
    PetscInt :: num_sources, num_cells
    PetscInt, allocatable :: cell_natural_index(:), cell_local_index(:)
    PetscReal, pointer, contiguous :: fluid_data(:), source_data(:)
    PetscSection :: fluid_section, source_section
    PetscInt, allocatable :: num_field_components(:), field_dim(:)
    PetscInt, parameter :: max_field_name_length = 40
    character(max_field_name_length), allocatable :: field_names(:)
    PetscErrorCode :: ierr

    total_num_cell_sources = 0
    call source%init(eos)
    call source_list%init(owner = PETSC_TRUE)

    if (fson_has_mpi(json, "source")) then

       call fson_get_mpi(json, "source", sources_json)
       num_sources = fson_value_count_mpi(sources_json, ".")
       source_json => fson_value_children_mpi(sources_json)
       do spec_index = 0, num_sources - 1
          call get_cells(source_json, dm, ao, cell_natural_index, &
               cell_local_index, num_cells)
          if (num_cells > 0) then
             allocate(spec)
             spec%spec_index = spec_index
             spec%cell_natural_index = cell_natural_index
             spec%cell_local_index = cell_local_index
             spec%json => source_json
             call source_list%append(spec)
             total_num_cell_sources = total_num_cell_sources + num_cells
          end if
          source_json => fson_value_next_mpi(source_json)
       end do

    end if

    call create_path_dm(total_num_cell_sources, dm_source)

    allocate(num_field_components(source%dof), field_dim(source%dof), &
         field_names(source%dof))
    num_field_components = 1
    field_dim = 1
    field_names(1: num_source_variables - 1) = &
         source_variable_names(1: num_source_variables - 1) ! scalar fields
    i = num_source_variables
    ! array fields (flow):
    do j = 1, eos%num_components
       field_names(i) = trim(eos%component_names(j)) // '_' // &
            trim(source_variable_names(num_source_variables))
       i = i + 1
    end do
    if (.not. eos%isothermal) field_names(i) = 'heat_' // &
         trim(source_variable_names(num_source_variables))
    call set_dm_data_layout(dm_source, num_field_components, field_dim, &
         field_names)
    call DMCreateGlobalVector(dm_source, source_vector, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(source_vector, "source", ierr); CHKERRQ(ierr)
    call global_vec_range_start(source_vector, source_range_start)

    call source_controls%init(owner = PETSC_TRUE)
    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayReadF90(fluid_vector, fluid_data, ierr); CHKERRQ(ierr)
    call global_vec_section(source_vector, source_section)
    call VecGetArrayF90(source_vector, source_data, ierr); CHKERRQ(ierr)

    if (fson_has_mpi(json, "source")) then
       local_source_index = 0
       call source_list%traverse(source_setup_iterator)
    else
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_INFO, "input", "no_sources")
       end if
    end if

    call VecRestoreArrayF90(source_vector, source_data, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(fluid_vector, fluid_data, ierr)
    CHKERRQ(ierr)
    call source%destroy()
    call source_list%destroy()

  contains

    subroutine source_setup_iterator(node, stopped)
      !! Iterator for setting up cell sources for a source specification.

      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      character(len=64) :: srcstr
      character(len=12) :: istr
      PetscBool :: can_inject
      PetscInt :: injection_component, production_component
      PetscReal :: initial_rate, initial_enthalpy
      PetscInt :: i, s, source_offset
      PetscInt, allocatable :: source_indices(:)

      select type (spec => node%data)
      type is (source_spec_type)
         write(istr, '(i0)') spec%spec_index
         srcstr = 'source[' // trim(istr) // '].'
         call get_components(spec%json, eos, &
              injection_component, production_component, logfile)
         call get_initial_rate(spec%json, initial_rate, can_inject)
         call get_initial_enthalpy(spec%json, eos, can_inject, &
              injection_component, initial_enthalpy)
         associate(num_sources => size(spec%cell_natural_index))
           allocate(source_indices(num_sources))
           do i = 1, num_sources
              s = spec%spec_index
              call global_section_offset(source_section, s, &
                   source_range_start, source_offset, ierr); CHKERRQ(ierr)
              call source%assign(source_data, source_offset)
              call source%setup(s, spec%cell_natural_index(i), &
                   spec%cell_local_index(i), initial_rate, initial_enthalpy, &
                   injection_component, production_component)
              source_indices(i) = local_source_index
              local_source_index = local_source_index + 1
           end do
         end associate
         call setup_inline_source_controls(spec%json, eos, thermo, &
              start_time, source_data, source_section, source_range_start, &
              fluid_data, fluid_section, fluid_range_start, srcstr, &
              source_indices, source_controls, logfile, err)
         deallocate(source_indices)
         stopped = (err > 0)
      end select

    end subroutine source_setup_iterator

  end subroutine setup_sources

!------------------------------------------------------------------------

  PetscInt function get_component(source_json, eos, name) result(component)
    !! Gets a specified component type for the source.

    type(fson_value), pointer, intent(in) :: source_json
    class(eos_type), intent(in) :: eos
    character(len=*), intent(in) :: name
    ! Locals:
    character(max_component_name_length) :: component_str
    PetscInt :: component_type

    component_type = fson_type_mpi(source_json, name)
    if (component_type == TYPE_INTEGER) then
       call fson_get_mpi(source_json, name, val = component)
    else if (component_type == TYPE_STRING) then
       call fson_get_mpi(source_json, name, val = component_str)
       component = eos%component_index(component_str)
    end if

  end function get_component

!------------------------------------------------------------------------

  subroutine get_components(source_json, eos, &
       injection_component, production_component, logfile)
    !! Gets injection and production components for the source.

    type(fson_value), pointer, intent(in) :: source_json
    class(eos_type), intent(in) :: eos
    PetscInt, intent(out) :: injection_component
    PetscInt, intent(out) :: production_component
    type(logfile_type), intent(in out), optional :: logfile

    if (fson_has_mpi(source_json, "component")) then
       injection_component = get_component(source_json, eos, &
            "component")
    else
       injection_component = default_source_component
    end if

    if (fson_has_mpi(source_json, "production_component")) then
       production_component = get_component(source_json, eos, &
            "production_component")
    else
       associate(np => eos%num_primary_variables)
         if ((.not. eos%isothermal) .and. (injection_component == np)) then
            production_component = injection_component
         else
            production_component = default_source_production_component
         end if
       end associate
    end if

  end subroutine get_components

!------------------------------------------------------------------------

  subroutine get_initial_rate(source_json, initial_rate, can_inject)
    !! Gets initial flow rate and whether the source can inject. These
    !! can only be determined for constant-rate sources. For other
    !! source types, a default initial rate is assigned and it is
    !! assumed that the source may inject.

    type(fson_value), pointer, intent(in) :: source_json
    PetscReal, intent(out) :: initial_rate
    PetscBool, intent(out) :: can_inject
    ! Locals:
    PetscInt :: rate_type

    if (fson_has_mpi(source_json, "rate")) then

       rate_type = fson_type_mpi(source_json, "rate")

       select case (rate_type)
       case (TYPE_REAL, TYPE_INTEGER)
          call fson_get_mpi(source_json, "rate", val = initial_rate)
          can_inject = (initial_rate > 0._dp)
       case (TYPE_ARRAY)
          initial_rate = default_source_rate
          can_inject = PETSC_TRUE
       end select

    else
       initial_rate = default_source_rate
       can_inject = PETSC_TRUE
    end if

  end subroutine get_initial_rate

!------------------------------------------------------------------------

  subroutine get_initial_enthalpy(source_json, eos, can_inject, &
       injection_component, enthalpy)
    !! Gets initial injection enthalpy for the source, if needed.

    type(fson_value), pointer, intent(in) :: source_json
    class(eos_type), intent(in) :: eos
    PetscBool, intent(in) :: can_inject
    PetscInt, intent(in) :: injection_component
    PetscReal, intent(out) :: enthalpy
    ! Locals:
    PetscInt :: enthalpy_type

    associate(np => eos%num_primary_variables)

      if (can_inject .and. (injection_component < np)) then

         if (fson_has_mpi(source_json, "enthalpy")) then

            enthalpy_type = fson_type_mpi(source_json, "enthalpy")

            select case (enthalpy_type)
            case(TYPE_REAL, TYPE_INTEGER)
               call fson_get_mpi(source_json, "enthalpy", &
                    val = enthalpy)
            case default
               enthalpy = default_source_injection_enthalpy
            end select

         else
            enthalpy = default_source_injection_enthalpy
         end if

      else
         enthalpy = 0._dp
      end if

    end associate

  end subroutine get_initial_enthalpy

!------------------------------------------------------------------------

  subroutine get_cells(source_json, dm, ao, cell_natural_index, &
       cell_local_index, num_cells)
    !! Gets arrays of cell natural and local indices for the
    !! source. The "cell" key can be used to specify a single
    !! cell. The "cells" key can be used to specify either a single
    !! cell or an array of cells. If both are present, the "cells" key
    !! is used.

    use dm_utils_module, only: natural_to_local_cell_index

    type(fson_value), pointer, intent(in) :: source_json
    DM, intent(in) :: dm
    AO, intent(in) :: ao
    PetscInt, allocatable, intent(out) :: cell_natural_index(:)
    PetscInt, allocatable, intent(out) :: cell_local_index(:)
    PetscInt, intent(out) :: num_cells
    ! Locals:
    ISLocalToGlobalMapping :: l2g
    PetscInt :: cell_type, cell
    PetscInt :: i, ghost
    DMLabel :: ghost_label
    PetscBool, allocatable :: filter(:)
    PetscErrorCode :: ierr

    if (fson_has_mpi(source_json, "cell")) then
       call fson_get_mpi(source_json, "cell", val = cell)
       cell_natural_index = [cell]
    end if

    if (fson_has_mpi(source_json, "cells")) then

       cell_type = fson_type_mpi(source_json, "cells")

       if (cell_type == TYPE_INTEGER) then
          call fson_get_mpi(source_json, "cells", val = cell)
          cell_natural_index = [cell]
       else if (cell_type == TYPE_ARRAY) then
          call fson_get_mpi(source_json, "cells", &
               val = cell_natural_index)
       end if

    end if

    if (allocated(cell_natural_index)) then

       call DMGetLocalToGlobalMapping(dm, l2g, ierr); CHKERRQ(ierr)
       cell_local_index = natural_to_local_cell_index(ao, l2g, &
            cell_natural_index)

       ! Flag ghost cells:
       call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
       do i = 1, size(cell_local_index)
          associate(c => cell_local_index(i))
            if (c >= 0) then
               call DMLabelGetValue(ghost_label, c, ghost, ierr)
               CHKERRQ(ierr)
               if (ghost >= 0) c = -1
            end if
          end associate
       end do

       ! Filter out off-process and ghost cells:
       filter = cell_local_index >= 0
       cell_local_index = pack(cell_local_index, filter)
       cell_natural_index = pack(cell_natural_index, filter)
       num_cells = size(cell_natural_index)

    else
       num_cells = 0
    end if

  end subroutine get_cells

!------------------------------------------------------------------------

  subroutine setup_inline_source_controls(source_json, eos, thermo, &
       start_time, source_data, source_section, source_range_start, &
       fluid_data, fluid_section, fluid_range_start, &
       srcstr, source_indices, source_controls, logfile, err)
    !! Sets up any 'inline' source controls for the source specification,
    !! i.e. controls defined implicitly in the specification.

    use interpolation_module, only: interpolation_type_from_str, &
         averaging_type_from_str, max_interpolation_str_length, &
         max_averaging_str_length, default_interpolation_str, &
         default_averaging_str

    type(fson_value), pointer, intent(in) :: source_json
    class(eos_type), intent(in) :: eos
    class(thermodynamics_type), intent(in) :: thermo
    PetscReal, intent(in) :: start_time
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    character(len = *), intent(in) :: srcstr
    PetscInt, intent(in) :: source_indices(:)
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: interpolation_type, averaging_type
    character(max_interpolation_str_length) :: interpolation_str
    character(max_averaging_str_length) :: averaging_str

    call fson_get_mpi(source_json, "interpolation", &
         default_interpolation_str, interpolation_str)
    interpolation_type = interpolation_type_from_str(interpolation_str)

    call fson_get_mpi(source_json, "averaging", &
         default_averaging_str, averaging_str)
    averaging_type = averaging_type_from_str(averaging_str)

    call setup_table_source_control(source_json, srcstr, interpolation_type, &
         averaging_type, source_indices, source_controls, logfile, err)

    if (err == 0) then

       call setup_deliverability_source_controls(source_json, srcstr, &
            start_time, source_data, source_section, source_range_start, &
            fluid_data, fluid_section, fluid_range_start, &
            interpolation_type, averaging_type, source_indices, eos, &
            source_controls, logfile, err)

       if (err == 0) then

          call setup_recharge_source_controls(source_json, srcstr, &
               source_data, source_section, source_range_start, &
               fluid_data, fluid_section, fluid_range_start, &
               interpolation_type, averaging_type, source_indices, eos, &
               source_controls, logfile, err)

          if (err == 0) then

             call setup_limiter_source_controls(source_json, srcstr, thermo, &
                  source_indices, source_controls, logfile)

             call setup_direction_source_control(source_json, srcstr, thermo, &
                  source_indices, source_controls, logfile)

          end if

       end if

    end if

  end subroutine setup_inline_source_controls

!------------------------------------------------------------------------

  subroutine setup_table_source_control(source_json, srcstr, &
       interpolation_type, averaging_type, source_indices, &
       source_controls, logfile, err)
    !! Set up rate or enthalpy table source controls.

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    PetscInt, intent(in) :: interpolation_type, averaging_type
    PetscInt, intent(in) :: source_indices(:)
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(source_control_rate_table_type), pointer :: rate_control
    type(source_control_enthalpy_table_type), pointer :: enthalpy_control
    PetscInt :: variable_type
    PetscReal, allocatable :: rate_data_array(:,:), enthalpy_data_array(:,:)
    type(fson_value), pointer :: table

    ! Rate table:
    if (fson_has_mpi(source_json, "rate")) then
       variable_type = fson_type_mpi(source_json, "rate")
       if (variable_type == TYPE_ARRAY) then
          call fson_get_mpi(source_json, "rate", val = rate_data_array)
       else if (variable_type == TYPE_OBJECT) then
          call fson_get_mpi(source_json, "rate", table)
          if (fson_has_mpi(table, "time")) then
             call fson_get_mpi(table, "time", val = rate_data_array)
          end if
       end if
    end if

    if (allocated(rate_data_array)) then
       allocate(rate_control)
       call rate_control%init(rate_data_array, interpolation_type, &
            averaging_type, source_indices, err)
       if (err == 0) call source_controls%append(rate_control)
    end if

    if (err == 0) then

       ! Enthalpy table:
       if (fson_has_mpi(source_json, "enthalpy")) then
          variable_type = fson_type_mpi(source_json, "enthalpy")
          if (variable_type == TYPE_ARRAY) then
             call fson_get_mpi(source_json, "enthalpy", val = enthalpy_data_array)
          else if (variable_type == TYPE_OBJECT) then
             call fson_get_mpi(source_json, "enthalpy", table)
             if (fson_has_mpi(table, "time")) then
                call fson_get_mpi(table, "time", val = enthalpy_data_array)
             end if
          end if
       end if

       if (allocated(enthalpy_data_array)) then
          allocate(enthalpy_control)
          call enthalpy_control%init(enthalpy_data_array, interpolation_type, &
               averaging_type, source_indices, err)
          if (err == 0) call source_controls%append(enthalpy_control)
       end if

       if (err > 0) then
          call logfile%write(LOG_LEVEL_ERR, "input", "unsorted_array", &
               real_array_key = trim(srcstr) // "enthalpy", &
               real_array_value = enthalpy_data_array(:, 1))
       end if

    else
       call logfile%write(LOG_LEVEL_ERR, "input", "unsorted_array", &
            real_array_key = trim(srcstr) // "rate", &
            real_array_value = rate_data_array(:, 1))
    end if

    if (allocated(rate_data_array)) deallocate(rate_data_array)
    if (allocated(enthalpy_data_array)) deallocate(enthalpy_data_array)

  end subroutine setup_table_source_control

!------------------------------------------------------------------------

  subroutine get_reference_pressure(json, srcstr, source_type, &
       reference_pressure_array, calculate_reference_pressure, &
       pressure_table_coordinate, logfile)
    !! Get reference pressure for deliverability or recharge source control.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: json
    character(len = *), intent(in) :: srcstr, source_type
    PetscReal, allocatable, intent(out) :: reference_pressure_array(:,:)
    PetscBool, intent(out) :: calculate_reference_pressure
    PetscInt, intent(out) :: pressure_table_coordinate
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscReal :: reference_pressure
    PetscInt :: pressure_type
    PetscInt, parameter :: max_pressure_str_length = 8
    character(max_pressure_str_length) :: pressure_str
    type(fson_value), pointer :: pressure_json
    PetscReal, parameter :: default_time = 0._dp

    calculate_reference_pressure = PETSC_FALSE
    reference_pressure = default_deliverability_reference_pressure
    reference_pressure_array = reshape([default_time, &
         reference_pressure], [1,2])
    pressure_table_coordinate = default_source_pressure_table_coordinate

    if (fson_has_mpi(json, "pressure")) then
       pressure_type = fson_type_mpi(json, "pressure")
       select case (pressure_type)
       case (TYPE_REAL, TYPE_INTEGER)
          call fson_get_mpi(json, "pressure", &
               val = reference_pressure)
          reference_pressure_array(1,2) = reference_pressure
          pressure_table_coordinate = SRC_PRESSURE_TABLE_COORD_TIME
       case (TYPE_ARRAY)
          call fson_get_mpi(json, "pressure", &
               val = reference_pressure_array)
          pressure_table_coordinate = SRC_PRESSURE_TABLE_COORD_TIME
       case (TYPE_OBJECT)
          call fson_get_mpi(json, "pressure", pressure_json)
          if (fson_has_mpi(pressure_json, "time")) then
             call fson_get_mpi(pressure_json, "time", &
                  val = reference_pressure_array)
             pressure_table_coordinate = SRC_PRESSURE_TABLE_COORD_TIME
          else if (fson_has_mpi(pressure_json, "enthalpy")) then
             call fson_get_mpi(pressure_json, "enthalpy", &
                  val = reference_pressure_array)
             pressure_table_coordinate = SRC_PRESSURE_TABLE_COORD_ENTHALPY
          end if
       case (TYPE_STRING)
          call fson_get_mpi(json, "pressure", &
               val = pressure_str)
          if (trim(str_to_lower(pressure_str)) == 'initial') then
             calculate_reference_pressure = PETSC_TRUE
             pressure_table_coordinate = SRC_PRESSURE_TABLE_COORD_TIME
          else
             if (present(logfile) .and. logfile%active) then
                call logfile%write(LOG_LEVEL_WARN, 'input', 'unrecognised', &
                     str_key = trim(srcstr) // source_type // ".pressure", &
                     str_value = pressure_str)
             end if
          end if
       end select
    else
       if (present(logfile) .and. logfile%active) then
          call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
               real_keys = [trim(srcstr) // source_type // ".pressure"], &
               real_values = [reference_pressure])
       end if
    end if

  end subroutine get_reference_pressure

!------------------------------------------------------------------------
  
  subroutine get_deliverability_productivity(json, source_json, srcstr, &
       productivity_array, calculate_PI_from_rate, logfile)
    !! Gets productivity index for deliverability source control.

    type(fson_value), pointer, intent(in) :: json, source_json
    character(len = *), intent(in) :: srcstr
    PetscReal, allocatable, intent(out) :: productivity_array(:,:)
    PetscBool, intent(out) :: calculate_PI_from_rate
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: PI_type
    PetscReal :: productivity
    type(fson_value), pointer :: PI_json
    PetscReal, parameter :: default_time = 0._dp

    calculate_PI_from_rate = PETSC_FALSE
    productivity = default_deliverability_productivity
    productivity_array = reshape([default_time, productivity], [1,2])

    if (fson_has_mpi(json, "productivity")) then

       PI_type = fson_type_mpi(json, "productivity")

       select case(PI_type)
       case (TYPE_REAL, TYPE_INTEGER)
          call fson_get_mpi(json, "productivity", val = productivity)
          productivity_array = reshape([default_time, &
               productivity], [1,2])
       case (TYPE_ARRAY)
          call fson_get_mpi(json, "productivity", &
               val = productivity_array)
       case (TYPE_OBJECT)
          call fson_get_mpi(json, "productivity", PI_json)
          if (fson_has_mpi(PI_json, "time")) then
             call fson_get_mpi(PI_json, "time", &
                  val = productivity_array)
          end if
       case default
          if (present(logfile) .and. logfile%active) then
             call logfile%write(LOG_LEVEL_WARN, 'input', 'unrecognised', &
                  str_key = trim(srcstr) // "deliverability.productivity", &
                  str_value = "...")
          end if
       end select

    else if (fson_has_mpi(source_json, "rate")) then

       calculate_PI_from_rate = PETSC_TRUE

    else
       if (present(logfile) .and. logfile%active) then
          call logfile%write(LOG_LEVEL_INFO, 'input', 'default', real_keys = &
               [trim(srcstr) // "deliverability.productivity"], &
               real_values = [productivity])
       end if
    end if

  end subroutine get_deliverability_productivity

!------------------------------------------------------------------------

  subroutine setup_deliverability_source_controls(source_json, srcstr, &
       start_time, source_data, source_section, source_range_start, &
       fluid_data, fluid_section, fluid_range_start, &
       interpolation_type, averaging_type, source_indices, eos, &
       source_controls, logfile, err)
    !! Set up deliverability source controls. Deliverability controls
    !! can control only one source, so if multiple cells are
    !! specified, multiple corresponding deliverability controls are
    !! created.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    PetscReal, intent(in) :: start_time
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    PetscInt, intent(in) :: interpolation_type, averaging_type
    PetscInt, intent(in) :: source_indices(:)
    class(eos_type), intent(in) :: eos
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(fson_value), pointer :: deliv_json
    PetscReal, allocatable :: reference_pressure_array(:,:)
    type(source_control_deliverability_type), pointer :: deliv
    PetscBool :: calculate_reference_pressure
    PetscBool :: calculate_PI_from_rate
    PetscReal :: initial_rate, threshold
    PetscReal, allocatable :: productivity_array(:,:)
    PetscInt :: pressure_table_coordinate, i, s
    PetscReal, parameter :: default_rate = 0._dp
    PetscReal, parameter :: default_threshold = -1._dp

    if (fson_has_mpi(source_json, "deliverability")) then

       call fson_get_mpi(source_json, "deliverability", deliv_json)

       call fson_get_mpi(deliv_json, "threshold", default_threshold, threshold)

       if (threshold <= 0._dp) then
          call fson_get_mpi(source_json, "rate", default_rate, initial_rate)
       end if

       call get_reference_pressure(deliv_json, srcstr, "deliverability", &
            reference_pressure_array, calculate_reference_pressure, &
            pressure_table_coordinate, logfile)

       call get_deliverability_productivity(deliv_json, source_json, &
            srcstr, productivity_array, calculate_PI_from_rate, logfile)

       do i = 1, size(source_indices)

          s = source_indices(i)
          allocate(deliv)
          call deliv%init(productivity_array, interpolation_type, &
               averaging_type, reference_pressure_array, &
               pressure_table_coordinate, threshold, s, err)

          select case (err)
          case (0)
             if (calculate_reference_pressure) then
                call deliv%set_reference_pressure_initial(source_data, &
                     source_section, source_range_start, fluid_data, &
                     fluid_section, fluid_range_start, eos)
             end if
             if (calculate_PI_from_rate) then
                call deliv%calculate_PI_from_rate(start_time, initial_rate, &
                     source_data, source_section, source_range_start, &
                     fluid_data, fluid_section, fluid_range_start, eos, &
                     deliv%productivity%val(1, 1))
             end if
             if (deliv%threshold > 0._dp) then
                deliv%threshold_productivity = &
                     deliv%productivity%interpolate(start_time, 1)
             end if
             call source_controls%append(deliv)
          case (1)
             call deliv%destroy()
             deallocate(deliv)
             call logfile%write(LOG_LEVEL_ERR, "input", "unsorted_array", &
                  real_array_key = trim(srcstr) // "deliverability.productivity", &
                  real_array_value = productivity_array(:, 1))
          case (2)
             call deliv%destroy()
             deallocate(deliv)
             call logfile%write(LOG_LEVEL_ERR, "input", "unsorted_array", &
                  real_array_key = trim(srcstr) // "deliverability.pressure", &
                  real_array_value = reference_pressure_array(:, 1))
          end select

       end do

    end if

  end subroutine setup_deliverability_source_controls

!------------------------------------------------------------------------

  subroutine get_recharge_coefficient(json, source_json, srcstr, &
       recharge_array, logfile)
    !! Gets recharge coefficient for recharge source control.

    type(fson_value), pointer, intent(in) :: json, source_json
    character(len = *), intent(in) :: srcstr
    PetscReal, allocatable :: recharge_array(:,:)
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: coef_type
    PetscReal :: recharge_coefficient
    type(fson_value), pointer :: coef_json
    PetscReal, parameter :: default_time = 0._dp

    recharge_coefficient = default_recharge_coefficient
    recharge_array = reshape([default_time, &
         recharge_coefficient], [1,2])

    if (fson_has_mpi(json, "coefficient")) then

       coef_type = fson_type_mpi(json, "coefficient")

       select case(coef_type)
       case (TYPE_REAL, TYPE_INTEGER)
          call fson_get_mpi(json, "coefficient", &
               val = recharge_coefficient)
          recharge_array(1,2) = recharge_coefficient
       case (TYPE_ARRAY)
          call fson_get_mpi(json, "coefficient", &
               val = recharge_array)
       case (TYPE_OBJECT)
          call fson_get_mpi(json, "coefficient", coef_json)
          if (fson_has_mpi(coef_json, "time")) then
             call fson_get_mpi(coef_json, "time", &
                  val = recharge_array)
          end if
       case default
          if (present(logfile) .and. logfile%active) then
             call logfile%write(LOG_LEVEL_WARN, 'input', 'unrecognised', &
                  str_key = trim(srcstr) // "recharge.coefficient", &
                  str_value = "...")
          end if
       end select

    else
       if (present(logfile) .and. logfile%active) then
          call logfile%write(LOG_LEVEL_INFO, 'input', 'default', real_keys = &
               [trim(srcstr) // "recharge.coefficient"], &
               real_values = [recharge_coefficient])
       end if
    end if

  end subroutine get_recharge_coefficient

!------------------------------------------------------------------------

  subroutine setup_recharge_source_controls(source_json, srcstr, &
       source_data, source_section, source_range_start, &
       fluid_data, fluid_section, fluid_range_start, &
       interpolation_type, averaging_type, source_indices, eos, &
       source_controls, logfile, err)
    !! Set up recharge source controls. Recharge controls
    !! can control only one source, so if multiple cells are
    !! specified, multiple corresponding recharge controls are
    !! created.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    PetscInt, intent(in) :: interpolation_type, averaging_type
    PetscInt, intent(in) :: source_indices(:)
    class(eos_type), intent(in) :: eos
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(fson_value), pointer :: recharge_json
    PetscReal, allocatable :: reference_pressure_array(:,:)
    type(source_control_recharge_type), pointer :: recharge
    PetscBool :: calculate_reference_pressure
    PetscReal, allocatable :: recharge_array(:,:)
    PetscInt :: pressure_table_coordinate, i, s

    if (fson_has_mpi(source_json, "recharge")) then

       call fson_get_mpi(source_json, "recharge", recharge_json)

       call get_reference_pressure(recharge_json, srcstr, "recharge", &
            reference_pressure_array, calculate_reference_pressure, &
            pressure_table_coordinate, logfile)

       if (pressure_table_coordinate == SRC_PRESSURE_TABLE_COORD_TIME) then

          call get_recharge_coefficient(recharge_json, source_json, &
               srcstr, recharge_array, logfile)

          do i = 1, size(source_indices)

             s = source_indices(i)
             allocate(recharge)
             call recharge%init(recharge_array, interpolation_type, &
                  averaging_type, reference_pressure_array, s, err)

             select case (err)
             case (0)
                if (calculate_reference_pressure) then
                   call recharge%set_reference_pressure_initial(source_data, &
                        source_section, source_range_start, fluid_data, &
                        fluid_section, fluid_range_start, eos)
                end if
                call source_controls%append(recharge)
             case (1)
                call recharge%destroy()
                deallocate(recharge)
                call logfile%write(LOG_LEVEL_ERR, "input", "unsorted_array", &
                     real_array_key = trim(srcstr) // "recharge.coefficient", &
                     real_array_value = recharge_array(:, 1))
             case (2)
                call recharge%destroy()
                deallocate(recharge)
                call logfile%write(LOG_LEVEL_ERR, "input", "unsorted_array", &
                     real_array_key = trim(srcstr) // "recharge.pressure", &
                     real_array_value = reference_pressure_array(:, 1))
             end select

          end do

       else
          if (present(logfile) .and. logfile%active) then
             call logfile%write(LOG_LEVEL_WARN, 'input', 'not_supported', &
                  str_key = trim(srcstr) // "recharge.pressure", &
                  str_value = "...")
          end if
       end if

    end if

  end subroutine setup_recharge_source_controls

!------------------------------------------------------------------------

  subroutine setup_limiter_source_controls(source_json, srcstr, &
       thermo, source_indices, source_controls, logfile)
    !! Set up limiter source control for each cell source. If
    !! separated water or steam is to be limited, first set up
    !! corresponding separator for each cell source.

    use utils_module, only: str_to_lower, str_array_index

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    class(thermodynamics_type), intent(in) :: thermo
    PetscInt, intent(in) :: source_indices(:)
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    type(fson_value), pointer :: limiter_json
    character(max_phase_name_length) :: limiter_type_str
    PetscInt :: limiter_type
    PetscReal :: limit, separator_pressure
    PetscInt :: i, s
    type(source_control_limiter_type), pointer :: limiter
    type(source_control_separator_type), pointer :: separator

    if (fson_has_mpi(source_json, "limiter")) then

       call fson_get_mpi(source_json, "limiter", limiter_json)

       call fson_get_mpi(limiter_json, "type", &
            default_source_control_limiter_type_str, &
            limiter_type_str, logfile, srcstr)
       limiter_type_str = str_to_lower(limiter_type_str)

       select case (limiter_type_str)
       case ("total")
          limiter_type = SRC_CONTROL_LIMITER_TYPE_TOTAL
       case ("water")
          limiter_type = SRC_CONTROL_LIMITER_TYPE_WATER
       case ("steam")
          limiter_type = SRC_CONTROL_LIMITER_TYPE_STEAM
       case default
          limiter_type = SRC_CONTROL_LIMITER_TYPE_TOTAL
       end select

       if (limiter_type /= SRC_CONTROL_LIMITER_TYPE_TOTAL) then
          call fson_get_mpi(limiter_json, "separator_pressure", &
               default_source_control_separator_pressure, &
               separator_pressure, logfile, srcstr)
       end if

       call fson_get_mpi(limiter_json, "limit", &
            default_source_control_limiter_limit, limit, logfile, srcstr)

       do i = 1, size(source_indices)

          s = source_indices(i)
          allocate(limiter)
          if (limiter_type == SRC_CONTROL_LIMITER_TYPE_TOTAL) then
             call limiter%init(limiter_type, s, limit, [s])
          else
             allocate(separator)
             call separator%init(s, thermo, separator_pressure)
             call source_controls%append(separator)
             call limiter%init(limiter_type, separator, limit, [s])
          end if
          call source_controls%append(limiter)

       end do

    end if

  end subroutine setup_limiter_source_controls

!------------------------------------------------------------------------

  subroutine setup_direction_source_control(source_json, srcstr, &
       thermo, source_indices, source_controls, logfile)
    !! Set up direction source control. This can control multiple
    !! sources, so only one is created.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    class(thermodynamics_type), intent(in) :: thermo
    PetscInt, intent(in) :: source_indices(:)
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: direction
    type(source_control_direction_type), pointer :: direction_control
    PetscInt, parameter :: max_direction_str_length = 16
    character(max_direction_str_length) :: direction_str

    if (fson_has_mpi(source_json, "direction")) then

       call fson_get_mpi(source_json, "direction", val = direction_str)
       select case (trim(str_to_lower(direction_str)))
       case ("production", "out")
          direction = SRC_DIRECTION_PRODUCTION
       case ("injection", "in")
          direction = SRC_DIRECTION_INJECTION
       case ("both")
          direction = SRC_DIRECTION_BOTH
       case default
          if (present(logfile) .and. logfile%active) then
             call logfile%write(LOG_LEVEL_WARN, 'input', 'unrecognised', &
                  str_key = trim(srcstr) // "direction", &
                  str_value = direction_str)
          end if
          direction = default_source_direction
       end select

       if (direction /= SRC_DIRECTION_BOTH) then
         allocate(direction_control)
         call direction_control%init(direction, source_indices)
         call source_controls%append(direction_control)
       end if

    end if

  end subroutine setup_direction_source_control

!------------------------------------------------------------------------

end module source_setup_module
