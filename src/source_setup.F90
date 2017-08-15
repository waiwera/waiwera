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

  subroutine setup_sources(json, dm, eos, thermo, start_time, fluid_vector, &
       fluid_range_start, sources, source_controls, logfile)
    !! Sets up lists of sinks / sources and source controls.

    use cell_order_module, only: cell_order_label_name
    use dm_utils_module, only: dm_order_local_index, global_vec_section

    type(fson_value), pointer, intent(in) :: json !! JSON file object
    DM, intent(in) :: dm !! Mesh DM
    class(eos_type), intent(in) :: eos !! Equation of state
    class(thermodynamics_type), intent(in) :: thermo !! Thermodynamics formulation
    PetscReal, intent(in) :: start_time
    Vec, intent(in) :: fluid_vector !! Fluid vector
    PetscInt, intent(in) :: fluid_range_start
    type(list_type), intent(in out) :: sources !! List of sources
    type(list_type), intent(in out) :: source_controls !! List of source controls
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    PetscInt :: icell, c, isrc, cell_order
    PetscInt :: injection_component, production_component
    type(fson_value), pointer :: sources_json, source_json
    PetscInt :: num_sources, num_cells
    type(source_type), pointer :: source
    PetscReal :: initial_rate, initial_enthalpy
    character(max_source_name_length) :: name
    character(max_source_name_length) :: default_name = ""
    DMLabel :: ghost_label
    character(len=64) :: srcstr
    character(len=12) :: istr
    PetscBool :: can_inject
    PetscInt, allocatable :: cells(:)
    type(list_type) :: cell_sources
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscSection :: fluid_section
    PetscErrorCode :: ierr

    call sources%init(owner = PETSC_TRUE)
    call source_controls%init(owner = PETSC_TRUE)

    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayReadF90(fluid_vector, fluid_data, ierr); CHKERRQ(ierr)

    if (fson_has_mpi(json, "source")) then

       call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)

       call fson_get_mpi(json, "source", sources_json)
       num_sources = fson_value_count_mpi(sources_json, ".")

       do isrc = 1, num_sources

          write(istr, '(i0)') isrc - 1
          srcstr = 'source[' // trim(istr) // '].'
          source => null()

          source_json => fson_value_get_mpi(sources_json, isrc)

          call fson_get_mpi(source_json, "name", default_name, name)
          call get_components(source_json, eos, &
               injection_component, production_component, logfile)
          call get_initial_rate(source_json, initial_rate, can_inject)
          call get_initial_enthalpy(source_json, eos, can_inject, &
               injection_component, initial_enthalpy)

          call get_cells(source_json, cells, num_cells)
          call cell_sources%init()

          do icell = 1, num_cells

             cell_order = cells(icell)
             c = dm_order_local_index(dm, cell_order, &
                  cell_order_label_name, ghost_label)
             if (c >= 0) then
                allocate(source)
                call source%init(cell_order, c, eos, &
                     initial_rate, initial_enthalpy, &
                     injection_component, production_component)
                call sources%append(source, name)
                call cell_sources%append(source, name)
             end if

          end do

          call setup_inline_source_controls(source_json, eos, thermo, &
               start_time, fluid_data, fluid_section, fluid_range_start, srcstr, &
               num_cells, cell_sources, source_controls, logfile)
          call cell_sources%destroy()

       end do

    else
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_INFO, "input", "no_sources")
       end if
    end if

    call VecRestoreArrayReadF90(fluid_vector, fluid_data, ierr)
    CHKERRQ(ierr)

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

  subroutine get_cells(source_json, cells, num_cells)
    !! Gets array of cell indices for the source. These are global
    !! indices (i.e. natural orders).  The "cell" key can be used to
    !! specify a single cell. The "cells" key can be used to specify
    !! either a single cell or an array of cells. If both are present,
    !! the "cells" key is used.

    type(fson_value), pointer, intent(in) :: source_json
    PetscInt, allocatable, intent(out) :: cells(:)
    PetscInt, intent(out) :: num_cells
    ! Locals:
    PetscInt :: cell_type, cell

    if (fson_has_mpi(source_json, "cell")) then
       call fson_get_mpi(source_json, "cell", val = cell)
       cells = [cell]
    end if

    if (fson_has_mpi(source_json, "cells")) then

       cell_type = fson_type_mpi(source_json, "cells")

       if (cell_type == TYPE_INTEGER) then
          call fson_get_mpi(source_json, "cells", val = cell)
          cells = [cell]
       else if (cell_type == TYPE_ARRAY) then
          call fson_get_mpi(source_json, "cells", val = cells)
       end if

    end if

    if (allocated(cells)) then
       num_cells = size(cells)
    else
       num_cells = 0
    end if

  end subroutine get_cells

!------------------------------------------------------------------------

  subroutine setup_inline_source_controls(source_json, eos, thermo, &
       start_time, fluid_data, fluid_section, fluid_range_start, &
       srcstr, num_cells, cell_sources, source_controls, logfile)
    !! Sets up any 'inline' source controls for the source,
    !! i.e. controls defined implicitly in the specification of the source.

    use interpolation_module, only: interpolation_type_from_str, &
         averaging_type_from_str, max_interpolation_str_length, &
         max_averaging_str_length, default_interpolation_str, &
         default_averaging_str

    type(fson_value), pointer, intent(in) :: source_json
    class(eos_type), intent(in) :: eos
    class(thermodynamics_type), intent(in) :: thermo
    PetscReal, intent(in) :: start_time
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    character(len = *), intent(in) :: srcstr
    PetscInt, intent(in) :: num_cells
    type(list_type), intent(in out) :: cell_sources
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
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

    call setup_table_source_control(source_json, interpolation_type, &
         averaging_type, cell_sources, source_controls)

    call setup_deliverability_source_control(source_json, srcstr, &
         start_time, fluid_data, fluid_section, fluid_range_start, &
         interpolation_type, averaging_type, num_cells, cell_sources, &
         source_controls, logfile)

    call setup_recharge_source_control(source_json, srcstr, &
         fluid_data, fluid_section, fluid_range_start, &
         interpolation_type, averaging_type, num_cells, cell_sources, &
         source_controls, logfile)

    call setup_limiter_source_control(source_json, srcstr, thermo, &
         cell_sources, source_controls, logfile)

    call setup_direction_source_control(source_json, srcstr, thermo, &
         cell_sources, source_controls, logfile)

  end subroutine setup_inline_source_controls

!------------------------------------------------------------------------

  subroutine setup_table_source_control(source_json, interpolation_type, &
       averaging_type, cell_sources, source_controls)
    !! Set up rate or enthalpy table source controls.

    type(fson_value), pointer, intent(in) :: source_json
    PetscInt, intent(in) :: interpolation_type, averaging_type
    type(list_type), intent(in out) :: cell_sources
    type(list_type), intent(in out) :: source_controls
    ! Locals:
    type(source_control_rate_table_type), pointer :: rate_control
    type(source_control_enthalpy_table_type), pointer :: enthalpy_control
    PetscInt :: variable_type
    PetscReal, allocatable :: data_array(:,:)
    type(fson_value), pointer :: table

    ! Rate table:
    if (fson_has_mpi(source_json, "rate")) then
       variable_type = fson_type_mpi(source_json, "rate")
       if (variable_type == TYPE_ARRAY) then
          call fson_get_mpi(source_json, "rate", val = data_array)
       else if (variable_type == TYPE_OBJECT) then
          call fson_get_mpi(source_json, "rate", table)
          if (fson_has_mpi(table, "time")) then
             call fson_get_mpi(table, "time", val = data_array)
          end if
       end if
    end if

    if ((cell_sources%count > 0) .and. (allocated(data_array))) then
       allocate(rate_control)
       call rate_control%init(data_array, interpolation_type, &
            averaging_type, cell_sources)
       call source_controls%append(rate_control)
       deallocate(data_array)
    end if

    ! Enthalpy table:
    if (fson_has_mpi(source_json, "enthalpy")) then
       variable_type = fson_type_mpi(source_json, "enthalpy")
       if (variable_type == TYPE_ARRAY) then
          call fson_get_mpi(source_json, "enthalpy", val = data_array)
       else if (variable_type == TYPE_OBJECT) then
          call fson_get_mpi(source_json, "enthalpy", table)
          if (fson_has_mpi(table, "time")) then
             call fson_get_mpi(table, "time", val = data_array)
          end if
       end if
    end if

    if ((cell_sources%count > 0) .and. (allocated(data_array))) then
       allocate(enthalpy_control)
       call enthalpy_control%init(data_array, interpolation_type, &
            averaging_type, cell_sources)
       call source_controls%append(enthalpy_control)
       deallocate(data_array)
    end if

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
       num_cells, cell_sources, productivity_array, &
       calculate_PI_from_rate, logfile)
    !! Gets productivity index for deliverability source control.

    type(fson_value), pointer, intent(in) :: json, source_json
    character(len = *), intent(in) :: srcstr
    PetscInt, intent(in) :: num_cells
    type(list_type), intent(in) :: cell_sources
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

    else if ((fson_has_mpi(source_json, "rate") .and. &
         (num_cells == 1) .and. (cell_sources%count == 1))) then

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

  subroutine setup_deliverability_source_control(source_json, srcstr, &
       start_time, fluid_data, fluid_section, fluid_range_start, &
       interpolation_type, averaging_type, num_cells, &
       cell_sources, source_controls, logfile)
    !! Set up deliverability source control.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    PetscReal, intent(in) :: start_time
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    PetscInt, intent(in) :: interpolation_type, averaging_type
    PetscInt, intent(in) :: num_cells
    type(list_type), intent(in out) :: cell_sources
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    type(fson_value), pointer :: deliv_json
    PetscReal, allocatable :: reference_pressure_array(:,:)
    type(source_control_deliverability_type), pointer :: deliv
    PetscBool :: calculate_reference_pressure
    PetscBool :: calculate_PI_from_rate
    PetscReal :: initial_rate
    PetscReal, allocatable :: productivity_array(:,:)
    PetscInt :: pressure_table_coordinate
    PetscReal, parameter :: default_rate = 0._dp

    if (fson_has_mpi(source_json, "deliverability")) then

       call fson_get_mpi(source_json, "rate", default_rate, initial_rate)

       call fson_get_mpi(source_json, "deliverability", deliv_json)

       call get_reference_pressure(deliv_json, srcstr, "deliverability", &
            reference_pressure_array, calculate_reference_pressure, &
            pressure_table_coordinate, logfile)

       call get_deliverability_productivity(deliv_json, source_json, &
            srcstr, num_cells, cell_sources, productivity_array, &
            calculate_PI_from_rate, logfile)

       if (cell_sources%count > 0) then

          allocate(deliv)
          call deliv%init(productivity_array, interpolation_type, &
               averaging_type, reference_pressure_array, &
               pressure_table_coordinate, cell_sources)

          if (calculate_reference_pressure) then
             call deliv%set_reference_pressure_initial(fluid_data, &
                  fluid_section, fluid_range_start)
          end if

          if (calculate_PI_from_rate) then
             call deliv%calculate_PI_from_rate(start_time, initial_rate, &
                  fluid_data, fluid_section, fluid_range_start)
          end if

          call source_controls%append(deliv)

       end if

    end if

  end subroutine setup_deliverability_source_control

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

  subroutine setup_recharge_source_control(source_json, srcstr, &
       fluid_data, fluid_section, fluid_range_start, &
       interpolation_type, averaging_type, num_cells, &
       cell_sources, source_controls, logfile)
    !! Set up recharge source control.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    PetscInt, intent(in) :: interpolation_type, averaging_type
    PetscInt, intent(in) :: num_cells
    type(list_type), intent(in out) :: cell_sources
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    type(fson_value), pointer :: recharge_json
    PetscReal, allocatable :: reference_pressure_array(:,:)
    type(source_control_recharge_type), pointer :: recharge
    PetscBool :: calculate_reference_pressure
    PetscReal, allocatable :: recharge_array(:,:)
    PetscInt :: pressure_table_coordinate

    if (fson_has_mpi(source_json, "recharge")) then

       call fson_get_mpi(source_json, "recharge", recharge_json)

       call get_reference_pressure(recharge_json, srcstr, "recharge", &
            reference_pressure_array, calculate_reference_pressure, &
            pressure_table_coordinate, logfile)

       if (pressure_table_coordinate == SRC_PRESSURE_TABLE_COORD_TIME) then

          call get_recharge_coefficient(recharge_json, source_json, &
               srcstr, recharge_array, logfile)

          if (cell_sources%count > 0) then

             allocate(recharge)
             call recharge%init(recharge_array, interpolation_type, &
                  averaging_type, reference_pressure_array, cell_sources)

             if (calculate_reference_pressure) then
                call recharge%set_reference_pressure_initial(fluid_data, &
                     fluid_section, fluid_range_start)
             end if

             call source_controls%append(recharge)

          end if

       else
          if (present(logfile) .and. logfile%active) then
             call logfile%write(LOG_LEVEL_WARN, 'input', 'not_supported', &
                  str_key = trim(srcstr) // "recharge.pressure", &
                  str_value = "...")
          end if
       end if

    end if

  end subroutine setup_recharge_source_control

!------------------------------------------------------------------------

  subroutine setup_limiter_source_control(source_json, srcstr, &
       thermo, cell_sources, source_controls, logfile)
    !! Set up limiter source control for each cell source. If
    !! separated water or steam is to be limited, first set up
    !! corresponding separator for each cell source.

    use utils_module, only: str_to_lower, str_array_index

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    class(thermodynamics_type), intent(in) :: thermo
    type(list_type), intent(in out) :: cell_sources
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    type(fson_value), pointer :: limiter_json
    character(max_phase_name_length) :: limiter_type_str
    PetscInt :: limiter_type
    PetscReal :: limit, separator_pressure

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

       call cell_sources%traverse(setup_limiter_iterator)

    end if

  contains

    subroutine setup_limiter_iterator(node, stopped)
      !! Sets up limiter (and separator if needed) at a source list node.

      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      type(source_control_limiter_type), pointer :: limiter
      type(source_control_separator_type), pointer :: separator

      select type (source => node%data)
      type is (source_type)

         allocate(limiter)

         if (limiter_type == SRC_CONTROL_LIMITER_TYPE_TOTAL) then

            call limiter%init(limiter_type, source%rate, limit, &
                 cell_sources)

         else

            allocate(separator)
            call separator%init(source, thermo, separator_pressure)
            call source_controls%append(separator)

            if (limiter_type == SRC_CONTROL_LIMITER_TYPE_WATER) then
               call limiter%init(limiter_type, separator%water_flow_rate, &
                    limit, cell_sources)
            else
               call limiter%init(limiter_type, separator%steam_flow_rate, &
                    limit, cell_sources)
            end if

         end if

         call source_controls%append(limiter)

      end select

    end subroutine setup_limiter_iterator

  end subroutine setup_limiter_source_control

!------------------------------------------------------------------------

  subroutine setup_direction_source_control(source_json, srcstr, &
       thermo, cell_sources, source_controls, logfile)
    !! Set up direction source control for each cell source.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    class(thermodynamics_type), intent(in) :: thermo
    type(list_type), intent(in out) :: cell_sources
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: direction
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
          call cell_sources%traverse(setup_direction_iterator)
       end if

    end if

  contains

    subroutine setup_direction_iterator(node, stopped)
      !! Sets up direction control at a source list node.

      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      type(source_control_direction_type), pointer :: direction_control

      select type (source => node%data)
      type is (source_type)

         allocate(direction_control)
         call direction_control%init(direction, cell_sources)
         call source_controls%append(direction_control)

      end select

    end subroutine setup_direction_iterator

  end subroutine setup_direction_source_control

!------------------------------------------------------------------------

end module source_setup_module
