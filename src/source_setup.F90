module source_setup_module
  !! Module for setting up sources and sinks.

  use kinds_module
  use fson
  use fson_value_m, only: TYPE_INTEGER, TYPE_REAL, TYPE_ARRAY, TYPE_STRING
  use fson_mpi_module
  use list_module
  use logfile_module
  use eos_module
  use thermodynamics_module, only: thermodynamics_type
  use source_module
  use source_control_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  public :: setup_sources

contains

!------------------------------------------------------------------------

  subroutine setup_sources(json, dm, eos, thermo, sources, &
       source_controls, logfile)
    !! Sets up lists of sinks / sources and source controls.

    use mesh_module, only: cell_order_label_name
    use dm_utils_module, only: dm_order_local_index

    type(fson_value), pointer, intent(in) :: json !! JSON file object
    DM, intent(in) :: dm !! Mesh DM
    class(eos_type), intent(in) :: eos
    class(thermodynamics_type), intent(in) :: thermo
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
    PetscErrorCode :: ierr

    call sources%init(owner = PETSC_TRUE)
    call source_controls%init(owner = PETSC_TRUE)

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
               srcstr, cell_sources, source_controls, logfile)
          call cell_sources%destroy()

       end do

    else
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_INFO, "input", "no_sources")
       end if
    end if

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

       if (rate_type == TYPE_REAL) then
          call fson_get_mpi(source_json, "rate", val = initial_rate)
          can_inject = (initial_rate > 0._dp)
       else if (rate_type == TYPE_ARRAY) then
          initial_rate = default_source_rate
          can_inject = PETSC_TRUE
       end if

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

            if (enthalpy_type == TYPE_REAL) then
               call fson_get_mpi(source_json, "enthalpy", &
                    val = enthalpy)
            else
               enthalpy = default_source_injection_enthalpy
            end if

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
       srcstr, cell_sources, source_controls, logfile)
    !! Sets up any 'inline' source controls for the source,
    !! i.e. controls defined implicitly in the specification of the source.

    type(fson_value), pointer, intent(in) :: source_json
    class(eos_type), intent(in) :: eos
    class(thermodynamics_type), intent(in) :: thermo
    character(len = *), intent(in) :: srcstr
    type(list_type), intent(in out) :: cell_sources
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile

    call setup_table_source_control(source_json, cell_sources, &
         source_controls)

    call setup_deliverability_source_control(source_json, srcstr, &
         cell_sources, source_controls, logfile)

    call setup_limiter_source_control(source_json, srcstr, thermo, &
         cell_sources, source_controls, logfile)

  end subroutine setup_inline_source_controls

!------------------------------------------------------------------------

  subroutine setup_table_source_control(source_json, cell_sources, &
       source_controls)
    !! Set up rate or enthalpy table source controls.

    use interpolation_module, only: interpolation_type_from_str, &
         averaging_type_from_str
    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    type(list_type), intent(in out) :: cell_sources
    type(list_type), intent(in out) :: source_controls
    ! Locals:
    type(source_control_rate_table_type), pointer :: rate_control
    type(source_control_enthalpy_table_type), pointer :: enthalpy_control
    PetscInt :: variable_type, interpolation_type, averaging_type
    PetscReal, allocatable :: data_array(:,:)
    PetscInt, parameter :: max_interpolation_str_length = 16
    character(max_interpolation_str_length), parameter :: &
         default_interpolation_str = "linear"
    character(max_interpolation_str_length) :: interpolation_str
    PetscInt, parameter :: max_averaging_str_length = 16
    character(max_averaging_str_length), parameter :: &
         default_averaging_str = "integrate"
    character(max_averaging_str_length) :: averaging_str

    ! Rate table:
    if (fson_has_mpi(source_json, "rate")) then
       variable_type = fson_type_mpi(source_json, "rate")
       if (variable_type == TYPE_ARRAY) then

          call fson_get_mpi(source_json, "rate", val = data_array)

          call fson_get_mpi(source_json, "interpolation", &
               default_interpolation_str, interpolation_str)
          interpolation_type = interpolation_type_from_str(interpolation_str)

          call fson_get_mpi(source_json, "averaging", &
               default_averaging_str, averaging_str)
          averaging_type = averaging_type_from_str(averaging_str)

          if (cell_sources%count > 0) then
             allocate(source_control_rate_table_type :: rate_control)
             call rate_control%init(data_array, interpolation_type, &
                  averaging_type, cell_sources)
             call source_controls%append(rate_control)
          end if

          deallocate(data_array)

       end if
    end if

    ! Enthalpy table:
    if (fson_has_mpi(source_json, "enthalpy")) then
       variable_type = fson_type_mpi(source_json, "enthalpy")
       if (variable_type == TYPE_ARRAY) then

          call fson_get_mpi(source_json, "enthalpy", val = data_array)

          call fson_get_mpi(source_json, "interpolation", &
               default_interpolation_str, interpolation_str)
          interpolation_type = interpolation_type_from_str(interpolation_str)

          call fson_get_mpi(source_json, "averaging", &
               default_averaging_str, averaging_str)
          averaging_type = averaging_type_from_str(averaging_str)

          if (cell_sources%count > 0) then
             allocate(source_control_enthalpy_table_type :: enthalpy_control)
             call enthalpy_control%init(data_array, interpolation_type, &
                  averaging_type, cell_sources)
             call source_controls%append(enthalpy_control)
          end if

          deallocate(data_array)

       end if
    end if

  end subroutine setup_table_source_control

!------------------------------------------------------------------------

  subroutine setup_deliverability_source_control(source_json, srcstr, &
       cell_sources, source_controls, logfile)
    !! Set up deliverability source control.

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    type(list_type), intent(in out) :: cell_sources
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    type(fson_value), pointer :: deliv_json
    PetscReal :: productivity_index, bottomhole_pressure
    type(source_control_deliverability_type), pointer :: deliv

    if (fson_has_mpi(source_json, "deliverability")) then

       call fson_get_mpi(source_json, "deliverability", deliv_json)

       call fson_get_mpi(deliv_json, "productivity_index", &
            default_deliverability_productivity_index, productivity_index, &
            logfile, srcstr)

       call fson_get_mpi(deliv_json, "bottomhole_pressure", &
            default_deliverability_bottomhole_pressure, bottomhole_pressure, &
            logfile, srcstr)

       if (cell_sources%count > 0) then
          allocate(source_control_deliverability_type :: deliv)
          call deliv%init(productivity_index, bottomhole_pressure, cell_sources)
          call source_controls%append(deliv)
       end if

    end if

  end subroutine setup_deliverability_source_control

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

         allocate(source_control_limiter_type :: limiter)

         if (limiter_type == SRC_CONTROL_LIMITER_TYPE_TOTAL) then

            call limiter%init(limiter_type, source%rate, limit, &
                 cell_sources)

         else

            allocate(source_control_separator_type :: separator)
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

end module source_setup_module
