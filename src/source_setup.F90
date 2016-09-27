module source_setup_module
  !! Module for setting up sources and sinks.

  use kinds_module
  use fson
  use fson_value_m, only: TYPE_INTEGER, TYPE_REAL, TYPE_ARRAY
  use fson_mpi_module
  use list_module
  use logfile_module
  use source_module
  use source_control_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  public :: setup_sources

contains

!------------------------------------------------------------------------

  subroutine setup_sources(json, dm, np, isothermal, sources, &
       source_controls, logfile)
    !! Sets up lists of sinks / sources and source controls.

    use mesh_module, only: cell_order_label_name
    use dm_utils_module, only: dm_order_local_index

    type(fson_value), pointer, intent(in) :: json !! JSON file object
    DM, intent(in) :: dm !! Mesh DM
    PetscInt, intent(in) :: np !! Number of primary variables
    PetscBool, intent(in) :: isothermal !! Whether EOS is isothermal
    type(list_type), intent(in out) :: sources !! List of sources
    type(list_type), intent(in out) :: source_controls !! List of source controls
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: icell, c, isrc, cell_order
    PetscInt :: injection_component, production_component
    type(fson_value), pointer :: sources_json, source_json
    PetscInt :: num_sources, num_cells
    type(source_type), pointer :: source
    PetscReal :: rate, enthalpy
    character(max_source_name_length) :: name
    character(max_source_name_length) :: default_name = ""
    DMLabel :: ghost_label
    character(len=64) :: srcstr
    character(len=12) :: istr
    PetscBool :: can_inject
    PetscInt, allocatable :: cells(:)
    type(list_type) :: cell_sources

    call sources%init(delete_deallocates = PETSC_TRUE)
    call source_controls%init(delete_deallocates = PETSC_TRUE)

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
          call get_rate(source_json, rate, can_inject)
          call get_components(source_json, isothermal, np, srcstr, &
               injection_component, production_component, logfile)
          call get_enthalpy(source_json, isothermal, can_inject, &
               injection_component, np, srcstr, enthalpy, logfile)

          call get_cells(source_json, cells)
          num_cells = size(cells)
          call cell_sources%init()

          do icell = 1, num_cells

             cell_order = cells(icell)
             c = dm_order_local_index(dm, cell_order, &
                  cell_order_label_name, ghost_label)
             if (c >= 0) then
                allocate(source)
                call source%init(cell_order, c, np, rate, enthalpy, &
                     injection_component, production_component)
                call sources%append(source, name)
                call cell_sources%append(source, name)
             end if

          end do

          if (cell_sources%count > 0) then
             call setup_inline_source_controls(source_json, cell_sources, &
                  source_controls)
          end if
          call cell_sources%destroy()

       end do

    else
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_INFO, "input", "no_sources")
       end if
    end if

  end subroutine setup_sources

!------------------------------------------------------------------------

  subroutine get_rate(source_json, rate, can_inject)
    !! Gets flow rate and whether the source can inject. These can
    !! only be determined for constant-rate sources. For other source
    !! types, a default initial rate is assigned and it is assumed
    !! that the source may inject.

    type(fson_value), pointer, intent(in) :: source_json
    PetscReal, intent(out) :: rate
    PetscBool, intent(out) :: can_inject
    ! Locals:
    PetscInt :: rate_type

    if (fson_has_mpi(source_json, "rate")) then

       rate_type = fson_type_mpi(source_json, "rate")

       if (rate_type == TYPE_REAL) then
          call fson_get_mpi(source_json, "rate", val = rate)
          can_inject = (rate > 0._dp)
       else if (rate_type == TYPE_ARRAY) then
          rate = default_source_rate
          can_inject = PETSC_TRUE
       end if

    else
       rate = default_source_rate
       can_inject = PETSC_TRUE
    end if

  end subroutine get_rate

!------------------------------------------------------------------------

  subroutine get_components(source_json, isothermal, &
       np, srcstr, injection_component, production_component, &
       logfile)
    !! Gets injection and production components for the source.

    type(fson_value), pointer, intent(in) :: source_json
    PetscBool, intent(in) :: isothermal
    PetscInt, intent(in) :: np
    character(len=*) :: srcstr
    PetscInt, intent(out) :: injection_component
    PetscInt, intent(out) :: production_component
    type(logfile_type), intent(in out), optional :: logfile

    call fson_get_mpi(source_json, "component", default_source_component, &
         injection_component, logfile, trim(srcstr) // "component")

    if (fson_has_mpi(source_json, "production_component")) then

       call fson_get_mpi(source_json, "production_component", &
            val = production_component)

    else

       if ((.not. isothermal) .and. (injection_component == np)) then
          production_component = injection_component
       else
          production_component = default_source_production_component
       end if

    end if

  end subroutine get_components

!------------------------------------------------------------------------

  subroutine get_enthalpy(source_json, isothermal, can_inject, &
       injection_component, np, srcstr, enthalpy, logfile)
    !! Gets production enthalpy for the source, if needed.

    type(fson_value), pointer, intent(in) :: source_json
    PetscBool, intent(in) :: isothermal, can_inject
    PetscInt, intent(in) :: injection_component, np
    character(len=*) :: srcstr
    PetscReal, intent(out) :: enthalpy
    type(logfile_type), intent(in out), optional :: logfile

    if (can_inject .and. (injection_component < np)) then

       call fson_get_mpi(source_json, "enthalpy", &
            default_source_injection_enthalpy, &
            enthalpy, logfile, trim(srcstr) // "enthalpy")

    else
       enthalpy = 0._dp
    end if

  end subroutine get_enthalpy

!------------------------------------------------------------------------

  subroutine get_cells(source_json, cells)
    !! Gets array of cell indices for the source. These are global
    !! indices (i.e. natural orders).

    type(fson_value), pointer, intent(in) :: source_json
    PetscInt, allocatable, intent(out) :: cells(:)
    ! Locals:
    PetscInt :: cell_type, cell

    if (fson_has_mpi(source_json, "cell")) then

       cell_type = fson_type_mpi(source_json, "cell")

       if (cell_type == TYPE_INTEGER) then
          call fson_get_mpi(source_json, "cell", val = cell)
          cells = [cell]
       else if (cell_type == TYPE_ARRAY) then
          call fson_get_mpi(source_json, "cell", val = cells)
       end if

    end if

  end subroutine get_cells

!------------------------------------------------------------------------

  subroutine setup_inline_source_controls(source_json, cell_sources, &
       source_controls)
    !! Sets up any 'inline' source controls for the source,
    !! i.e. controls defined implicitly in the specification of the source.

    use interpolation_module, only: INTERP_LINEAR, INTERP_AVERAGING_INTEGRATE

    type(fson_value), pointer, intent(in) :: source_json
    type(list_type), intent(in out) :: cell_sources
    type(list_type), intent(in out) :: source_controls
    ! Locals:
    type(source_control_rate_table_type), pointer :: rate_control
    PetscInt :: rate_type, interpolation_type, averaging_type
    PetscReal, allocatable :: rate_array(:,:)
    PetscInt, parameter :: default_interpolation_type = INTERP_LINEAR
    PetscInt, parameter :: default_averaging_type = INTERP_AVERAGING_INTEGRATE

    ! Rate table:
    if (fson_has_mpi(source_json, "rate")) then
       rate_type = fson_type_mpi(source_json, "rate")
       if (rate_type == TYPE_ARRAY) then

          call fson_get_mpi(source_json, "rate", val = rate_array)
          call fson_get_mpi(source_json, "interpolation", &
               default_interpolation_type, interpolation_type)
          call fson_get_mpi(source_json, "averaging", &
               default_averaging_type, averaging_type)

          allocate(source_control_rate_table_type :: rate_control)
          call rate_control%init()
          call rate_control%table%init(rate_array, interpolation_type, &
               averaging_type)
          call rate_control%sources%add(cell_sources)
          call source_controls%append(rate_control)

       end if
    end if

    ! TODO: other inline controls (enthalpy table, limiter, etc.)

  end subroutine setup_inline_source_controls

!------------------------------------------------------------------------

end module source_setup_module
