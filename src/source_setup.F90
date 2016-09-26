module source_setup_module
  !! Module for setting up sources and sinks.

  use kinds_module
  use fson
  use fson_value_m, only: TYPE_REAL, TYPE_ARRAY
  use fson_mpi_module
  use list_module
  use logfile_module
  use mesh_module, only: cell_order_label_name
  use source_module
  use source_control_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  public :: setup_sources

contains

!------------------------------------------------------------------------

  subroutine setup_sources(json, dm, np, isothermal, sources_list, &
       source_controls_list, logfile)
    !! Sets up lists of sinks / sources and source controls.

    type(fson_value), pointer, intent(in) :: json !! JSON file object
    DM, intent(in) :: dm !! Mesh DM
    PetscInt, intent(in) :: np !! Number of primary variables
    PetscBool, intent(in) :: isothermal !! Whether EOS is isothermal
    type(list_type), intent(in out) :: sources_list !! List of sources
    type(list_type), intent(in out) :: source_controls_list !! List of source controls
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: c, isrc, injection_component, production_component
    type(fson_value), pointer :: sources_json, source_json
    PetscInt :: num_sources, cell, num_cells, ghost
    type(source_type), pointer :: source
    class(source_control_type), pointer :: control
    PetscInt, pointer :: cells(:)
    PetscReal :: rate, enthalpy
    character(max_source_name_length) :: name
    character(max_source_name_length) :: default_name = ""
    IS :: cell_IS
    DMLabel :: ghost_label
    character(len=64) :: srcstr
    character(len=12) :: istr
    PetscBool :: can_inject

    call sources_list%init(delete_deallocates = PETSC_TRUE)
    call source_controls_list%init(delete_deallocates = PETSC_TRUE)

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
          call fson_get_mpi(source_json, "cell", val = cell)

          call get_rate(source_json, rate, can_inject)

          call get_components(source_json, isothermal, np, srcstr, &
               injection_component, production_component, logfile)

          call get_enthalpy(source_json, isothermal, can_inject, &
               injection_component, np, srcstr, enthalpy, logfile)

          call DMGetStratumIS(dm, cell_order_label_name, &
               cell, cell_IS, ierr); CHKERRQ(ierr)
          if (cell_IS /= 0) then

             call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
             num_cells = size(cells)
             if ((num_cells > 1) .and. present(logfile)) then
                call logfile%write(LOG_LEVEL_WARN, 'source', &
                     'multiple cells found for source: ' // trim(name))
             end if

             c = cells(1)
             call DMLabelGetValue(ghost_label, c, ghost, ierr)
             CHKERRQ(ierr)
             if (ghost < 0) then
                allocate(source)
                call source%init(cell, c, np, rate, enthalpy, &
                     injection_component, production_component)
                call sources_list%append(source, name)
             end if
             call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)

          end if
          call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)

          if (associated(source)) then
             ! Set up any inline controls for this source:

          end if

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

end module source_setup_module
