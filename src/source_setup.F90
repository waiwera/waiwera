module source_setup_module
  !! Module for setting up sources and sinks.

  use source_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  public :: setup_sources

contains

!------------------------------------------------------------------------

  subroutine setup_sources(json, dm, np, isothermal, sources_list, &
       source_controls_list, logfile)
    !! Sets up lists of sinks / sources and source controls.

    use kinds_module
    use fson
    use fson_mpi_module
    use list_module
    use logfile_module
    use mesh_module, only: cell_order_label_name

    type(fson_value), pointer, intent(in) :: json !! JSON file object
    DM, intent(in) :: dm !! Mesh DM
    PetscInt, intent(in) :: np !! Number of primary variables
    PetscBool, intent(in) :: isothermal !! Whether EOS is isothermal
    type(list_type), intent(in out) :: sources_list !! List of sources
    type(list_type), intent(in out) :: source_controls_list !! List of source controls
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: c, isrc, i, injection_component, production_component
    type(fson_value), pointer :: sources, src
    PetscInt :: num_sources, cell, num_cells, ghost
    type(source_type), pointer :: s
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

       call fson_get_mpi(json, "source", sources)
       num_sources = fson_value_count_mpi(sources, ".")

       do isrc = 1, num_sources

          write(istr, '(i0)') isrc - 1
          srcstr = 'source[' // trim(istr) // '].'
          src => fson_value_get_mpi(sources, isrc)
          call fson_get_mpi(src, "name", default_name, name)
          call fson_get_mpi(src, "cell", val = cell)

          if (fson_has_mpi(src, "rate")) then
             call fson_get_mpi(src, "rate", val = rate)
             can_inject = (rate > 0._dp)
          else
             rate = default_source_rate
             can_inject = PETSC_TRUE
          end if

          call fson_get_mpi(src, "component", default_source_component, &
               injection_component, logfile, trim(srcstr) // "component")
          if (fson_has_mpi(src, "production_component")) then
             call fson_get_mpi(src, "production_component", &
                  val = production_component)
          else
             if ((.not. isothermal) .and. (injection_component == np)) then
                production_component = injection_component
             else
                production_component = default_source_production_component
             end if
          end if

          if (can_inject .and. (injection_component < np)) then
             call fson_get_mpi(src, "enthalpy", default_source_injection_enthalpy, &
                  enthalpy, logfile, trim(srcstr) // "enthalpy")
          else
             enthalpy = 0._dp
          end if

          call DMGetStratumIS(dm, cell_order_label_name, &
               cell, cell_IS, ierr); CHKERRQ(ierr)
          if (cell_IS /= 0) then

             call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
             num_cells = size(cells)
             if (num_cells > 1) then
                num_cells = 1
                if (present(logfile)) then
                   call logfile%write(LOG_LEVEL_WARN, 'source', &
                        'multiple cells found for source: ' // trim(name))
                end if
             end if
             do i = 1, num_cells
                c = cells(i)
                call DMLabelGetValue(ghost_label, c, ghost, ierr)
                CHKERRQ(ierr)
                if (ghost < 0) then
                   allocate(s)
                   call s%init(cell, c, np, rate, enthalpy, &
                        injection_component, production_component)
                   call sources_list%append(s, name)
                end if
             end do
             call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)

          end if
          call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)

       end do

    else
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_INFO, "input", "no_sources")
       end if
    end if

  end subroutine setup_sources

!------------------------------------------------------------------------

end module source_setup_module
