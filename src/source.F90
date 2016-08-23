module source_module
 !! Module for handling sinks and sources.

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  public :: setup_source_vector

contains

!------------------------------------------------------------------------

  subroutine setup_source_vector(json, dm, np, isothermal, &
       source, range_start, logfile)
    !! Sets up sinks and sources. Source strengths are stored (for
    !! now) in the source vector, with values for all components in
    !! all cells.

    use kinds_module
    use fson
    use fson_mpi_module
    use logfile_module
    use mesh_module, only: cell_order_label_name
    use dm_utils_module, only: global_section_offset

    type(fson_value), pointer, intent(in) :: json !! JSON file object
    DM, intent(in) :: dm !! Mesh DM 
    PetscInt, intent(in) :: np !! Number of primary variables
    PetscBool, intent(in) :: isothermal !! Whether EOS is isothermal
    PetscInt, intent(in) :: range_start !! Range start for global source vector
    Vec, intent(in out) :: source !! Global source vector
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: c, isrc, i, component
    type(fson_value), pointer :: sources, src
    PetscInt :: num_sources, cell, offset, num_cells, ghost
    PetscInt, pointer :: cells(:)
    PetscReal :: q, enthalpy
    PetscReal, pointer, contiguous :: source_array(:), cell_source(:)
    PetscBool :: mass_inject
    PetscSection :: section
    IS :: cell_IS
    DMLabel :: ghost_label
    character(len=64) :: srcstr
    character(len=12) :: istr
    PetscInt, parameter ::  default_component = 1
    PetscReal, parameter :: default_enthalpy = 83.9e3

    call DMCreateGlobalVector(dm, source, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(source, "source", ierr); CHKERRQ(ierr)
    call VecSet(source, 0._dp, ierr); CHKERRQ(ierr)

    if (fson_has_mpi(json, "source")) then

       call DMGetDefaultGlobalSection(dm, section, ierr); CHKERRQ(ierr)
       call VecGetArrayF90(source, source_array, ierr); CHKERRQ(ierr)
       call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)

       call fson_get_mpi(json, "source", sources)
       num_sources = fson_value_count_mpi(sources, ".")
       do isrc = 1, num_sources
          write(istr, '(i0)') isrc - 1
          srcstr = 'source[' // trim(istr) // '].'
          src => fson_value_get_mpi(sources, isrc)
          call fson_get_mpi(src, "cell", val = cell)
          call fson_get_mpi(src, "component", default_component, &
               component, logfile, trim(srcstr) // "component")
          call fson_get_mpi(src, "value", val = q)
          mass_inject = ((.not.(isothermal)) .and. (component < np) &
               .and. (q > 0._dp))
          if (mass_inject) then
             call fson_get_mpi(src, "enthalpy", default_enthalpy, &
                  enthalpy, logfile, trim(srcstr) // "enthalpy")
          else
             enthalpy = 0._dp
          end if
          call DMGetStratumIS(dm, cell_order_label_name, &
               cell, cell_IS, ierr); CHKERRQ(ierr)
          if (cell_IS /= 0) then
             call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
             num_cells = size(cells)
             do i = 1, num_cells
                c = cells(i)
                call DMLabelGetValue(ghost_label, c, ghost, ierr)
                CHKERRQ(ierr)
                if (ghost < 0) then
                   call global_section_offset(section, c, range_start, &
                        offset, ierr)
                   cell_source => source_array(offset : offset + np - 1)
                   cell_source(component) = cell_source(component) + q
                   if (mass_inject) then
                      cell_source(np) = cell_source(np) + enthalpy * q
                   end if
                end if
             end do
             call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
          end if
          call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
       end do
       call VecRestoreArrayF90(source, source_array, ierr); CHKERRQ(ierr)

    else
       call logfile%write(LOG_LEVEL_INFO, "input", "no_sources")
    end if

  end subroutine setup_source_vector

!------------------------------------------------------------------------

end module source_module
