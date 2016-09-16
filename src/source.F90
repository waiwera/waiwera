module source_module
 !! Module for handling sinks and sources.

  use list_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  type, public :: source_type
     !! Type for mass / energy source, applying specified values of
     !! generation to each equation in a particular cell at the
     !! current time.
     private
     PetscReal, allocatable, public :: flow(:) !! Generation flow rate for each equation
     PetscInt, public :: cell_natural_index !! Natural index of cell the source is in
     PetscInt, public :: cell_index !! Local index of cell the source is in
   contains
     private
     procedure, public :: init => source_init
     procedure, public :: destroy => source_destroy
     ! procedure, public :: energy_production => source_energy_production
  end type source_type

  type, public, abstract :: source_control_type
     !! Abstract type for mass / energy source control, controlling
     !! generation values over time in one or more sources.
     private
     type(list_type), public :: sources
   contains
     procedure(source_control_init_procedure), public, deferred :: init
     procedure(source_control_destroy_procedure), public, deferred :: destroy
     procedure(source_control_update_procedure), public, deferred :: update
  end type source_control_type

  abstract interface

     subroutine source_control_init_procedure(self)
       !! Initialises source control object
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
     end subroutine source_control_init_procedure

     subroutine source_control_destroy_procedure(self)
       !! Destroys source control object
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
     end subroutine source_control_destroy_procedure

     subroutine source_control_update_procedure(self, time)
       !! Updates sources at the specified time.
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
       PetscReal, intent(in) :: time
     end subroutine source_control_update_procedure

  end interface

  public :: setup_source_vector

contains

!------------------------------------------------------------------------

  subroutine source_init(self, num_primary_variables, cell_natural_index)
    !! Initialises a source object.

    class(source_type), intent(in out) :: self
    PetscInt, intent(in) :: num_primary_variables
    PetscInt, intent(in) :: cell_natural_index

    allocate(self%flow(num_primary_variables))
    self%cell_natural_index = cell_natural_index
    ! Determine local cell index? or determine before calling this
    ! routine, and only init gener if cell on-process? (or in
    ! coordinated distributed source)

  end subroutine source_init

!------------------------------------------------------------------------

  subroutine source_destroy(self)
    !! Destroys a source object.
    class(source_type), intent(in out) :: self

    deallocate(self%flow)

  end subroutine source_destroy

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
