module source_module
 !! Module for handling sinks and sources.

  use kinds_module
  use list_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  PetscInt, parameter, public :: max_source_name_length = 32
  PetscInt, parameter, public :: default_source_component = 0
  PetscInt, parameter, public :: default_source_injection_component = 1
  PetscInt, parameter, public :: default_source_production_component = 0
  PetscReal, parameter, public :: default_source_rate = 0._dp
  PetscReal, parameter, public :: default_source_injection_enthalpy = 83.9e3

  type, public :: source_type
     !! Type for mass / energy source, applying specified values of
     !! generation to each equation in a particular cell at the
     !! current time.
     private
     PetscInt, public :: cell_natural_index !! Natural index of cell the source is in
     PetscInt, public :: cell_index !! Local index of cell the source is in
     PetscInt, public :: component !! Which mass (or energy) component is being produced or injected
     PetscReal, public :: rate !! Flow rate
     PetscReal, public :: injection_enthalpy !! Enthalpy to apply for injection
     PetscReal, allocatable, public :: flow(:) !! Flows in each mass and energy component
     PetscReal, public :: enthalpy !! Enthalpy of produced or injected fluid
   contains
     private
     procedure :: update_injection_mass_flow => source_update_injection_mass_flow
     procedure :: update_production_mass_flow => source_update_production_mass_flow
     procedure :: update_energy_flow => source_update_energy_flow
     procedure, public :: init => source_init
     procedure, public :: update_flow => source_update_flow
     procedure, public :: destroy => source_destroy
  end type source_type

  type, public, abstract :: source_control_type
     !! Abstract type for source control, controlling source
     !! parameters over time, for one or more sources.
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

  public :: setup_sources

contains

!------------------------------------------------------------------------

  subroutine source_init(self, cell_natural_index, cell_index, &
       num_primary, component, rate, injection_enthalpy)
    !! Initialises a source object.

    class(source_type), intent(in out) :: self
    PetscInt, intent(in) :: cell_natural_index !! natural index of cell the source is in
    PetscInt, intent(in) :: cell_index !! local index of cell the source is in
    PetscInt, intent(in) :: num_primary !! Number of primary variables
    PetscInt, intent(in) :: component !! mass (or energy) component the source is applied to
    PetscReal, intent(in) :: rate !! source flow rate
    PetscReal, intent(in) :: injection_enthalpy !! enthalpy for injection

    self%cell_natural_index = cell_natural_index
    self%cell_index = cell_index
    allocate(self%flow(num_primary))
    self%rate = rate
    self%component = component
    self%injection_enthalpy = injection_enthalpy

  end subroutine source_init

!------------------------------------------------------------------------

  subroutine source_destroy(self)
    !! Destroys a source object.

    class(source_type), intent(in out) :: self

    deallocate(self%flow)

  end subroutine source_destroy

!------------------------------------------------------------------------

  subroutine source_update_injection_mass_flow(self, isothermal)
    !! Updates the mass components of the flow array (and the
    !! enthalpy) for injection. Only to be called if self%rate >= 0.
    
    class(source_type), intent(in out) :: self
    PetscBool, intent(in) :: isothermal

    self%flow = 0._dp
    if (self%component > 0) then
       self%enthalpy = self%injection_enthalpy
       self%flow(self%component) = self%rate
    end if
    
  end subroutine source_update_injection_mass_flow

!------------------------------------------------------------------------

  subroutine source_update_production_mass_flow(self, isothermal, fluid)
    !! Updates the mass components of the flow array (and the
    !! enthalpy) for production. Only to be called if self%rate < 0.

    use fluid_module, only: fluid_type

    class(source_type), intent(in out) :: self
    PetscBool, intent(in) :: isothermal
    type(fluid_type), intent(in) :: fluid
    ! Locals:
    PetscReal :: flow_fractions(fluid%num_phases)

    self%flow = 0._dp

    associate(np => size(self%flow), nc => fluid%num_components)

      if (self%component < np) then
         flow_fractions = fluid%flow_fractions()
         if (.not. isothermal) then
            self%enthalpy = fluid%specific_enthalpy(flow_fractions)
         end if
      end if
      
      if (self%component <= 0) then
         ! distribute production over all mass components:
         self%flow(1: nc) = self%rate * &
              fluid%component_flow_fractions(flow_fractions)
      else
         ! produce only specified mass or energy component:
         self%flow(self%component) = self%rate
      end if

    end associate

  end subroutine source_update_production_mass_flow

!------------------------------------------------------------------------

  subroutine source_update_energy_flow(self)
    !! Updates energy flow from mass production or injection.

    class(source_type), intent(in out) :: self

    associate(np => size(self%flow))
      if (self%component < np) then
         self%flow(np) = self%flow(np) + self%enthalpy * self%rate
      end if
    end associate

  end subroutine source_update_energy_flow

!------------------------------------------------------------------------

  subroutine source_update_flow(self, fluid, isothermal)
    !! Updates the flow array, according to the source parameters and
    !! fluid conditions.

    use fluid_module, only: fluid_type

    class(source_type), intent(in out) :: self
    type(fluid_type), intent(in) :: fluid
    PetscBool, intent(in) :: isothermal

    if (self%rate >= 0._dp) then
       call self%update_injection_mass_flow(isothermal)
    else
       call self%update_production_mass_flow(isothermal, fluid)
    end if

    if (.not.(isothermal)) then
       call self%update_energy_flow()
    end if

  end subroutine source_update_flow

!------------------------------------------------------------------------

  subroutine setup_sources(json, dm, np, isothermal, sources_list, logfile)
    !! Sets up list of sinks and sources.

    use kinds_module
    use fson
    use fson_mpi_module
    use logfile_module
    use mesh_module, only: cell_order_label_name

    type(fson_value), pointer, intent(in) :: json !! JSON file object
    DM, intent(in) :: dm !! Mesh DM 
    PetscInt, intent(in) :: np !! Number of primary variables
    PetscBool, intent(in) :: isothermal !! Whether EOS is isothermal
    type(list_type), intent(in out) :: sources_list !! List of sources
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: c, isrc, i, component
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
    PetscReal, parameter :: default_rate = 0._dp
    PetscInt, parameter ::  default_injection_component = 1
    PetscReal, parameter :: default_enthalpy = 83.9e3
    ! for production, default is to produce all components:
    PetscInt, parameter ::  default_production_component = 0

    call sources_list%init(delete_deallocates = PETSC_TRUE)

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
          call fson_get_mpi(src, "rate", default_rate, rate, logfile, &
               trim(srcstr) // "rate")
          if (rate >= 0._dp) then
             call fson_get_mpi(src, "component", default_injection_component, &
                  component, logfile, trim(srcstr) // "component")
             if (component < np) then
                call fson_get_mpi(src, "enthalpy", default_enthalpy, &
                     enthalpy, logfile, trim(srcstr) // "enthalpy")
             else
                enthalpy = 0._dp
             end if
          else
             call fson_get_mpi(src, "component", default_production_component, &
                  component, logfile, trim(srcstr) // "component")
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
                   call s%init(cell, c, np, component, rate, enthalpy)
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

end module source_module
