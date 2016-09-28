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
     PetscInt, public :: component !! Mass (or energy) component being produced or injected
     PetscInt, public :: injection_component !! Component for injection
     PetscInt, public :: production_component !! Component for production (default 0 means all)
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

contains

!------------------------------------------------------------------------

  subroutine source_init(self, cell_natural_index, cell_index, &
       num_primary, rate, injection_enthalpy, &
       injection_component, production_component)
    !! Initialises a source object.

    class(source_type), intent(in out) :: self
    PetscInt, intent(in) :: cell_natural_index !! natural index of cell the source is in
    PetscInt, intent(in) :: cell_index !! local index of cell the source is in
    PetscInt, intent(in) :: num_primary !! Number of primary variables
    PetscReal, intent(in) :: rate !! source flow rate
    PetscReal, intent(in) :: injection_enthalpy !! enthalpy for injection
    PetscInt, intent(in) :: injection_component !! mass (or energy) component for injection
    PetscInt, intent(in) :: production_component !! mass (or energy) component for production

    self%cell_natural_index = cell_natural_index
    self%cell_index = cell_index
    allocate(self%flow(num_primary))
    self%rate = rate
    self%injection_enthalpy = injection_enthalpy
    self%injection_component = injection_component
    self%production_component = production_component

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

    if (self%injection_component <= 0) then
       self%component = default_source_injection_component
    else
       self%component = self%injection_component
    end if

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

    if (self%production_component <= 0) then
       self%component = default_source_production_component
    else
       self%component = self%production_component
    end if

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

end module source_module
