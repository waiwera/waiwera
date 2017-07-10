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

module source_module
 !! Module for handling sinks and sources.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fluid_module

  implicit none
  private

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
     type(fluid_type), public :: fluid !! Fluid properties in cell (for production)
     PetscBool :: fluid_updated !! Whether fluid has been updated
     PetscBool :: isothermal !! Whether equation of state is isothermal
   contains
     private
     procedure :: update_injection_mass_flow => source_update_injection_mass_flow
     procedure :: update_production_mass_flow => source_update_production_mass_flow
     procedure :: update_energy_flow => source_update_energy_flow
     procedure, public :: init => source_init
     procedure, public :: update_fluid => source_update_fluid
     procedure :: update_component_production => source_update_component_production
     procedure :: update_component_injection => source_update_component_injection
     procedure, public :: update_component => source_update_component
     procedure, public :: update_flow => source_update_flow
     procedure, public :: destroy => source_destroy
  end type source_type

contains

!------------------------------------------------------------------------

  subroutine source_init(self, cell_natural_index, cell_index, &
       eos, rate, injection_enthalpy, &
       injection_component, production_component)
    !! Initialises a source object.

    use eos_module, only: eos_type

    class(source_type), intent(in out) :: self
    PetscInt, intent(in) :: cell_natural_index !! natural index of cell the source is in
    PetscInt, intent(in) :: cell_index !! local index of cell the source is in
    class(eos_type), intent(in) :: eos !! Equation of state
    PetscReal, intent(in) :: rate !! source flow rate
    PetscReal, intent(in) :: injection_enthalpy !! enthalpy for injection
    PetscInt, intent(in) :: injection_component !! mass (or energy) component for injection
    PetscInt, intent(in) :: production_component !! mass (or energy) component for production

    self%cell_natural_index = cell_natural_index
    self%cell_index = cell_index
    allocate(self%flow(eos%num_primary_variables))
    call self%fluid%init(eos%num_components, eos%num_phases)
    self%fluid_updated = PETSC_FALSE
    self%isothermal = eos%isothermal
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
    call self%fluid%destroy()

  end subroutine source_destroy

!------------------------------------------------------------------------

  subroutine source_update_fluid(self, local_fluid_data, local_fluid_section)
    !! Updates fluid object from given data array and section, and
    !! calculates the fluid phase flow fractions.

    use dm_utils_module, only: section_offset

    class(source_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    ! Locals:
    PetscInt :: fluid_offset
    PetscErrorCode :: ierr

    if (.not. self%fluid_updated) then

       call section_offset(local_fluid_section, self%cell_index, &
            fluid_offset, ierr); CHKERRQ(ierr)
       call self%fluid%assign(local_fluid_data, fluid_offset)

       self%fluid_updated = PETSC_TRUE

    end if

  end subroutine source_update_fluid

!------------------------------------------------------------------------

  subroutine source_update_component_production(self)
    !! Update the component being produced.

    class(source_type), intent(in out) :: self
    
    if (self%production_component <= 0) then
       self%component = default_source_production_component
    else
       self%component = self%production_component
    end if

  end subroutine source_update_component_production

!------------------------------------------------------------------------

  subroutine source_update_component_injection(self)
    !! Update the component being injected.

    class(source_type), intent(in out) :: self

    if (self%injection_component <= 0) then
       self%component = default_source_injection_component
    else
       self%component = self%injection_component
    end if

  end subroutine source_update_component_injection

!------------------------------------------------------------------------

  subroutine source_update_component(self)
    !! Updates the component being injected or produced.

    class(source_type), intent(in out) :: self

    if (self%rate >= 0._dp) then
       call self%update_component_injection()
    else
       call self%update_component_production()
    end if

  end subroutine source_update_component

!------------------------------------------------------------------------

  subroutine source_update_injection_mass_flow(self)
    !! Updates the mass components of the flow array (and the
    !! enthalpy) for injection. Only to be called if self%rate >= 0.

    class(source_type), intent(in out) :: self

    self%flow = 0._dp

    if (self%component > 0) then
       self%enthalpy = self%injection_enthalpy
       self%flow(self%component) = self%rate
    end if
    
  end subroutine source_update_injection_mass_flow

!------------------------------------------------------------------------

  subroutine source_update_production_mass_flow(self)
    !! Updates the mass components of the flow array (and the
    !! enthalpy) for production. Only to be called if self%rate <
    !! 0. Should always be preceded by a call to update_fluid().

    use fluid_module, only: fluid_type

    class(source_type), intent(in out) :: self
    ! Locals:
    PetscReal :: phase_flow_fractions(self%fluid%num_phases)

    self%flow = 0._dp
    self%enthalpy = 0._dp

    associate(np => size(self%flow), nc => self%fluid%num_components)

      if (self%component < np) then
         phase_flow_fractions = self%fluid%phase_flow_fractions()
         if (.not. self%isothermal) then
            self%enthalpy = self%fluid%specific_enthalpy( &
                 phase_flow_fractions)
         end if
      end if
      
      if (self%component <= 0) then
         ! distribute production over all mass components:
         self%flow(1: nc) = self%rate * &
              self%fluid%component_flow_fractions(phase_flow_fractions)
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

  subroutine source_update_flow(self, local_fluid_data, local_fluid_section)
    !! Updates the flow array, according to the source parameters and
    !! fluid conditions.

    use fluid_module, only: fluid_type

    class(source_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section

    if (self%rate >= 0._dp) then
       call self%update_component_injection()
       call self%update_injection_mass_flow()
    else
       call self%update_component_production()
       call self%update_fluid(local_fluid_data, local_fluid_section)
       call self%update_production_mass_flow()
    end if

    if (.not.(self%isothermal)) then
       call self%update_energy_flow()
    end if

    self%fluid_updated = PETSC_FALSE

  end subroutine source_update_flow

!------------------------------------------------------------------------

end module source_module
