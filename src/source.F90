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
  use source_network_module
  use fluid_module
  use thermodynamics_module
  use hdf5io_module, only: max_field_name_length

  implicit none
  private

  PetscInt, parameter, public :: default_source_component = 0
  PetscInt, parameter, public :: default_source_injection_component = 1
  PetscInt, parameter, public :: default_source_production_component = 0
  PetscReal, parameter, public :: default_source_rate = 0._dp
  PetscReal, parameter, public :: default_source_injection_enthalpy = 83.9e3

  PetscInt, parameter, public :: num_source_scalar_variables = 13
  PetscInt, parameter, public :: max_source_variable_name_length = 24
  character(max_source_variable_name_length), public :: &
       source_scalar_variable_names(num_source_scalar_variables) = [ &
       "rate                ", "enthalpy            ", &
       "separator_pressure  ", &
       "ref_water_enthalpy  ", "ref_steam_enthalpy  ", &
       "steam_fraction      ", &
       "water_rate          ", "water_enthalpy      ", &
       "steam_rate          ", "steam_enthalpy      ", &
       "source_index        ", "natural_cell_index  ", &
       "component           "]
  PetscInt, parameter, public :: num_source_array_variables = 3
  character(max_source_variable_name_length), public :: &
       source_array_variable_names(num_source_array_variables) = [ &
       "flow          ", &
       "injection_rate", &
       "flow          "]
  character(max_field_name_length), parameter, public :: required_output_source_fields(0) = [&
       character(max_field_name_length)::]
  character(max_field_name_length), parameter, public :: default_output_source_fields(3) = [&
       "component         ",  "rate              ", "enthalpy          "]

  character(len = 6), public :: source_label_name = "source" !! Name of DMLabel for identifying source locations

  type, public, extends(source_network_node_type) :: source_type
     !! Type for mass / energy source, applying specified values of
     !! generation to each equation in a particular cell at the
     !! current time.
     private
     PetscInt, public :: local_source_index !! Index of source in local part of source vector
     PetscInt, public :: local_cell_index !! Local index of cell the source is in
     PetscReal, public :: injection_enthalpy !! Enthalpy to apply for injection
     PetscInt, public :: injection_component !! Component for injection
     PetscInt, public :: production_component !! Component for production (default 0 means all)
     PetscInt, public :: dof !! Number of degrees of freedom
     PetscInt, public :: num_primary_variables !! Number of primary thermodynamic variables
     PetscInt, public :: num_tracers !! Number of tracers
     PetscBool, public :: isothermal !! Whether equation of state is isothermal
     PetscReal, pointer, public :: source_index !! Index of source in input
     PetscReal, pointer, public :: natural_cell_index !! Natural index of cell the source is in
     PetscReal, pointer, public :: component !! Mass (or energy) component being produced or injected
     PetscReal, pointer, contiguous, public :: flow(:) !! Flows in each mass and energy component
     type(fluid_type), public :: fluid !! Fluid properties in cell (for production)
     PetscReal, pointer, contiguous, public :: tracer_injection_rate(:) !! Tracer injection rates
     PetscReal, pointer, contiguous, public :: tracer_flow(:) !! Tracer flow rates
   contains
     private
     procedure :: update_injection_mass_flow => source_update_injection_mass_flow
     procedure :: update_production_mass_flow => source_update_production_mass_flow
     procedure :: update_energy_flow => source_update_energy_flow
     procedure :: update_component_production => source_update_component_production
     procedure :: update_component_injection => source_update_component_injection
     procedure, public :: init => source_init
     procedure, public :: assign => source_assign
     procedure, public :: assign_fluid_local => source_assign_fluid_local
     procedure, public :: assign_fluid => source_assign_fluid
     procedure, public :: setup => source_setup
     procedure, public :: update_component => source_update_component
     procedure, public :: update_flow => source_update_flow
     procedure, public :: update_tracer_flow => source_update_tracer_flow
     procedure, public :: destroy => source_destroy
  end type source_type

contains

!------------------------------------------------------------------------

  subroutine source_init(self, name, eos, local_source_index, &
       local_cell_index, injection_enthalpy, injection_component, &
       production_component, num_tracers)
    !! Initialises a source object.

    use eos_module, only: eos_type

    class(source_type), intent(in out) :: self
    character(*), intent(in) :: name !! Source name
    class(eos_type), intent(in) :: eos !! Equation of state
    PetscInt, intent(in) :: local_source_index !! Source index on local process
    PetscInt, intent(in) :: local_cell_index !! Cell index on local process
    PetscReal, intent(in) :: injection_enthalpy !! Enthalpy for injection
    PetscInt, intent(in) :: injection_component !! Component for injection
    PetscInt, intent(in) :: production_component !! Component for production
    PetscInt, intent(in), optional :: num_tracers !! Number of tracers

    self%name = name
    self%num_primary_variables = eos%num_primary_variables
    call self%fluid%init(eos%num_components, eos%num_phases)
    self%local_source_index = local_source_index
    self%local_cell_index = local_cell_index
    self%injection_enthalpy = injection_enthalpy
    self%injection_component = injection_component
    self%production_component = production_component
    self%isothermal = eos%isothermal

    if (present(num_tracers)) then
       self%num_tracers = num_tracers
    else
       self%num_tracers = 0
    end if

    self%dof = num_source_scalar_variables + self%num_primary_variables + &
         self%num_tracers * 2
    self%heat = (self%production_component == self%num_primary_variables)

  end subroutine source_init

!------------------------------------------------------------------------

  subroutine source_assign(self, data, offset)
    !! Assigns pointers in source object to elements in the data
    !! array, starting from the specified offset.

    class(source_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! source data array
    PetscInt, intent(in) :: offset  !! source array offset
    ! Locals:
    PetscInt :: iflow_start, iflow_end

    call self%source_network_node_type%assign(data, offset)

    self%source_index => data(offset + 8)
    self%natural_cell_index => data(offset + 9)
    self%component => data(offset + 10)
    iflow_start = offset + 11
    iflow_end = iflow_start + self%num_primary_variables - 1
    self%flow => data(iflow_start: iflow_end)
    if (self%num_tracers > 0) then
       self%tracer_injection_rate => data(iflow_end + 1: &
            iflow_end + self%num_tracers)
       self%tracer_flow => data(iflow_end + self%num_tracers + 1: &
            iflow_end + self%num_tracers * 2)
    else
       self%tracer_injection_rate => null()
       self%tracer_flow => null()
    end if

  end subroutine source_assign

!------------------------------------------------------------------------

  subroutine source_assign_fluid_local(self, local_fluid_data, local_fluid_section)
    !! Updates fluid object from given data array and local section.

    use dm_utils_module, only: section_offset

    class(source_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    ! Locals:
    PetscInt :: fluid_offset

    fluid_offset = section_offset(local_fluid_section, self%local_cell_index)
    call self%fluid%assign(local_fluid_data, fluid_offset)

  end subroutine source_assign_fluid_local

!------------------------------------------------------------------------

  subroutine source_assign_fluid(self, fluid_data, fluid_section, fluid_range_start)
    !! Updates fluid object from given data array and global section.

    use dm_utils_module, only: global_section_offset

    class(source_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    ! Locals:
    PetscInt :: fluid_offset

    fluid_offset = global_section_offset(fluid_section, &
         self%local_cell_index, fluid_range_start)
    call self%fluid%assign(fluid_data, fluid_offset)

  end subroutine source_assign_fluid

!------------------------------------------------------------------------

  subroutine source_setup(self, source_index, natural_cell_index, rate, &
       tracer_injection_rate, separator_pressure, thermo)
    !! Sets up main parameters of a source object.

    class(source_type), intent(in out) :: self
    PetscInt, intent(in) :: source_index !! index of source in input
    PetscInt, intent(in) :: natural_cell_index !! natural index of cell the source is in
    PetscReal, intent(in) :: rate !! source flow rate
    PetscReal, intent(in) :: tracer_injection_rate(:) !! tracer injection rates
    PetscReal, intent(in) :: separator_pressure !! Separator pressure (-1 for no separator)
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics

    self%source_index = source_index
    self%natural_cell_index = natural_cell_index
    if (self%num_tracers > 0) then
       self%tracer_injection_rate = tracer_injection_rate
    end if
    if (separator_pressure > 0._dp) then
       call self%separator%init(separator_pressure, thermo)
    end if
    call self%set_rate(rate)

  end subroutine source_setup

!------------------------------------------------------------------------

  subroutine source_destroy(self)
    !! Destroys a source object.

    class(source_type), intent(in out) :: self

    call self%source_network_node_type%destroy()
    self%source_index => null()
    self%natural_cell_index => null()
    self%component => null()
    self%flow => null()
    self%tracer_injection_rate => null()
    self%tracer_flow => null()

    call self%fluid%destroy()

  end subroutine source_destroy

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

    if (nint(self%component) > 0) then
       self%enthalpy = self%injection_enthalpy
       self%flow(nint(self%component)) = self%rate
    end if
    
  end subroutine source_update_injection_mass_flow

!------------------------------------------------------------------------

  subroutine source_update_production_mass_flow(self)
    !! Updates the mass components of the flow array (and the
    !! enthalpy) for production. Only to be called if self%rate <
    !! 0. Should always be preceded by a call to assign_fluid_local().

    use fluid_module, only: fluid_type

    class(source_type), intent(in out) :: self
    ! Locals:
    PetscReal :: phase_flow_fractions(self%fluid%num_phases)

    self%flow = 0._dp
    self%enthalpy = 0._dp

    associate(np => size(self%flow), nc => self%fluid%num_components)

      if (nint(self%component) < np) then
         phase_flow_fractions = self%fluid%phase_flow_fractions()
         if (.not. self%isothermal) then
            self%enthalpy = self%fluid%specific_enthalpy( &
                 phase_flow_fractions)
         end if
      end if
      
      if (nint(self%component) <= 0) then
         ! distribute production over all mass components:
         self%flow(1: nc) = self%rate * &
              self%fluid%component_flow_fractions(phase_flow_fractions)
      else
         ! produce only specified mass or energy component:
         self%flow(nint(self%component)) = self%rate
      end if

    end associate

  end subroutine source_update_production_mass_flow

!------------------------------------------------------------------------

  subroutine source_update_energy_flow(self)
    !! Updates energy flow from mass production or injection.

    class(source_type), intent(in out) :: self

    associate(np => self%num_primary_variables)
      if (nint(self%component) < np) then
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
       call self%assign_fluid_local(local_fluid_data, local_fluid_section)
       call self%update_production_mass_flow()
    end if

    if (.not.(self%isothermal)) then
       call self%update_energy_flow()
    end if

  end subroutine source_update_flow

!------------------------------------------------------------------------

  subroutine source_update_tracer_flow(self, fluid_data, fluid_section, &
       fluid_range_start, tracer_data, tracer_section, &
       tracer_range_start, tracer_phase_index)
    !! Updates tracer flow rates from given global fluid and tracer
    !! arrays.

    use dm_utils_module, only: global_section_offset

    class(source_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    PetscReal, pointer, contiguous, intent(in) :: tracer_data(:)
    PetscSection, intent(in) :: tracer_section
    PetscInt, intent(in) :: tracer_range_start
    PetscInt, intent(in) :: tracer_phase_index(:)
    ! Locals:
    PetscInt :: tracer_offset, i
    PetscReal, pointer, contiguous :: tracer_mass_fraction(:)
    PetscReal :: phase_flow_fractions(self%fluid%num_phases)

    if (self%rate >= 0._dp) then
       self%tracer_flow = self%tracer_injection_rate
    else
       call self%assign_fluid(fluid_data, fluid_section, &
            fluid_range_start)
       tracer_offset = global_section_offset(tracer_section, &
            self%local_cell_index, tracer_range_start)
       tracer_mass_fraction => tracer_data(tracer_offset: &
            tracer_offset + self%num_tracers - 1)
       phase_flow_fractions = self%fluid%phase_flow_fractions()
       do i = 1, self%num_tracers
          self%tracer_flow(i) = tracer_mass_fraction(i) * &
               phase_flow_fractions(tracer_phase_index(i)) * &
               self%rate
       end do
    end if

  end subroutine source_update_tracer_flow

!------------------------------------------------------------------------

end module source_module
