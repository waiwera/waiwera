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
  use source_network_node_module
  use separator_module, only: num_separator_variables, separator_variable_names
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

  PetscInt, parameter, public :: num_source_scalar_variables = &
       num_source_network_node_variables + num_separator_variables + 3
  PetscInt, parameter, public :: max_source_variable_name_length = 24
  character(max_source_variable_name_length), parameter, public :: &
       source_scalar_variable_names(num_source_scalar_variables) = [ &
       source_network_variable_names, separator_variable_names, [ &
       "source_index        ", "natural_cell_index  ", &
       "component           "]]
  PetscInt, parameter, public :: num_source_constant_integer_variables = 2
  character(max_source_variable_name_length), parameter, public :: &
       source_constant_integer_variables(num_source_constant_integer_variables) = [ &
       "source_index        ", "natural_cell_index  "]
  PetscInt, parameter, public :: num_source_array_variables = 3
  character(max_source_variable_name_length), public :: &
       source_array_variable_names(num_source_array_variables) = [ &
       "flow          ", &
       "injection_rate", &
       "flow          "]
  character(max_field_name_length), parameter, public :: required_output_source_fields(0) = [&
       character(max_field_name_length)::]
  character(max_field_name_length), parameter, public :: default_output_source_fields(4) = [&
       "natural_cell_index", "component         ", &
       "rate              ", "enthalpy          "]

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
     PetscBool, public :: rate_specified !! Whether rate specified, by value or controls
     PetscReal, public :: specified_rate !! Rate specified (may be overridden by reinjectors)
     PetscBool, public :: enthalpy_specified !! Whether enthalpy specified, by value or controls
     PetscReal, public :: specified_enthalpy !! Enthalpy specified (may be overridden by reinjectors)
     PetscInt, public :: production_component !! Component for production (default 0 means all)
     PetscInt, public :: dof !! Number of degrees of freedom
     PetscInt, public :: num_primary_variables !! Number of primary thermodynamic variables
     PetscInt, public :: num_tracers !! Number of tracers
     PetscBool, public :: isothermal !! Whether equation of state is isothermal
     PetscBool, public :: fluid_dependent !! Whether flow rate / enthalpy depend on fluid state
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
     procedure, public :: init_data => source_init_data
     procedure, public :: set_rate => source_set_rate
     procedure, public :: specified_injection_rate => source_specified_injection_rate
     procedure, public :: set_enthalpy => source_set_enthalpy
     procedure, public :: update_component => source_update_component
     procedure, public :: update_flow => source_update_flow
     procedure, public :: update_tracer_flow => source_update_tracer_flow
     procedure, public :: destroy => source_destroy
  end type source_type

  type, public :: source_dependency_type
     !! Type for dependencies between sources.
     private
     PetscInt, public :: equation !! Global (block) equation index
     PetscInt, public :: cell !! Global (block) cell index
  end type source_dependency_type

contains

!------------------------------------------------------------------------

  subroutine source_init(self, name, eos, local_source_index, &
       local_cell_index, injection_enthalpy, injection_component, &
       production_component, rate_specified, specified_rate, &
       enthalpy_specified, specified_enthalpy, num_tracers)
    !! Initialises a source object. Only values stored in the object
    !! itself are initialised, not those in the source data vector
    !! accesssed via pointers.

    use eos_module, only: eos_type

    class(source_type), intent(in out) :: self
    character(*), intent(in) :: name !! Source name
    class(eos_type), intent(in) :: eos !! Equation of state
    PetscInt, intent(in) :: local_source_index !! Source index on local process
    PetscInt, intent(in) :: local_cell_index !! Cell index on local process
    PetscReal, intent(in) :: injection_enthalpy !! Enthalpy for injection
    PetscInt, intent(in) :: injection_component !! Component for injection
    PetscInt, intent(in) :: production_component !! Component for production
    PetscBool, intent(in) :: rate_specified !! Whether rate specified
    PetscReal, intent(in) :: specified_rate !! Specified rate
    PetscBool, intent(in) :: enthalpy_specified !! Whether enthalpy specified
    PetscReal, intent(in) :: specified_enthalpy !! Specified enthalpy
    PetscInt, intent(in), optional :: num_tracers !! Number of tracers

    self%name = name
    self%num_primary_variables = eos%num_primary_variables
    call self%fluid%init(eos%num_components, eos%num_mobile_phases)
    self%local_source_index = local_source_index
    self%local_cell_index = local_cell_index
    self%injection_enthalpy = injection_enthalpy
    self%injection_component = injection_component
    self%rate_specified = rate_specified
    self%specified_rate = specified_rate
    self%enthalpy_specified = enthalpy_specified
    self%specified_enthalpy = specified_enthalpy
    self%production_component = production_component
    self%isothermal = eos%isothermal
    self%fluid_dependent = PETSC_FALSE ! default - can be overridden by source controls

    if (present(num_tracers)) then
       self%num_tracers = num_tracers
    else
       self%num_tracers = 0
    end if

    self%dof = num_source_scalar_variables + self%num_primary_variables + &
         self%num_tracers * 2
    self%heat = (self%production_component == self%num_primary_variables)
    self%link_index = -1

  end subroutine source_init

!------------------------------------------------------------------------

  subroutine source_assign(self, data, offset)
    !! Assigns pointers in source object to elements in the data
    !! array, starting from the specified offset.

    class(source_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! source data array
    PetscInt, intent(in) :: offset  !! source array offset
    ! Locals:
    PetscInt :: iflow_start, iflow_end, source_offset

    call self%source_network_node_type%assign(data, offset)

    source_offset = offset + num_source_network_node_variables + &
         num_separator_variables

    self%source_index => data(source_offset)
    self%natural_cell_index => data(source_offset + 1)
    self%component => data(source_offset + 2)
    iflow_start = source_offset + 3
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

  subroutine source_init_data(self, source_index, natural_cell_index, rate, &
       tracer_injection_rate, separator_pressure, thermo)
    !! Sets up main data stored in the source data vector. The source
    !! assign() method must be called first.

    class(source_type), intent(in out) :: self
    PetscInt, intent(in) :: source_index !! index of source in input
    PetscInt, intent(in) :: natural_cell_index !! natural index of cell the source is in
    PetscReal, intent(in) :: rate !! source flow rate
    PetscReal, intent(in) :: tracer_injection_rate(:) !! tracer injection rates
    PetscReal, intent(in) :: separator_pressure(:) !! Separator pressures ([-1] for no separator)
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics

    self%source_index = source_index
    self%natural_cell_index = natural_cell_index
    if (self%num_tracers > 0) then
       self%tracer_injection_rate = tracer_injection_rate
    end if
    call self%separator%init(separator_pressure, thermo)
    call self%set_rate(rate)

  end subroutine source_init_data

!------------------------------------------------------------------------

  subroutine source_set_rate(self, rate)
    !! Sets source flow rate to specified value. For injection, also
    !! sets specified injection rate.

    class(source_type), intent(in out) :: self
    PetscReal, intent(in) :: rate !! Flow rate

    call self%source_network_node_type%set_rate(rate)
    if (self%rate_specified) then
       self%specified_rate = rate
    end if

  end subroutine source_set_rate

!------------------------------------------------------------------------

  PetscReal function source_specified_injection_rate(self) result(rate)
    !! Returns specified rate if rate is specified, or -1 otherwise.

    class(source_type), intent(in) :: self

    if (self%rate_specified) then
       rate = self%specified_rate
    else
       rate = -1._dp
    end if

  end function source_specified_injection_rate

!------------------------------------------------------------------------

  subroutine source_set_enthalpy(self, enthalpy)
    !! Sets source injection enthalpy and specified enthalpy to specified value.

    class(source_type), intent(in out) :: self
    PetscReal, intent(in) :: enthalpy !! Injection enthalpy

    self%injection_enthalpy = enthalpy
    if (self%enthalpy_specified) then
       self%specified_enthalpy = enthalpy
    end if

  end subroutine source_set_enthalpy

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

    if ((self%rate > 0._dp) .or. (self%natural_cell_index < 0)) then
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
