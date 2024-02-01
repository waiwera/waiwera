!   Copyright 2022 University of Auckland.

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

module source_network_module
  !! Module for source networks, including sources, controls, groups
  !! and reinjectors.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use list_module
  use source_module
  use source_network_group_module
  use source_network_reinjector_module
  use control_module
  use source_control_module

  implicit none
  private

  type, public :: source_network_type
     !! Type for source network.
     private
     type(list_type), public :: sources !! List of source objects
     type(list_type), public :: groups !! Source network groups
     type(list_type), public :: separated_sources !! Sources with separators
     type(list_type), public :: source_controls !! Controls on sources/sinks
     type(list_type), public :: network_controls !! Controls on source network nodes
     type(list_type), public :: reinjectors !! Reinjectors distributing flows to injection sources
     Vec, public :: source !! Vector for source/sink data
     Vec, public :: group !! Vector for source network group data
     Vec, public :: reinjector !! Vector for source network reinjector data
     PetscInt, public :: source_range_start !! Range start for source vector
     PetscInt, public :: group_range_start !! Range start for source group vector
     PetscInt, public :: reinjector_range_start !! Range start for source reinjector vector
     PetscInt, public :: num_sources !! Total number of source/sink terms on all processes
     PetscInt, public :: num_groups !! Total number of source network groups on all processes
     PetscInt, public :: num_reinjectors !! Total number of source network reinjectors on all processes
     IS, public :: source_index !! Index set defining natural to global source ordering
     IS, public :: group_index !! Index set defining natural to global source network group ordering
     IS, public :: reinjector_index !! Index set defining natural to global source network reinjector ordering
   contains
     private
     procedure, public :: init => source_network_init
     procedure, public :: update => source_network_update
     procedure, public :: assemble_cell_inflows => source_network_assemble_cell_inflows
     procedure, public :: destroy => source_network_destroy
  end type source_network_type

contains

!------------------------------------------------------------------------

  subroutine source_network_init(self)
    !! Initialises source network. Only the main lists of source
    !! network objects are initialised here.

    class(source_network_type), intent(in out) :: self

       call self%sources%init(owner = PETSC_TRUE)
       call self%groups%init(owner = PETSC_TRUE)
       call self%separated_sources%init(owner = PETSC_FALSE)
       call self%source_controls%init(owner = PETSC_TRUE)
       call self%network_controls%init(owner = PETSC_TRUE)
       call self%reinjectors%init(owner = PETSC_TRUE)

  end subroutine source_network_init

!------------------------------------------------------------------------

  subroutine source_network_update(self, t, interval, fluid_data, fluid_section, &
       update_reinjection)
    !! Updates flows through source network, applying controls and
    !! updating source rates. If update_reinjection is true then
    !! reinjector inputs and capacities are updated (this is not done
    !! for perturbed primary variables during Jacobian calculations,
    !! as it results in poor non-linear solver performance).

    use dm_utils_module, only: global_vec_section, global_section_offset

    class(source_network_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time (s)
    PetscReal, intent(in) :: interval(2) !! time interval bounds
    PetscReal, pointer, contiguous, intent(in out) :: fluid_data(:) !! array on fluid vector
    PetscSection, intent(in out) :: fluid_section !! fluid section
    PetscBool, intent(in) :: update_reinjection !! Whether to update reinjection inputs and capacities
    ! Locals:
    PetscSection :: source_section, group_section, reinjector_section
    PetscReal, pointer, contiguous :: source_data(:), group_data(:), reinjector_data(:)
    PetscErrorCode :: ierr

    call global_vec_section(self%source, source_section)
    call VecGetArrayF90(self%source, source_data, ierr); CHKERRQ(ierr)
    call global_vec_section(self%group, group_section)
    call VecGetArrayF90(self%group, group_data, ierr); CHKERRQ(ierr)
    call global_vec_section(self%reinjector, reinjector_section)
    call VecGetArrayF90(self%reinjector, reinjector_data, ierr); CHKERRQ(ierr)

    call self%sources%traverse(source_assign_iterator)
    call self%groups%traverse(group_assign_iterator)
    call self%reinjectors%traverse(reinjector_assign_iterator)

    call self%separated_sources%traverse(source_separator_iterator)
    call self%source_controls%traverse(control_iterator)
    call self%groups%traverse(group_iterator)
    call self%network_controls%traverse(control_iterator)
    if (update_reinjection) then
       call self%reinjectors%traverse(reinjector_capacity_iterator)
    end if
    call self%reinjectors%traverse(reinjector_iterator, backwards = PETSC_TRUE)

    call VecRestoreArrayF90(self%reinjector, reinjector_data, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%group, group_data, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%source, source_data, ierr); CHKERRQ(ierr)

  contains

!........................................................................

    subroutine source_assign_iterator(node, stopped)
      !! Assigns data pointers for all sources.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: s, source_offset

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         s = source%local_source_index
         source_offset = global_section_offset(source_section, &
              s, self%source_range_start)
         call source%assign(source_data, source_offset)
         call source%update_flow(fluid_data, fluid_section)
      end select
    end subroutine source_assign_iterator

!........................................................................

    subroutine group_assign_iterator(node, stopped)
      !! Assigns data pointers for all source network groups.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: g, group_offset

      stopped = PETSC_FALSE
      select type (group => node%data)
      class is (source_network_group_type)
         if (group%rank == 0) then
            g = group%local_group_index
            group_offset = global_section_offset( &
                 group_section, g, self%group_range_start)
            call group%assign(group_data, group_offset)
         end if
      end select

    end subroutine group_assign_iterator

!........................................................................

    subroutine reinjector_assign_iterator(node, stopped)
      !! Assigns data pointers for all source network reinjectors.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: r, reinjector_offset

      stopped = PETSC_FALSE
      select type (reinjector => node%data)
      class is (source_network_reinjector_type)
         if (reinjector%rank == 0) then
            r = reinjector%local_reinjector_index
            reinjector_offset = global_section_offset( &
                 reinjector_section, r, self%reinjector_range_start)
            call reinjector%assign(reinjector_data, reinjector_offset)
         end if
      end select

    end subroutine reinjector_assign_iterator

!........................................................................

    subroutine source_separator_iterator(node, stopped)
      !! Updates enthalpy and separated outputs from separators.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscReal, allocatable :: phase_flow_fractions(:)

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)

         allocate(phase_flow_fractions(source%fluid%num_phases))
         phase_flow_fractions = source%fluid%phase_flow_fractions()
         source%enthalpy = source%fluid%specific_enthalpy(phase_flow_fractions)
         deallocate(phase_flow_fractions)

         call source%get_separated_flows()

      end select

    end subroutine source_separator_iterator

!........................................................................

    subroutine control_iterator(node, stopped)
      !! Applies source network node controls.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (control => node%data)
      class is (integer_object_control_type)
         call control%update()
      class is (interval_update_object_control_type)
         call control%update(interval)
      class is (pressure_reference_source_control_type)
         call control%update(t, interval, fluid_data, &
              fluid_section)
      end select

    end subroutine control_iterator

!........................................................................

    subroutine group_iterator(node, stopped)
      !! Computes output for source network groups.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (group => node%data)
      class is (source_network_group_type)
         call group%sum()
      end select

    end subroutine group_iterator

!........................................................................

    subroutine reinjector_capacity_iterator(node, stopped)
      !! Calculates output capacity of a reinjector.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (reinjector => node%data)
      class is (source_network_reinjector_type)
         call reinjector%capacity()
      end select

    end subroutine reinjector_capacity_iterator

!........................................................................

    subroutine reinjector_iterator(node, stopped)
      !! Updates reinjection output flows to injection sources (or
      !! other reinjectors).

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (reinjector => node%data)
      class is (source_network_reinjector_type)
         call reinjector%distribute(update_reinjection)
      end select

    end subroutine reinjector_iterator

  end subroutine source_network_update

!------------------------------------------------------------------------

  subroutine source_network_assemble_cell_inflows(self, eos, inflow_data, &
       inflow_section, inflow_range_start, fluid_data, fluid_section, &
       cell_geom_data, cell_geom_section)
    !! Assembles contributions from sources into array on inflow vector.

    use eos_module, only: eos_type
    use dm_utils_module, only: section_offset, global_section_offset
    use cell_module, only: cell_type

    class(source_network_type), intent(in out) :: self
    class(eos_type), intent(in out) :: eos
    PetscReal, pointer, contiguous, intent(in out) :: inflow_data(:) !! array on vector to assemble cell inflows into
    PetscSection, intent(in out) :: inflow_section !! section for vector to assemble into
    PetscInt, intent(in) :: inflow_range_start !! range start for vector to assemble into
    PetscReal, pointer, contiguous, intent(in out) :: fluid_data(:) !! array on fluid vector
    PetscSection, intent(in out) :: fluid_section !! section for fluid vector
    PetscReal, pointer, contiguous, intent(in out) :: cell_geom_data(:) !! array on cell geometry vector
    PetscSection, intent(in out) :: cell_geom_section !! section for cell geometry vector
    ! Locals:
    type(cell_type) :: cell

    call cell%init(eos%num_components, eos%num_mobile_phases)
    call self%sources%traverse(source_assembly_iterator)
    call cell%destroy()

  contains

    subroutine source_assembly_iterator(node, stopped)
      !! Assembles contributions from sources to global RHS array.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: c, cell_geom_offset, inflow_offset
      PetscReal, pointer, contiguous :: inflow(:)

      stopped = PETSC_FALSE
      select type (source => node%data)
      class is (source_type)

         c = source%local_cell_index
         if (c >= 0) then

            inflow_offset = global_section_offset(inflow_section, c, &
                 inflow_range_start)
            inflow => inflow_data(inflow_offset : inflow_offset + &
                 eos%num_primary_variables - 1)

            cell_geom_offset = section_offset(cell_geom_section, c)
            call cell%assign_geometry(cell_geom_data, cell_geom_offset)

            call source%update_flow(fluid_data, fluid_section)
            inflow = inflow + source%flow / cell%volume

         end if
      end select

    end subroutine source_assembly_iterator

  end subroutine source_network_assemble_cell_inflows

!------------------------------------------------------------------------

  subroutine source_network_destroy(self)
    !! Destroys source network.

    use source_network_node_module, only: source_network_node_list_node_data_destroy
    use control_module, only: object_control_list_node_data_destroy

    class(source_network_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call VecDestroy(self%source, ierr); CHKERRQ(ierr)
    call VecDestroy(self%group, ierr); CHKERRQ(ierr)
    call VecDestroy(self%reinjector, ierr); CHKERRQ(ierr)

    call ISDestroy(self%source_index, ierr); CHKERRQ(ierr)
    call ISDestroy(self%group_index, ierr); CHKERRQ(ierr)
    call ISDestroy(self%reinjector_index, ierr); CHKERRQ(ierr)

    call self%separated_sources%destroy()
    call self%source_controls%destroy(object_control_list_node_data_destroy, &
         reverse = PETSC_TRUE)
    call self%groups%destroy( &
         source_network_node_list_node_data_destroy, reverse = PETSC_TRUE)
    call self%network_controls%destroy( &
         object_control_list_node_data_destroy, reverse = PETSC_TRUE)
    call self%reinjectors%destroy( &
         source_network_node_list_node_data_destroy, reverse = PETSC_TRUE)
    call self%sources%destroy(source_network_node_list_node_data_destroy)

  end subroutine source_network_destroy
  
!------------------------------------------------------------------------

end module source_network_module
  
