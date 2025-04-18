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
     type(list_type), public :: dependencies !! Dependencies between sources
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
     procedure, public :: identify_source_dependencies => &
          source_network_identify_source_dependencies
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
    call self%dependencies%init(owner = PETSC_TRUE)

  end subroutine source_network_init

!------------------------------------------------------------------------

  subroutine source_network_update(self, t, interval, fluid_data, fluid_section)
    !! Updates flows through source network, applying controls and
    !! updating source rates.

    use dm_utils_module, only: global_vec_section, global_section_offset

    class(source_network_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time (s)
    PetscReal, intent(in) :: interval(2) !! time interval bounds
    PetscReal, pointer, contiguous, intent(in out) :: fluid_data(:) !! array on fluid vector
    PetscSection, intent(in out) :: fluid_section !! fluid section
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
    call self%reinjectors%traverse(reinjector_capacity_iterator)
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
         call reinjector%distribute()
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

  subroutine source_network_identify_source_dependencies(self)
    !! Identifies dependencies between sources which need to be added
    !! to the Jacobian matrix.

    class(source_network_type), intent(in out) :: self
    ! Locals:

    call self%reinjectors%traverse(reinjection_production_dependency_iterator)
    call self%network_controls%traverse(group_dependency_iterator)
    call self%reinjectors%traverse(fluid_dep_reinjection_dependency_iterator)

  contains

!........................................................................

    subroutine reinjection_production_dependency_iterator(node, stopped)
      !! Adds dependencies between reinjection and production sources.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt, allocatable :: production_cell_indices(:), &
           reinjector_cell_indices(:)
      PetscInt :: ir, ip
      type(source_dependency_type), pointer :: dep

      stopped = PETSC_FALSE
      select type (reinjector => node%data)
      class is (source_network_reinjector_type)
         if (associated(reinjector%in)) then

            select type (input => reinjector%in)
            type is (source_type)
               production_cell_indices = [nint(input%natural_cell_index)]
            type is (source_network_group_type)
               production_cell_indices = input%source_cell_indices
            end select

            if (allocated(production_cell_indices)) then
               reinjector_cell_indices = [reinjector%water_source_cell_indices, &
                    reinjector%steam_source_cell_indices]
               do ip = 1, size(production_cell_indices)
                  do ir = 1, size(reinjector_cell_indices)
                     allocate(dep)
                     dep%equation = reinjector_cell_indices(ir)
                     dep%cell = production_cell_indices(ip)
                     call self%dependencies%append(dep)
                  end do
               end do
               deallocate(reinjector_cell_indices)
            end if
         end if
      end select

    end subroutine reinjection_production_dependency_iterator

!........................................................................

    subroutine group_dependency_iterator(node, stopped)
      !! Adds dependencies between sources in limited network groups.

      use source_network_control_module, only: &
           limiter_table_source_network_control_type

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: i1, i2
      type(source_dependency_type), pointer :: dep

      stopped = PETSC_FALSE
      select type (control => node%data)
      type is (limiter_table_source_network_control_type)
         select type (group => control%objects%head%data)
         class is (source_network_group_type)
            if (group%rank == 0) then
               do i1 = 1, size(group%source_cell_indices)
                  do i2 = 1, size(group%source_cell_indices)
                     if (i1 /= i2) then
                        allocate(dep)
                        dep%equation = group%source_cell_indices(i1)
                        dep%cell = group%source_cell_indices(i2)
                        call self%dependencies%append(dep)
                     end if
                  end do
               end do
            end if
         end select
      end select

    end subroutine group_dependency_iterator

!........................................................................

    subroutine fluid_dep_reinjection_dependency_iterator(node, stopped)
      !! Adds dependencies between reinjection sources and any
      !! fluid-dependent reinjection sources in their reinjectors.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: i1, i2, row, col
      type(source_dependency_type), pointer :: dep

      stopped = PETSC_FALSE
      select type (reinjector => node%data)
      class is (source_network_reinjector_type)
         if (reinjector%rank == 0) then

            do i1 = 1, size(reinjector%water_fluid_dep_source_cell_indices)
               do i2 = 1, size(reinjector%water_source_cell_indices)
                  col = reinjector%water_fluid_dep_source_cell_indices(i1)
                  row = reinjector%water_source_cell_indices(i2)
                  if (row /= col) then
                     allocate(dep)
                     dep%equation = row
                     dep%cell = col
                     call self%dependencies%append(dep)
                  end if
               end do
            end do

            do i1 = 1, size(reinjector%steam_fluid_dep_source_cell_indices)
               do i2 = 1, size(reinjector%steam_source_cell_indices)
                  col = reinjector%steam_fluid_dep_source_cell_indices(i1)
                  row = reinjector%steam_source_cell_indices(i2)
                  if (row /= col) then
                     allocate(dep)
                     dep%equation = row
                     dep%cell = col
                     call self%dependencies%append(dep)
                  end if
               end do
            end do

         end if
      end select

    end subroutine fluid_dep_reinjection_dependency_iterator

  end subroutine source_network_identify_source_dependencies

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
    call self%dependencies%destroy()

  end subroutine source_network_destroy
  
!------------------------------------------------------------------------

end module source_network_module
  
