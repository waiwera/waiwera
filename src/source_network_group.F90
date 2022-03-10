!   Copyright 2021 University of Auckland.

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

module source_network_group_module
  !! Module for source groups.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use source_network_module
  use source_module, only: source_type
  use separator_module, only: num_separator_variables, separator_variable_names
  use list_module, only: list_type, list_node_type
  use hdf5io_module, only: max_field_name_length
  use thermodynamics_module

  implicit none
  private

  PetscInt, parameter, public :: num_source_network_group_variables = &
       num_source_network_node_variables + num_separator_variables + 1
  PetscInt, parameter, public :: max_source_network_group_variable_name_length = 24
  character(max_source_network_group_variable_name_length), parameter, public :: &
       source_network_group_variable_names(num_source_network_group_variables) = [ &
       source_network_variable_names, separator_variable_names, [ &
       "group_index         "]]
  PetscInt, parameter, public :: num_source_network_group_constant_integer_variables = 1
  character(max_source_network_group_variable_name_length), parameter, public :: &
       source_network_group_constant_integer_variables( &
       num_source_network_group_constant_integer_variables) = [ &
       "group_index       "]
  character(max_field_name_length), parameter, public :: &
       required_output_source_network_group_fields(0) = [&
       character(max_field_name_length)::]
  character(max_field_name_length), parameter, public :: &
       default_output_source_network_group_fields(2) = [&
       "rate              ", "enthalpy          "]

  type, public, extends(source_network_node_type) :: source_network_group_type
     !! Type for group of source network nodes, e.g. multi-feed well
     !! or group of makeup wells.
     private
     type(list_type), public :: nodes !! List of nodes in group
     MPI_Comm :: comm !! MPI communicator for group
     PetscInt, public :: rank !! Rank of group in its own communicator
     PetscInt, public :: local_group_index !! Index of group in local part of source group vector (-1 if not a root group)
     PetscReal, pointer, public :: group_index !! Index of source group in input
   contains
     private
     procedure, public :: init => source_network_group_init
     procedure, public :: init_comm => source_network_group_init_comm
     procedure, public :: assign => source_network_group_assign
     procedure, public :: init_data => source_network_group_init_data
     procedure, public :: set_rate => source_network_group_set_rate
     procedure, public :: scale_rate => source_network_group_scale_rate
     procedure, public :: sum => source_network_group_sum
     procedure, public :: default_separated_flows => source_network_group_default_separated_flows
     procedure, public :: get_separated_flows => source_network_group_get_separated_flows
     procedure, public :: add_separated_flows => source_network_group_add_separated_flows
     procedure, public :: add_flows => source_network_group_add_flows
     procedure, public :: destroy => source_network_group_destroy
  end type source_network_group_type

contains

!------------------------------------------------------------------------

  subroutine source_network_group_init(self, name)
    !! Initialises a source network group. Only variables stored in
    !! the object itself are initialised, not those stored in the
    !! group data vector and accessed via pointers.

    class(source_network_group_type), intent(in out) :: self
    character(*), intent(in) :: name !! Group name

    self%name = name
    call self%nodes%init(owner = PETSC_FALSE)

  end subroutine source_network_group_init

!------------------------------------------------------------------------

  subroutine source_network_group_init_comm(self)
    !! Initialises MPI communicator for the group. It is assumed that
    !! the nodes list has already been populated.

    class(source_network_group_type), intent(in out) :: self
    ! Locals:
    PetscInt :: colour, group_rank
    PetscErrorCode :: ierr

    colour = MPI_UNDEFINED
    call self%nodes%traverse(group_comm_iterator)

    call MPI_comm_split(PETSC_COMM_WORLD, colour, 0, self%comm, ierr)

    if (self%comm /= MPI_COMM_NULL) then
       call MPI_COMM_RANK(self%comm, self%rank, ierr)
    else
       self%rank = -1
    end if

  contains

    subroutine group_comm_iterator(node, stopped)
      !! Stops and sets colour if there are any local sources, or any
      !! groups with this rank as root, in the group node list.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (n => node%data)
      type is (source_type)
         colour = 1
         stopped = PETSC_TRUE
      type is (source_network_group_type)
         if (n%rank == 0) then
            colour = 1
            stopped = PETSC_TRUE
         end if
      end select

    end subroutine group_comm_iterator

  end subroutine source_network_group_init_comm

!------------------------------------------------------------------------

  subroutine source_network_group_assign(self, data, offset)
    !! Assigns pointers in source network group object to elements in
    !! the data array, starting from the specified offset.

    class(source_network_group_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! source data array
    PetscInt, intent(in) :: offset  !! source array offset
    ! Locals:
    PetscInt :: group_offset

    call self%source_network_node_type%assign(data, offset)

    group_offset = offset + num_source_network_node_variables + &
         num_separator_variables

    self%group_index => data(group_offset)

  end subroutine source_network_group_assign

!------------------------------------------------------------------------

  subroutine source_network_group_init_data(self, group_index, separator_pressure, &
       thermo)
    !! Initialised source network group variables accessed via
    !! pointers to the source group vector. The group assign() method
    !! must be called first.

    class(source_network_group_type), intent(in out) :: self
    PetscInt, intent(in) :: group_index !! Index of group in input
    PetscReal, intent(in) :: separator_pressure(:) !! Separator pressures ([-1] for no separator)
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics

    self%group_index = dble(group_index)
    call self%separator%init(separator_pressure, thermo)

  end subroutine source_network_group_init_data

!------------------------------------------------------------------------

  subroutine source_network_group_set_rate(self, rate)
    !! Sets source network group flow rate to specified value, on the
    !! root process of the group.

    class(source_network_group_type), intent(in out) :: self
    PetscReal, intent(in) :: rate !! Flow rate

    if (self%rank == 0) then
       self%rate = rate
    end if
    call self%get_separated_flows()

  end subroutine source_network_group_set_rate

!------------------------------------------------------------------------

  subroutine source_network_group_scale_rate(self, scale)
    !! Scales source network group flow rate by specified scale
    !! factor, by scaling flows in group nodes.

    class(source_network_group_type), intent(in out) :: self
    PetscReal, intent(in) :: scale !! Flow rate scale factor

    call self%nodes%traverse(group_scale_iterator)
    call self%sum()

  contains

    subroutine group_scale_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (s => node%data)
      class is (source_network_node_type)
         call s%scale_rate(scale)
      end select

    end subroutine group_scale_iterator

  end subroutine source_network_group_scale_rate

!------------------------------------------------------------------------

  subroutine source_network_group_sum(self)
    !! Computes total flow rate, enthalpy etc. in a source group. The
    !! results are stored only on the root rank of the group
    !! communicator.

    use dm_utils_module, only: global_section_offset

    class(source_network_group_type), intent(in out) :: self
    ! Locals:
    PetscReal :: local_q, local_qh
    PetscReal :: total_q, total_qh
    PetscErrorCode :: ierr
    PetscReal, parameter :: rate_tol = 1.e-9_dp

    local_q = 0._dp
    local_qh = 0._dp

    call self%nodes%traverse(group_sum_iterator)

    call MPI_reduce(local_q, total_q, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
    call MPI_reduce(local_qh, total_qh, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)

    if (self%rank == 0) then
       if (abs(total_q) > rate_tol) then
          self%enthalpy = total_qh / total_q
       else
          self%enthalpy = 0._dp
       end if
    end if
    call self%set_rate(total_q)

  contains

    subroutine group_sum_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (s => node%data)
      class is (source_network_node_type)
         call s%add_flows(local_q, local_qh)
      end select

    end subroutine group_sum_iterator

  end subroutine source_network_group_sum

!------------------------------------------------------------------------

  subroutine source_network_group_default_separated_flows(self)
      !! Gets default separated water and steam flows when the group
      !! does not have its own separator. These are calculated by
      !! summing the separated flows from the inputs.

    class(source_network_group_type), intent(in out) :: self

    ! Locals:
    PetscReal :: local_water_q, local_water_qh
    PetscReal :: total_water_q, total_water_qh
    PetscReal :: local_steam_q, local_steam_qh
    PetscReal :: total_steam_q, total_steam_qh
    PetscErrorCode :: ierr
    PetscReal, parameter :: rate_tol = 1.e-9_dp

    local_water_q = 0._dp
    local_water_qh = 0._dp
    local_steam_q = 0._dp
    local_steam_qh = 0._dp

    call self%nodes%traverse(group_sum_separated_iterator)

    call MPI_reduce(local_water_q, total_water_q, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
    call MPI_reduce(local_water_qh, total_water_qh, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
    call MPI_reduce(local_steam_q, total_steam_q, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
    call MPI_reduce(local_steam_qh, total_steam_qh, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)

    if (self%rank == 0) then
       self%water_rate = total_water_q
       self%steam_rate = total_steam_q
       if (abs(total_water_q) > rate_tol) then
          self%water_enthalpy = total_water_qh / total_water_q
       else
          self%water_enthalpy = 0._dp
       end if
       if (abs(total_steam_q) > rate_tol) then
          self%steam_enthalpy = total_steam_qh / total_steam_q
       else
          self%steam_enthalpy = 0._dp
       end if
    end if

  contains

    subroutine group_sum_separated_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (s => node%data)
      class is (source_network_node_type)
         call s%add_separated_flows(local_water_q, local_water_qh, &
              local_steam_q, local_steam_qh)
      end select

    end subroutine group_sum_separated_iterator

  end subroutine source_network_group_default_separated_flows

!------------------------------------------------------------------------

  subroutine source_network_group_get_separated_flows(self)
      !! Gets separated water and steam flows on the root process of
      !! the group.

    class(source_network_group_type), intent(in out) :: self
    ! Locals:
    PetscBool :: use_default
    PetscErrorCode :: ierr

    use_default = PETSC_FALSE

    if (self%rank == 0) then
       if (self%rate < 0._dp) then
          if (self%separator%on) then
             call self%separate()
          else
             use_default = PETSC_TRUE
          end if
       else
          call self%zero_separated_flows()
       end if
    end if

    call MPI_bcast(use_default, 1, MPI_LOGICAL, 0, self%comm, ierr)
    if (use_default) then
       call self%default_separated_flows()
    end if

  end subroutine source_network_group_get_separated_flows

!------------------------------------------------------------------------

  subroutine source_network_group_add_flows(self, mass_flow, energy_flow)
    !! Adds mass and energy flows from the group to the specified
    !! totals (only on the root process of the group).

    class(source_network_group_type), intent(in out) :: self
    PetscReal, intent(in out) :: mass_flow !! Total mass flow rate
    PetscReal, intent(in out) :: energy_flow !! Total energy flow rate

    if (self%rank == 0) then
       call self%source_network_node_type%add_flows(mass_flow, energy_flow)
    end if

  end subroutine source_network_group_add_flows

!------------------------------------------------------------------------

  subroutine source_network_group_add_separated_flows(self, water_mass_flow, &
       water_energy_flow, steam_mass_flow, steam_energy_flow)
    !! Adds mass and energy flows for separated water and steam from
    !! the group to the specified totals, on the group's root process
    !! only.

    class(source_network_group_type), intent(in out) :: self
    PetscReal, intent(in out) :: water_mass_flow !! Total water mass flow rate
    PetscReal, intent(in out) :: water_energy_flow !! Total water energy flow rate
    PetscReal, intent(in out) :: steam_mass_flow !! Total steam mass flow rate
    PetscReal, intent(in out) :: steam_energy_flow !! Total steam energy flow rate

    if (self%rank == 0) then
       call self%source_network_node_type%add_separated_flows(water_mass_flow, &
       water_energy_flow, steam_mass_flow, steam_energy_flow)
    end if

  end subroutine source_network_group_add_separated_flows

!------------------------------------------------------------------------

  subroutine source_network_group_destroy(self)
    !! Destroys a source network group.

    class(source_network_group_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call self%nodes%destroy()
    call MPI_comm_free(self%comm, ierr)
    call self%source_network_node_type%destroy()

  end subroutine source_network_group_destroy
    
!------------------------------------------------------------------------

end module source_network_group_module
