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
     PetscBool, public :: is_root !! Whether group is on root rank of its communicator
     PetscInt, public :: local_group_index !! Index of group in local part of source group vector (-1 if not a root group)
     PetscReal, pointer, public :: group_index !! Index of source group in input
   contains
     private
     procedure, public :: init => source_network_group_init
     procedure, public :: init_comm => source_network_group_init_comm
     procedure, public :: assign => source_network_group_assign
     procedure, public :: init_data => source_network_group_init_data
     procedure, public :: sum => source_network_group_sum
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
       call MPI_COMM_RANK(self%comm, group_rank, ierr)
       self%is_root = (group_rank == 0)
    else
       self%is_root = PETSC_FALSE
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
         if (n%is_root) then
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

  subroutine source_network_group_sum(self, source_data, source_section, &
       source_range_start, source_network_group_data, source_network_group_section, &
       source_network_group_range_start)
    !! Computes total flow rate, enthalpy etc. in a source group. The
    !! results are stored only on the root rank of the group
    !! communicator.

    use dm_utils_module, only: global_section_offset

    class(source_network_group_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in out) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: source_network_group_data(:)
    PetscSection, intent(in out) :: source_network_group_section
    PetscInt, intent(in) :: source_network_group_range_start
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

    if (self%is_root) then
       if (abs(total_q) > rate_tol) then
          self%enthalpy = total_qh / total_q
       else
          self%enthalpy = 0._dp
       end if
       call self%set_rate(total_q)
    end if

  contains

    subroutine group_sum_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: offset

      stopped = PETSC_FALSE
      select type (s => node%data)
      type is (source_type)
         offset = global_section_offset(source_section, &
              s%local_source_index, source_range_start)
         call s%assign(source_data, offset)
         local_q = local_q + s%rate
         local_qh = local_qh + s%rate * s%enthalpy
      type is (source_network_group_type)
         if (s%is_root) then
            offset = global_section_offset(source_network_group_section, &
                 s%local_group_index, source_network_group_range_start)
            call s%assign(source_network_group_data, offset)
            local_q = local_q + s%rate
            local_qh = local_qh + s%rate * s%enthalpy
         end if
      end select

    end subroutine group_sum_iterator

  end subroutine source_network_group_sum

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
