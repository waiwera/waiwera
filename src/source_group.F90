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

module source_group_module
  !! Module for source groups.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use source_network_module
  use source_module, only: source_type
  use separator_module, only: num_separator_variables, separator_variable_names
  use list_module, only: list_type, list_node_type

  implicit none
  private

  PetscInt, parameter, public :: num_source_group_variables = &
       num_source_network_node_variables + num_separator_variables + 1
  PetscInt, parameter, public :: max_source_group_variable_name_length = 24
  character(max_source_group_variable_name_length), parameter, public :: &
       source_group_variable_names(num_source_group_variables) = [ &
       source_network_variable_names, separator_variable_names, [ &
       "group_index         "]]
  PetscInt, parameter, public :: num_source_group_constant_integer_variables = 1
  character(max_source_group_variable_name_length), parameter, public :: &
       source_group_constant_integer_variables( &
       num_source_group_constant_integer_variables) = [ &
       "group_index       "]

  type, public, extends(source_network_node_type) :: source_group_type
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
     procedure, public :: init => source_group_init
     procedure, public :: init_comm => source_group_init_comm
     procedure, public :: assign => source_group_assign
     procedure, public :: init_data => source_group_init_data
     procedure, public :: sum => source_group_sum
     procedure, public :: destroy => source_group_destroy
  end type source_group_type

contains

!------------------------------------------------------------------------

  subroutine source_group_init(self, name)
    !! Initialises a source group. Only variables stored in the object
    !! itself are initialised, not those stored in the group data
    !! vector and accessed via pointers.

    class(source_group_type), intent(in out) :: self
    character(*), intent(in) :: name !! Group name

    self%name = name
    call self%nodes%init(owner = PETSC_FALSE)

  end subroutine source_group_init

!------------------------------------------------------------------------

  subroutine source_group_init_comm(self)
    !! Initialises MPI communicator for the group. It is assumed that
    !! the nodes list has already been populated.

    class(source_group_type), intent(in out) :: self
    ! Locals:
    PetscInt :: colour, group_rank
    PetscErrorCode :: ierr

    if (self%nodes%count > 0) then
       colour = 1
    else
       colour = MPI_UNDEFINED
    end if

    call MPI_comm_split(PETSC_COMM_WORLD, colour, 0, self%comm, ierr)

    if (self%comm /= MPI_COMM_NULL) then
       call MPI_COMM_RANK(self%comm, group_rank, ierr)
       self%is_root = (group_rank == 0)
    else
       self%is_root = PETSC_FALSE
    end if

  end subroutine source_group_init_comm

!------------------------------------------------------------------------

  subroutine source_group_assign(self, data, offset)
    !! Assigns pointers in source group object to elements in the data
    !! array, starting from the specified offset.

    class(source_group_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! source data array
    PetscInt, intent(in) :: offset  !! source array offset
    ! Locals:
    PetscInt :: group_offset

    call self%source_network_node_type%assign(data, offset)

    group_offset = offset + num_source_network_node_variables + &
         num_separator_variables

    self%group_index => data(group_offset)

  end subroutine source_group_assign

!------------------------------------------------------------------------

  subroutine source_group_init_data(self, group_index)
    !! Initialised source group variables accessed via pointers to the
    !! source group vector. The group assign() method must be called
    !! first.

    class(source_group_type), intent(in out) :: self
    PetscInt, intent(in) :: group_index !! Index of group in input

    self%group_index = dble(group_index)

  end subroutine source_group_init_data

!------------------------------------------------------------------------

  subroutine source_group_sum(self)
    !! Computes total flow rate, enthalpy etc. in a source group. The
    !! results are stored only on the root rank of the group
    !! communicator.

    class(source_group_type), intent(in out) :: self
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
       self%rate = total_q
       if (abs(total_q) > rate_tol) then
          self%enthalpy = total_qh / total_q
       else
          self%enthalpy = 0._dp
       end if
    end if

  contains

    subroutine group_sum_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (source => node%data)
      type is (source_type)
         local_q = local_q + source%rate
         local_qh = local_qh + source%rate * source%enthalpy
      end select

    end subroutine group_sum_iterator

  end subroutine source_group_sum

!------------------------------------------------------------------------

  subroutine source_group_destroy(self)
    !! Destroys a source group.

    class(source_group_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call self%nodes%destroy()
    call MPI_comm_free(self%comm, ierr)
    call self%source_network_node_type%destroy()

  end subroutine source_group_destroy
    
!------------------------------------------------------------------------

end module source_group_module
