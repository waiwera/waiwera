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

module source_network_module
  !! Module for source network nodes.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use separator_module

  implicit none
  private

  PetscInt, parameter, public :: max_source_network_node_name_length = 32
  PetscInt, parameter, public :: num_source_network_node_variables = 2
  PetscInt, parameter, public :: max_source_network_variable_name_length = 24
  character(max_source_network_variable_name_length), parameter, public :: &
       source_network_variable_names(num_source_network_node_variables) = [ &
       "rate                ", "enthalpy            "]

  type, public :: source_network_node_type
     !! Type for node in source network, e.g. source or source group.
     private
     character(max_source_network_node_name_length), public :: name !! Name of source network node
     PetscReal, pointer, public :: rate !! Flow rate
     PetscReal, pointer, public :: enthalpy !! Enthalpy of fluid
     type(separator_type), public :: separator !! Separator
     PetscBool, public :: heat !! Whether flow through node is heat-only (not mass)
   contains
     private
     procedure, public :: assign => source_network_node_assign
     procedure, public :: set_rate => source_network_node_set_rate
     procedure, public :: destroy => source_network_node_destroy
  end type source_network_node_type

contains

!------------------------------------------------------------------------

  subroutine source_network_node_assign(self, data, offset)
    !! Assigns pointers in network node object to elements in the data
    !! array, starting from the specified offset.

    class(source_network_node_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! source data array
    PetscInt, intent(in) :: offset  !! source array offset

    self%rate => data(offset)
    self%enthalpy => data(offset + 1)

    call self%separator%assign(data, offset + 2)

  end subroutine source_network_node_assign

!------------------------------------------------------------------------

  subroutine source_network_node_set_rate(self, rate)
    !! Sets network node flow rate to specified value.

    class(source_network_node_type), intent(in out) :: self
    PetscReal, intent(in) :: rate !! Flow rate

    self%rate = rate

    ! If have separator and producing mass:
    if ((self%separator%on) .and. (self%rate <= 0._dp) &
         .and. (.not. self%heat)) then
       call self%separator%separate(self%rate, self%enthalpy)
    else
       call self%separator%zero()
    end if

  end subroutine source_network_node_set_rate

!------------------------------------------------------------------------

  subroutine source_network_node_destroy(self)
    !! Destroys source network node.

    class(source_network_node_type), intent(in out) :: self

    self%rate => null()
    self%enthalpy => null()

  end subroutine source_network_node_destroy

!------------------------------------------------------------------------

end module source_network_module
