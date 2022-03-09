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

module source_network_control_module
  !! Module for source network controls- for controlling source network node parameters over time.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use control_module
  use source_network_module
  use list_module

  PetscInt, parameter, public :: max_limiter_type_length = 5
  character(max_limiter_type_length), parameter, public :: &
       default_source_control_limiter_type_str = "total"
    PetscReal, parameter, public :: default_source_control_limiter_limit = 1._dp

  type, public, extends(table_object_control_type) :: limiter_table_source_network_control_type
     !! Limits total flow through a source network node. Derived types limit
     !! separated water or steam flow rather than total flow.
   contains
     private
     procedure :: get_rate => limiter_table_source_network_control_get_rate
     procedure :: get_scale => limiter_table_source_network_control_get_scale
     procedure, public :: iterator => limiter_table_source_network_control_iterator
  end type limiter_table_source_network_control_type

  type, public, extends(limiter_table_source_network_control_type) :: &
       water_limiter_table_source_network_control_type
     !! Limits separated water flow (e.g. from a separator) through a source network node.
   contains
     private
     procedure :: get_rate => water_limiter_table_source_network_control_get_rate
  end type water_limiter_table_source_network_control_type

  type, public, extends(limiter_table_source_network_control_type) :: &
       steam_limiter_table_source_network_control_type
     !! Limits separated steam flow (e.g. from a separator) through a source network node.
   contains
     private
     procedure :: get_rate => steam_limiter_table_source_network_control_get_rate
  end type steam_limiter_table_source_network_control_type

contains

!------------------------------------------------------------------------
! Limiter table source network controls:
!------------------------------------------------------------------------

  PetscReal function limiter_table_source_network_control_get_rate(self, &
       network_node) result(rate)
    !! Get total flow rate from source network node.
    class(limiter_table_source_network_control_type), intent(in out) :: self
    class(source_network_node_type), intent(in) :: network_node

    rate = network_node%rate

  end function limiter_table_source_network_control_get_rate

!------------------------------------------------------------------------

  PetscReal function water_limiter_table_source_network_control_get_rate(self, &
       network_node) result(rate)
    !! Get separated water rate from source network node.
    class(water_limiter_table_source_network_control_type), intent(in out) :: self
    class(source_network_node_type), intent(in) :: network_node

    rate = network_node%water_rate

  end function water_limiter_table_source_network_control_get_rate

!------------------------------------------------------------------------

  PetscReal function steam_limiter_table_source_network_control_get_rate(self, &
       network_node) result(rate)
    !! Get separated steam rate from source network node.
    class(steam_limiter_table_source_network_control_type), intent(in out) :: self
    class(source_network_node_type), intent(in) :: network_node

    rate = network_node%steam_rate

  end function steam_limiter_table_source_network_control_get_rate

!------------------------------------------------------------------------

  PetscReal function limiter_table_source_network_control_get_scale(self, &
       rate, limit) result(scale)

    class(limiter_table_source_network_control_type), intent(in out) :: self
    PetscReal, intent(in) :: rate, limit
    ! Locals:
    PetscReal :: abs_rate
    PetscReal, parameter :: small = 1.e-6_dp

    abs_rate = abs(rate)
    if ((abs_rate > limit) .and. (abs_rate > small)) then
       scale = limit / abs_rate
    else
       scale = 1._dp
    end if

  end function limiter_table_source_network_control_get_scale

!------------------------------------------------------------------------

  subroutine limiter_table_source_network_control_iterator(self, node, &
       stopped)
    !! Update flow so limit is not exceeded.

    class(limiter_table_source_network_control_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node !! List node
    PetscBool, intent(out) :: stopped
    ! Locals:
    PetscReal :: rate, scale

    stopped = PETSC_FALSE
    select type(network_node => node%data)
    class is (source_network_node_type)
       associate(limit => self%value(1))
         rate = self%get_rate(network_node)
         scale = self%get_scale(rate, limit)
         call network_node%scale_rate(scale)
       end associate
    end select

  end subroutine limiter_table_source_network_control_iterator

!------------------------------------------------------------------------  

end module source_network_control_module
