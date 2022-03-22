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

    type, public, extends(multi_table_object_control_type) :: &
         limiter_table_source_network_control_type
     !! Limits flows (total, separated water or steam) through a source network node.
     private
     PetscInt, allocatable, public :: flow_type(:) !! Types of flow being limited - total, water or steam
   contains
     private
     procedure, public :: iterator => limiter_table_source_network_control_iterator
  end type limiter_table_source_network_control_type

contains

!------------------------------------------------------------------------
! Limiter table source network controls:
!------------------------------------------------------------------------

  subroutine limiter_table_source_network_control_iterator(self, node, &
       stopped)
    !! Update flow so limits are not exceeded.

    class(limiter_table_source_network_control_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node !! List node
    PetscBool, intent(out) :: stopped

    stopped = PETSC_FALSE
    select type(network_node => node%data)
    class is (source_network_node_type)
       associate(limit => self%value)
         call network_node%limit_rate(self%flow_type, limit)
       end associate
    end select

  end subroutine limiter_table_source_network_control_iterator

!------------------------------------------------------------------------  

end module source_network_control_module
