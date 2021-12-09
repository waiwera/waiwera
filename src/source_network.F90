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

  type, public :: source_network_node_type
     !! Type for node in source network, e.g. source or source group.
     private
     PetscReal, pointer, public :: rate !! Flow rate
     PetscReal, pointer, public :: enthalpy !! Enthalpy of fluid
     type(separator_type), public :: separator !! Separator
     PetscBool, public :: heat !! Whether flow through node is heat-only (not mass)
   contains
     private
     procedure, public :: destroy => source_network_node_destroy
  end type source_network_node_type

contains

!------------------------------------------------------------------------

  subroutine source_network_node_destroy(self)
    !! Destroys source network node.

    class(source_network_node_type), intent(in out) :: self

    self%rate => null()
    self%enthalpy => null()

  end subroutine source_network_node_destroy

!------------------------------------------------------------------------

end module source_network_module
