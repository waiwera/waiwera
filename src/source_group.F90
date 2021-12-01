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
  use list_module, only: list_type

  implicit none
  private

  type, public, extends(source_network_node_type) :: source_group_type
     !! Type for group of sources, e.g. multi-feed well or group of
     !! makeup wells.
     private
     type(list_type), public :: sources !! List of sources in group
   contains
     private
     procedure, public :: init => source_group_init
     procedure, public :: destroy => source_group_destroy
  end type source_group_type

contains

!------------------------------------------------------------------------

  subroutine source_group_init(self)
    !! Initialises a source group.

    class(source_group_type), intent(in out) :: self
    
  end subroutine source_group_init

!------------------------------------------------------------------------

  subroutine source_group_destroy(self)
    !! Destroys a source group.

    class(source_group_type), intent(in out) :: self

    call self%sources%destroy()

  end subroutine source_group_destroy
    
!------------------------------------------------------------------------

end module source_group_module
