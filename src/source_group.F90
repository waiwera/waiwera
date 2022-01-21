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
     !! Type for group of source network nodes, e.g. multi-feed well
     !! or group of makeup wells.
     private
     type(list_type), public :: nodes !! List of nodes in group
     MPI_Comm :: comm !! MPI communicator for group
   contains
     private
     procedure, public :: init => source_group_init
     procedure, public :: init_comm => source_group_init_comm
     procedure, public :: assign => source_group_assign
     procedure, public :: destroy => source_group_destroy
  end type source_group_type

contains

!------------------------------------------------------------------------

  subroutine source_group_init(self, name)
    !! Initialises a source group.

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
    PetscInt :: colour
    PetscErrorCode :: ierr

    if (self%nodes%count > 0) then
       colour = 1
    else
       colour = MPI_UNDEFINED
    end if

    call MPI_comm_split(PETSC_COMM_WORLD, colour, 0, self%comm, ierr)

  end subroutine source_group_init_comm

!------------------------------------------------------------------------

  subroutine source_group_assign(self, data, offset)
    !! Assigns pointers in source group object to elements in the data
    !! array, starting from the specified offset.

    class(source_group_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! source data array
    PetscInt, intent(in) :: offset  !! source array offset

    call self%source_network_node_type%assign(data, offset)

  end subroutine source_group_assign

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
