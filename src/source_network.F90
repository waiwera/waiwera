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

  implicit none
  private

  type, public :: source_network_type
     !! Type for source network.
     private
     PetscInt, public :: source_range_start !! Range start for source vector
     PetscInt, public :: group_range_start !! Range start for source group vector
     Vec, public :: source !! Vector for source/sink data
     Vec, public :: group !! Vector for source network group data
     type(list_type), public :: sources !! List of source objects
     PetscInt, public :: num_sources !! Total number of source/sink terms on all processes
     PetscInt, public :: num_groups !! Total number of source network groups on all processes
     IS, public :: source_index !! Index set defining natural to global source ordering
     IS, public :: group_index !! Index set defining natural to global source network group ordering
     type(list_type), public :: source_controls !! Controls on sources/sinks
     type(list_type), public :: groups !! Source network groups
     type(list_type), public :: network_controls !! Controls on source network nodes
     type(list_type), public :: separated_sources !! Sources with separators
   contains
     private
     procedure, public :: destroy => source_network_destroy
  end type source_network_type

contains

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

    call ISDestroy(self%source_index, ierr); CHKERRQ(ierr)
    call ISDestroy(self%group_index, ierr); CHKERRQ(ierr)

    call self%separated_sources%destroy()
    call self%source_controls%destroy(object_control_list_node_data_destroy, &
         reverse = PETSC_TRUE)
    call self%groups%destroy( &
         source_network_node_list_node_data_destroy, reverse = PETSC_TRUE)
    call self%network_controls%destroy( &
         object_control_list_node_data_destroy, reverse = PETSC_TRUE)
    call self%sources%destroy(source_network_node_list_node_data_destroy)

  contains

!........................................................................
    
    
  end subroutine source_network_destroy
  
!------------------------------------------------------------------------

end module source_network_module
  
