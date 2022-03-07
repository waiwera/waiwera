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

module control_module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use list_module
  use interpolation_module

  implicit none
  private

  type, public, abstract :: object_control_type
     !! Type for controlling how object parameters vary with time.
     private
     type(list_type), public :: objects
   contains
     private
     procedure(object_control_iterator_procedure), deferred, public :: iterator
     procedure(object_control_destroy_procedure), deferred, public :: destroy
  end type object_control_type

  type, public, extends(object_control_type) :: table_object_control_type
     !! Controls object parameters using an interpolation table of
     !! time-dependent values.
     private
     type(interpolation_table_type), public :: table !! Table of values vs. time
     PetscReal, allocatable, public :: value(:) !! Interpolated value
   contains
     private
     procedure, public :: init => table_object_control_init
     procedure, public :: update => table_object_control_update
     procedure, public :: iterator => table_object_control_iterator
     procedure, public :: destroy => table_object_control_destroy
  end type table_object_control_type

  abstract interface

     subroutine object_control_iterator_procedure(self, node, stopped)
       !! Wrapper for objects update iterator
       import :: object_control_type, list_node_type
       class(object_control_type), intent(in out) :: self
       type(list_node_type), pointer, intent(in out) :: node
       PetscBool, intent(out) :: stopped
     end subroutine object_control_iterator_procedure

     subroutine object_control_destroy_procedure(self)
       !! Object control destroy procedure
       import :: object_control_type
       class(object_control_type), intent(in out) :: self
     end subroutine object_control_destroy_procedure

  end interface

contains

!------------------------------------------------------------------------
! Objcet table control routines
!------------------------------------------------------------------------

  subroutine table_object_control_init(self, objects, data, interpolation_type, &
       averaging_type)
    !! Initialises table_object_control object.

    class(table_object_control_type), intent(in out) :: self
    type(list_type), intent(in) :: objects
    PetscReal, intent(in) :: data(:,:)
    PetscInt, intent(in) :: interpolation_type
    PetscInt, intent(in) :: averaging_type

    self%objects = objects
    call self%table%init(data, interpolation_type, averaging_type)
    allocate(self%value(self%table%dim))

  end subroutine table_object_control_init

!------------------------------------------------------------------------

  subroutine table_object_control_iterator(self, node, stopped)
    !! Update iterator for table object control list. Derived types
    !! will use this routine to update their objects using self%value.

    class(table_object_control_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node
    PetscBool, intent(out) :: stopped

    stopped = PETSC_FALSE
    
  end subroutine table_object_control_iterator

!------------------------------------------------------------------------

  subroutine table_object_control_update(self, interval)
    !! Updates table object control object properties using average of
    !! table over the specified time interval.

    class(table_object_control_type), intent(in out) :: self
    PetscReal, intent(in) :: interval(2)

    self%value = self%table%average(interval)
    call self%objects%traverse(iterator)

  contains

    subroutine iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      call self%iterator(node, stopped)

    end subroutine iterator
    
  end subroutine table_object_control_update

!------------------------------------------------------------------------

  subroutine table_object_control_destroy(self)
    !! Destroys a table object control.

    class(table_object_control_type), intent(in out) :: self

    call self%objects%destroy()
    call self%table%destroy()
    deallocate(self%value)

  end subroutine table_object_control_destroy

!------------------------------------------------------------------------

end module control_module
