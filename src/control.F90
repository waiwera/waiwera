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
     type(list_type), public :: objects !! List of objects to control
   contains
     private
     procedure(object_control_destroy_procedure), deferred, public :: destroy
  end type object_control_type

  type, public, extends(object_control_type) :: integer_object_control_type
     !! Controls object parameters based on the value of a single
     !! integer control parameter.
     private
     PetscInt, public :: value !! Integer control parameter
   contains
     procedure, public :: init => integer_object_control_init
     procedure, public :: update => integer_object_control_update
     procedure, public :: iterator => integer_object_control_iterator
     procedure, public :: destroy => integer_object_control_destroy
  end type integer_object_control_type

  type, public, abstract, extends(object_control_type) :: &
       interval_update_object_control_type
     !! Controls with an update method depending only on a time
     !! interval.
   contains
     procedure(interval_update_object_control_update_procedure), deferred, &
          public :: update
  end type

  type, public, extends(interval_update_object_control_type) :: table_object_control_type
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

  type, public, extends(interval_update_object_control_type) :: multi_table_object_control_type
     !! Controls object parameters using an array of interpolation tables of
     !! time-dependent values.
     private
     PetscInt, public :: num_tables !! Number of interpolation tables
     type(interpolation_table_type), allocatable, public :: table(:) !! Tables of values vs. time
     PetscReal, allocatable, public :: value(:) !! Interpolated scalar values for each table
   contains
     private
     procedure, public :: init => multi_table_object_control_init
     procedure, public :: init_table => multi_table_object_control_init_table
     procedure, public :: update => multi_table_object_control_update
     procedure, public :: iterator => multi_table_object_control_iterator
     procedure, public :: destroy => multi_table_object_control_destroy
  end type multi_table_object_control_type

  type, public, abstract :: vector_control_type
     !! Type for controlling how selected vector values vary with
     !! time.
     private
     PetscInt, allocatable, public :: indices(:) !! Array indices (block indices for block vector)
   contains
     private
     procedure(vector_control_destroy_procedure), deferred, public :: destroy
  end type vector_control_type

  type, public, extends(vector_control_type) :: table_vector_control_type
     !! Controls vector values using an interpolation table of
     !! time-dependent values.
     private
     type(interpolation_table_type), public :: table !! Table of values vs. time
   contains
     private
     procedure, public :: init => table_vector_control_init
     procedure, public :: update => table_vector_control_update
     procedure, public :: destroy => table_vector_control_destroy
  end type table_vector_control_type

  abstract interface

     subroutine object_control_destroy_procedure(self)
       !! Object control destroy procedure
       import :: object_control_type
       class(object_control_type), intent(in out) :: self
     end subroutine object_control_destroy_procedure

     subroutine interval_update_object_control_update_procedure(self, interval)
       !! Update procedure based on a time interval
       import :: interval_update_object_control_type
       class(interval_update_object_control_type), intent(in out) :: self
       PetscReal, intent(in) :: interval(2)
     end subroutine interval_update_object_control_update_procedure

     subroutine vector_control_destroy_procedure(self)
       !! Vector control destroy procedure
       import :: vector_control_type
       class(vector_control_type), intent(in out) :: self
     end subroutine vector_control_destroy_procedure

  end interface

  public :: object_control_list_node_data_destroy

contains

!------------------------------------------------------------------------

  subroutine object_control_list_node_data_destroy(node)
    ! Destroys object control in a list node.

    use list_module, only: list_node_type

    type(list_node_type), pointer, intent(in out) :: node

    select type (control => node%data)
    class is (object_control_type)
       call control%destroy()
    end select

  end subroutine object_control_list_node_data_destroy

!------------------------------------------------------------------------
! Integer object control routines
!------------------------------------------------------------------------

  subroutine integer_object_control_init(self, objects, value)
    !! Initialises integer_object_control object.

    class(integer_object_control_type), intent(in out) :: self
    type(list_type), intent(in) :: objects
    PetscInt, intent(in) :: value

    self%objects = objects
    self%value = value

  end subroutine integer_object_control_init

!------------------------------------------------------------------------

  subroutine integer_object_control_iterator(self, node, stopped)
    !! Update iterator for integer object control list. Derived types
    !! will use this routine to update their objects using self%value.

    class(integer_object_control_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node
    PetscBool, intent(out) :: stopped

    stopped = PETSC_FALSE

  end subroutine integer_object_control_iterator

!------------------------------------------------------------------------

  subroutine integer_object_control_update(self)
    !! Updates integer object control object properties using
    !! self%value.

    class(integer_object_control_type), intent(in out) :: self

    call self%objects%traverse(iterator)

  contains

    subroutine iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      call self%iterator(node, stopped)

    end subroutine iterator

  end subroutine integer_object_control_update

!------------------------------------------------------------------------

  subroutine integer_object_control_destroy(self)
    !! Destroys an integer object control.

    class(integer_object_control_type), intent(in out) :: self

    call self%objects%destroy()

  end subroutine integer_object_control_destroy

!------------------------------------------------------------------------
! Table object control routines
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
    if (allocated(self%value)) deallocate(self%value)
    call self%table%destroy()

  end subroutine table_object_control_destroy

!------------------------------------------------------------------------
! Mutil-table object control routines
!------------------------------------------------------------------------

  subroutine multi_table_object_control_init(self, objects, num_tables)

    !! Initialises table_object_control object.

    class(multi_table_object_control_type), intent(in out) :: self
    type(list_type), intent(in) :: objects !! Objects to control
    PetscInt, intent(in) :: num_tables !! Number of interpolation tables

    self%objects = objects
    self%num_tables = num_tables
    allocate(self%table(num_tables), self%value(num_tables))

  end subroutine multi_table_object_control_init

!------------------------------------------------------------------------

  subroutine multi_table_object_control_init_table(self, index, data, &
       interpolation_type, averaging_type)
    !! Initialises interpolation table with specified index.

    class(multi_table_object_control_type), intent(in out) :: self
    PetscInt, intent(in) :: index !! Index of table
    PetscReal, intent(in) :: data(:,:) !! Interpolation data for table
    PetscInt, intent(in) :: interpolation_type !! Table interpolation type
    PetscInt, intent(in) :: averaging_type !! Table averaging type

    call self%table(index)%init(data, interpolation_type, averaging_type)

  end subroutine multi_table_object_control_init_table

!------------------------------------------------------------------------

  subroutine multi_table_object_control_iterator(self, node, stopped)
    !! Update iterator for multi-table object control list. Derived
    !! types will use this routine to update their objects using
    !! self%value.

    class(multi_table_object_control_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node
    PetscBool, intent(out) :: stopped

    stopped = PETSC_FALSE

  end subroutine multi_table_object_control_iterator

!------------------------------------------------------------------------

  subroutine multi_table_object_control_update(self, interval)
    !! Updates multi-table object control object properties using
    !! average of tables over the specified time interval.

    class(multi_table_object_control_type), intent(in out) :: self
    PetscReal, intent(in) :: interval(2)
    ! Locals:
    PetscInt :: i

    do i = 1, self%num_tables
       self%value(i) = self%table(i)%average(interval, 1)
    end do
    call self%objects%traverse(iterator)

  contains

    subroutine iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      call self%iterator(node, stopped)

    end subroutine iterator

  end subroutine multi_table_object_control_update

!------------------------------------------------------------------------

  subroutine multi_table_object_control_destroy(self)
    !! Destroys a multi-table object control.

    class(multi_table_object_control_type), intent(in out) :: self
    ! Locals:
    PetscInt :: i

    call self%objects%destroy()
    if (allocated(self%value)) deallocate(self%value)
    if (allocated(self%table)) then
       do i = 1, self%num_tables
          call self%table(i)%destroy()
       end do
       deallocate(self%table)
    end if

  end subroutine multi_table_object_control_destroy

!------------------------------------------------------------------------
! Table vector control routines
!------------------------------------------------------------------------

  subroutine table_vector_control_init(self, data, indices, &
       interpolation_type)
    !! Initialises table_object_control object.

    class(table_vector_control_type), intent(in out) :: self
    PetscReal, intent(in) :: data(:,:) !! Data for interpolation table
    PetscInt, intent(in) :: indices(:) !! Vector indices
    PetscInt, intent(in) :: interpolation_type !! Interpolation type for data

    call self%table%init(data, interpolation_type)
    self%indices = indices

  end subroutine table_vector_control_init

!------------------------------------------------------------------------

  subroutine table_vector_control_update(self, time, &
       vector_array, section, range_start)
    !! Updates vector values by interpolating table at the specified
    !! time. Derived types override this routine to update the vector
    !! values.

    class(table_vector_control_type), intent(in out) :: self
    PetscReal, intent(in) :: time !! Time of update
    PetscReal, pointer, contiguous, intent(in) :: vector_array(:) !! Array on vector
    PetscSection, intent(in) :: section !! Global section for vector
    PetscInt, intent(in) :: range_start !! Range start for vector

  end subroutine table_vector_control_update

!------------------------------------------------------------------------

  subroutine table_vector_control_destroy(self)
    !! Destroys a table vector control.

    class(table_vector_control_type), intent(in out) :: self

    if (allocated(self%indices)) deallocate(self%indices)
    call self%table%destroy()

  end subroutine table_vector_control_destroy

!------------------------------------------------------------------------

end module control_module
