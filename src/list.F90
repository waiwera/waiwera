!   Copyright 2016 University of Auckland.

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

module list_module
  !! Module for linked lists with arbitrary data and optional character tags.

#include <petsc/finclude/petscsys.h>

  use petscsys

  implicit none
  private

  type, public :: list_node_type
     !! List node type.
     private
     character(:), allocatable, public :: tag !! Optional character tag
     class(*), pointer, public :: data => null() !! Node data
     type(list_node_type), pointer, public :: next => null() !! Pointer to next node
     type(list_node_type), pointer, public :: previous => null() !! Pointer to previous node
  end type list_node_type

  type, public :: list_type
     !! Linked list type.
     private
     PetscBool, public :: owner !! Whether the list 'owns' the data
     PetscInt, public :: count  !! Number of nodes in the list
     type(list_node_type), pointer, public :: head !! Node at start of list
     type(list_node_type), pointer, public :: tail !! Node at end of list
   contains
     private
     procedure :: list_init
     procedure :: list_init_default
     generic, public :: init => list_init, list_init_default
     procedure, public :: append => list_append
     procedure, public :: prepend => list_prepend
     procedure, public :: add => list_add
     procedure :: list_delete
     procedure :: list_delete_tag
     procedure :: list_delete_destroy
     procedure :: list_delete_destroy_tag
     generic, public :: delete => list_delete, list_delete_tag, &
          list_delete_destroy, list_delete_destroy_tag
     procedure :: list_traverse
     procedure :: list_traverse_default
     generic, public :: traverse => list_traverse, list_traverse_default
     procedure :: list_find
     procedure :: list_find_default
     generic, public :: find => list_find, list_find_default
     procedure, public :: get => list_get
     procedure :: destroy_default_forward => list_destroy_default_forward
     procedure :: destroy_default_reverse => list_destroy_default_reverse
     procedure :: destroy_default_direction => list_destroy_default_direction
     procedure :: destroy_proc_forward => list_destroy_proc_forward
     procedure :: destroy_proc_reverse => list_destroy_proc_reverse
     procedure :: destroy_proc_direction => list_destroy_proc_direction
     generic, public :: destroy => destroy_default_forward, &
          destroy_default_direction, destroy_proc_forward, &
          destroy_proc_direction
     procedure, public :: tags => list_tags
     procedure, public :: copy => list_copy
  end type list_type

  abstract interface

     subroutine list_iterator(node, stopped)
       import :: list_node_type
       type(list_node_type), pointer, intent(in out)  :: node
       PetscBool, intent(out) :: stopped
     end subroutine list_iterator

     subroutine list_node_data_destroy_procedure(node)
       import :: list_node_type
       type(list_node_type), pointer, intent(in out)  :: node
     end subroutine list_node_data_destroy_procedure

  end interface

  public :: list_node_data_destroy_procedure

contains

!------------------------------------------------------------------------

  subroutine list_init(self, owner)
    !! Initialises list, specifying whether the list 'owns' the data,
    !! i.e. whether the data should be destroyed when the list is destroyed.
    
    class(list_type), intent(in out) :: self
    PetscBool, intent(in) :: owner

    self%owner = owner
    self%count = 0
    self%head => null()
    self%tail => null()

  end subroutine list_init

  subroutine list_init_default(self)
    !! Initialises list with default owner = false.

    class(list_type), intent(in out) :: self

    call self%init(PETSC_FALSE)

  end subroutine list_init_default

!------------------------------------------------------------------------

  subroutine list_append(self, data, tag)
    !! Adds node with specified data to end of list, with an optional
    !! tag.

    class(list_type), intent(in out) :: self
    class(*), target, intent(in), optional :: data
    character(len = *), intent(in), optional :: tag
    ! Locals:
    type(list_node_type), pointer :: old_tail

    if (associated(self%tail)) then
       old_tail => self%tail
       allocate(self%tail%next)
       self%tail => self%tail%next
       self%tail%previous => old_tail
    else ! empty:
       allocate(self%tail)
       self%head => self%tail
    end if

    self%count = self%count + 1

    if (present(data)) then
       self%tail%data => data
    end if

    if (present(tag)) then
       if (tag /= "") then
          self%tail%tag = tag
       end if
    end if

  end subroutine list_append

!------------------------------------------------------------------------

  subroutine list_prepend(self, data, tag)
    !! Adds node with specified data to start of list, with an
    !! optional tag.

    class(list_type), intent(in out) :: self
    class(*), target, intent(in), optional :: data
    character(len = *), intent(in), optional :: tag
    ! Locals:
    type(list_node_type), pointer :: old_head

    if (associated(self%head)) then
       old_head => self%head
       allocate(self%head%previous)
       self%head => self%head%previous
       self%head%next => old_head
    else ! empty:
       allocate(self%head)
       self%tail => self%head
    end if

    self%count = self%count + 1

    if (present(data)) then
       self%head%data => data
    end if

    if (present(tag)) then
       if (tag /= "") then
          self%head%tag  = tag
       end if
    end if

  end subroutine list_prepend

!------------------------------------------------------------------------

  subroutine list_add(self, other)
    !! Appends nodes with data (and any tags) from the other list to self.

    class(list_type), intent(in out) :: self
    type(list_type), intent(in out) :: other

    call other%traverse(append_to_self_iterator)

  contains

    subroutine append_to_self_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      if (allocated(node%tag)) then
         call self%append(node%data, node%tag)
      else
         call self%append(node%data)
      end if

    end subroutine append_to_self_iterator

  end subroutine list_add

!------------------------------------------------------------------------

  subroutine list_delete(self, node)
    !! Deletes specified list node. If self%owner is true,
    !! this will also deallocate the data in the node.

    class(list_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node
    ! Locals:
    type(list_node_type), pointer :: previous, next

    if (associated(node)) then

       if (self%owner) then
          deallocate(node%data)
       end if

       previous => node%previous
       next => node%next

       node%data => null()
       node%next => null()
       node%previous => null()
       if (allocated(node%tag)) deallocate(node%tag)
       deallocate(node)
       node => null()

       if (associated(previous)) then
          previous%next => next
       else
          self%head => next
       end if

       if (associated(next)) then
          next%previous => previous
       else
          self%tail => previous
       end if

       self%count = self%count - 1

    end if

  end subroutine list_delete

  recursive subroutine list_delete_destroy(self, node, node_data_destroy_procedure)
    !! Delete specified list node and calls the given procedure to
    !! destroy its data.

    class(list_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node
    procedure(list_node_data_destroy_procedure) :: node_data_destroy_procedure

    if (self%owner) then
       call node_data_destroy_procedure(node)
    end if
    call self%delete(node)

  end subroutine list_delete_destroy

  subroutine list_delete_tag(self, tag)
    !! Deletes first list node with specified tag, if one exists.

    class(list_type), intent(in out) :: self
    character(len = *), intent(in) :: tag
    ! Locals:
    type(list_node_type), pointer :: node

    node => self%find(tag)
    call self%delete(node)
    
  end subroutine list_delete_tag

  subroutine list_delete_destroy_tag(self, tag, node_data_destroy_procedure)

    class(list_type), intent(in out) :: self
    character(len = *), intent(in) :: tag
    procedure(list_node_data_destroy_procedure) :: node_data_destroy_procedure
    ! Locals:
    type(list_node_type), pointer :: node

    node => self%find(tag)
    call self%delete(node, node_data_destroy_procedure)

  end subroutine list_delete_destroy_tag

!------------------------------------------------------------------------

  recursive subroutine list_traverse(self, iterator, backwards)
    !! Traverses the list, applying the iterator function to each
    !! node. If backwards is true, the list is traversed from tail to
    !! head, otherwise from head to tail.

    class(list_type), intent(in out) :: self
    procedure(list_iterator)  :: iterator
    PetscBool, intent(in) :: backwards
    ! Locals:
    type(list_node_type), pointer :: node
    PetscBool :: stopped

    stopped = PETSC_FALSE

    if (backwards) then
       node => self%tail
    else
       node => self%head
    end if
    
    if (associated(node)) then
       do
          call iterator(node, stopped)
          if (backwards) then
             node => node%previous
          else
             node => node%next
          end if
          if (stopped .or. .not. (associated(node))) then
             exit
          end if
       end do
    end if

  end subroutine list_traverse

  recursive subroutine list_traverse_default(self, iterator)
    !! Traverses the list in the default forward direction.

    class(list_type), intent(in out) :: self
    procedure(list_iterator)  :: iterator

    call self%traverse(iterator, PETSC_FALSE)

  end subroutine list_traverse_default

!------------------------------------------------------------------------

  function list_find(self, tag, backwards) result(found_node)
    !! Returns pointer to first list node with specified tag, if one
    !! exists, or a null pointer otherwise. If backwards is true, the
    !! list is searched from tail to head, otherwise from head to
    !! tail.

    class(list_type), intent(in out) :: self
    character(len = *), intent(in) :: tag
    PetscBool, intent(in) :: backwards
    type(list_node_type), pointer :: found_node

    found_node => null()
    call self%traverse(tag_find_iterator, backwards)

  contains

    subroutine tag_find_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      if (allocated(node%tag) .and. (node%tag == tag)) then
         stopped = PETSC_TRUE
         found_node => node
      else
         stopped = PETSC_FALSE
      end if
    end subroutine tag_find_iterator

  end function list_find

  function list_find_default(self, tag) result(found_node)
    !! Returns pointer to first list node with specified tag,
    !! searching in the default forward direction.

    class(list_type), intent(in out) :: self
    character(len = *), intent(in) :: tag
    type(list_node_type), pointer :: found_node

    found_node => self%find(tag, PETSC_FALSE)

  end function list_find_default

!------------------------------------------------------------------------

  function list_get(self, index) result(found_node)
    !! Returns pointer to list node with specified (zero-based) index,
    !! if it exists, or a null pointer otherwise.  index = 0
    !! corresponds to the head node, index = 1 the second node, index
    !! = -1 the tail node, index = -2 the second to last node, etc.

    class(list_type), intent(in out) :: self
    PetscInt, intent(in) :: index
    type(list_node_type), pointer :: found_node
    ! Locals:
    PetscInt :: count, target_count
    PetscBool :: backwards

    found_node => null()
    count = 0
    if (index >= 0) then
       backwards = PETSC_FALSE
       target_count = index
    else
       backwards = PETSC_TRUE
       target_count = -index - 1
    end if

    call self%traverse(get_iterator, backwards)

  contains

    subroutine get_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      if (count == target_count) then
         stopped = PETSC_TRUE
         found_node => node
      else
         stopped = PETSC_FALSE
         count = count + 1
      end if
    end subroutine get_iterator

  end function list_get

!------------------------------------------------------------------------

  subroutine list_destroy_default_forward(self)
    !! Destroys list, in forward order.

    class(list_type), intent(in out) :: self

    do while (associated(self%head))
       call self%delete(self%head)
    end do

  end subroutine list_destroy_default_forward

  subroutine list_destroy_default_reverse(self)
    !! Destroys list, in reverse order.

    class(list_type), intent(in out) :: self

    do while (associated(self%tail))
       call self%delete(self%tail)
    end do

  end subroutine list_destroy_default_reverse

  subroutine list_destroy_default_direction(self, reverse)
    !! Destroys list, with direction specified.

    class(list_type), intent(in out) :: self
    PetscBool, intent(in) :: reverse

    if (reverse) then
       call self%destroy_default_reverse()
    else
       call self%destroy_default_forward()
    end if

  end subroutine list_destroy_default_direction

  subroutine list_destroy_proc_forward(self, node_data_destroy_procedure)
    !! Destroys list, in forward order, applying specified procedure
    !! to destroy the data in the given node.

    class(list_type), intent(in out) :: self
    procedure(list_node_data_destroy_procedure) :: node_data_destroy_procedure

    do while (associated(self%head))
       call self%delete(self%head, node_data_destroy_procedure)
    end do

  end subroutine list_destroy_proc_forward

  subroutine list_destroy_proc_reverse(self, node_data_destroy_procedure)
    !! Destroys list, in reverse order, applying specified procedure
    !! to destroy the data in the given node.

    class(list_type), intent(in out) :: self
    procedure(list_node_data_destroy_procedure) :: node_data_destroy_procedure

    do while (associated(self%tail))
       call self%delete(self%tail, node_data_destroy_procedure)
    end do

  end subroutine list_destroy_proc_reverse

  subroutine list_destroy_proc_direction(self, &
       node_data_destroy_procedure, reverse)
    !! Destroys list, applying specified procedure to destroy the data
    !! in the given node, with direction specified.

    class(list_type), intent(in out) :: self
    procedure(list_node_data_destroy_procedure) :: node_data_destroy_procedure
    PetscBool, intent(in) :: reverse

    if (reverse) then
       call self%destroy_proc_reverse(node_data_destroy_procedure)
    else
       call self%destroy_proc_forward(node_data_destroy_procedure)
    end if

  end subroutine list_destroy_proc_direction

!------------------------------------------------------------------------

  subroutine list_tags(self, tags)
    !! Returns array of tags on each node of the list.

    class(list_type), intent(in out) :: self
    character(*), allocatable, intent(out) :: tags(:)
    ! Locals:
    PetscInt :: i

    allocate(tags(0: self%count - 1))
    i = 0
    call self%traverse(get_tag_iterator)

  contains

    subroutine get_tag_iterator(node, stopped)
       type(list_node_type), pointer, intent(in out)  :: node
       PetscBool, intent(out) :: stopped

       stopped = PETSC_FALSE
       tags(i) = node%tag
       i = i + 1

     end subroutine get_tag_iterator

  end subroutine list_tags

!------------------------------------------------------------------------

  type(list_type) function list_copy(self)
    !! Returns a copy of the list. The copy is assumed not to be the
    !! owner of any data at the list nodes.

    class(list_type), intent(in out) :: self
    ! Locals:

    call list_copy%init(owner = PETSC_FALSE)
    call self%traverse(copy_iterator)

  contains

    subroutine copy_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      if (allocated(node%tag)) then
         call list_copy%append(node%data, node%tag)
      else
         call list_copy%append(node%data)
      end if

    end subroutine copy_iterator

  end function list_copy

!------------------------------------------------------------------------

end module list_module
